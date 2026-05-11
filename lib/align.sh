#!/bin/bash
# lib/align.sh — CLIPittyClip v3.4
# Part of the CLIPittyClip pipeline. Source via lib/modules.sh or directly.
# Auto-split from modules.sh by build_v34_modules.sh.

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ── Demultiplexing ──────────────────────────────────────────────────────────
function run_demultiplexing {
    local input_fastq="$1"
    local barcode_file="$2"
    local run_sample_size="$3"
    local demux_out="${4:-demux_fastq}"
    # Note: dedup_mode parameter ignored - dedup is now handled by caller
    local mismatches="${DEMUX_MISMATCHES:-1}" # Default to 1 if unset/passed global

    log_info "Starting Demultiplexing with cutadapt..."
    log_info "Barcode File: $barcode_file"
    log_info "Allowed Mismatches: $mismatches"
    
    # Input is either the original pooled file or already-deduplicated file
    local work_input="$input_fastq"
    # We call the script relative to the repo root
    local checker_script="${SCRIPT_DIR}/check_barcodes.sh"
    
    if [[ ! -x "$checker_script" ]]; then
        # Fallback if not executable or found (chmod just in case)
        chmod +x "$checker_script" 2>/dev/null
    fi

    "$checker_script" -f "$barcode_file" -m "$mismatches"
    if [ $? -ne 0 ]; then
        log_error "Barcode collisions detected. Aborting pipeline to prevent data mix-up."
        exit 1
    fi
    log_info "Barcode safety check PASSED."

    # 2. Calculate average barcode length to set cutadapt error rate

    # 2. Calculate average barcode length to set cutadapt error rate
    # cutadapt -e is a rate (0.1 = 10%). 
    # rate = mismatches / length
    # We grep the first barcode to estimate length (assuming uniform)
    local first_seq=$(awk '{print $2; exit}' "$barcode_file")
    local bc_len=${#first_seq}
    
    # Calculate rate using awk for floating point
    local error_rate=$(awk "BEGIN {print $mismatches / $bc_len}")
    
    # Cap strictness if 0 errors
    if [[ "$mismatches" -eq 0 ]]; then
        error_rate=0
    fi
    
    log_info "Calculated cutadapt error rate: $error_rate ($mismatches errors in ${bc_len}bp)"

    # Create demux output directory
    mkdir -p "$demux_out"
    
    # Convert barcodes for cutadapt
    local fasta_barcodes="${demux_out}/barcodes.fasta"
    awk '{print ">"$1"\n"$2}' "$barcode_file" > "$fasta_barcodes"

    local cmd="cutadapt \
        -e $error_rate --no-indels \
        -m 1 \
        --action=none \
        -g file:$fasta_barcodes \
        -o \"$demux_out/{name}.fastq\" \
        $work_input \
        -j ${THREADS:-1}"
    
    log_info "Running demultiplexing..."
    execute_cmd "$cmd"

    # Check outputs
    count=$(ls "$demux_out"/*.fastq 2>/dev/null | wc -l)
    if [ "$count" -eq 0 ]; then
        log_error "Demultiplexing failed. No output files created."
        exit 1
    fi
    log_info "Demultiplexing complete. Created $count sample files in '$demux_out/'."

    # Cleanup barcodes FASTA (no longer needed after demux)
    rm -f "$fasta_barcodes"
    
    # Cleanup Deduplicated Temp File if it exists
    if [[ "$work_input" != "$input_fastq" ]] && [[ -f "$work_input" ]]; then
        log_info "Deleting temporary deduplicated file: $work_input"
        rm -f "$work_input"
    fi
}

# run_geo_demux — GEO mode: raw barcode-based splitting, no read modification
# Splits pooled FASTQ by barcode (cutadapt --action=none) and writes gzipped
# per-sample files. No dedup, no fastp. Reads written exactly as received.
# Barcode is found by cutadapt's unanchored 5' search, which tolerates a UMI
# prefix of UMI_LEN bases before the barcode (same behavior as standard demux).
#
# Args: $1 = input_fastq  (.fastq.gz or .fastq)
#       $2 = barcode_file (name<TAB>sequence format)
#       $3 = out_dir      (output directory for split files)
#       $4 = umi_offset   (UMI length; informational only for logging; default: 0)
#       $5 = allowed_mm   (mismatches; default: 1)
run_geo_demux() {
    local input_fastq="$1"
    local barcode_file="$2"
    local out_dir="$3"
    local umi_offset="${4:-0}"
    local allowed_mm="${5:-1}"

    log_info "GEO demux: raw split (no read modification)"
    log_info "UMI offset: ${umi_offset}bp | Mismatches: $allowed_mm"

    # Barcode collision check
    local checker_script="${SCRIPT_DIR}/check_barcodes.sh"
    if [[ -x "$checker_script" ]]; then
        "$checker_script" -f "$barcode_file" -m "$allowed_mm"
        if [[ $? -ne 0 ]]; then
            log_error "Barcode collisions detected. Aborting."
            exit 1
        fi
        log_info "Barcode safety check PASSED."
    fi

    # Calculate cutadapt error rate
    local first_seq
    first_seq=$(awk '!/^#/{print $2; exit}' "$barcode_file")
    local bc_len=${#first_seq}
    local error_rate
    error_rate=$(awk "BEGIN {print $allowed_mm / $bc_len}")
    [[ "$allowed_mm" -eq 0 ]] && error_rate=0
    log_info "Cutadapt error rate: $error_rate ($allowed_mm errors in ${bc_len}bp)"

    mkdir -p "$out_dir"

    # Convert barcodes to FASTA for cutadapt
    local fasta_barcodes="${out_dir}/.barcodes_geo.fasta"
    awk '!/^#/{print ">"$1"\n"$2}' "$barcode_file" > "$fasta_barcodes"

    # Run cutadapt: --action=none preserves reads exactly as received
    # Output is gzipped (.fastq.gz) since these are final GEO deposit files
    local cmd="cutadapt \
        -e $error_rate --no-indels \
        -m 1 \
        --action=none \
        -g file:$fasta_barcodes \
        -o \"${out_dir}/{name}.fastq.gz\" \
        $input_fastq \
        -j ${THREADS:-1}"

    log_info "Running GEO demux..."
    execute_cmd "$cmd"
    local demux_exit=$?
    rm -f "$fasta_barcodes"

    local count
    count=$(ls "${out_dir}"/*.fastq.gz 2>/dev/null | wc -l)
    if [[ $demux_exit -ne 0 || "$count" -eq 0 ]]; then
        log_error "GEO demux failed. No output files created."
        exit 1
    fi
    log_info "GEO demux complete. Created $count files in ${out_dir}/"

    log_info "Calculating MD5 checksums for GEO submission..."
    if command -v md5sum >/dev/null 2>&1; then
        (cd "${out_dir}" && md5sum *.fastq.gz > md5sums.txt)
    elif command -v md5 >/dev/null 2>&1; then
        (cd "${out_dir}" && md5 -r *.fastq.gz > md5sums.txt)
    else
        log_warning "Neither md5sum nor md5 found. Skipping MD5 calculation."
    fi

    # Print summary table
    local total_reads=0
    for f in "${out_dir}"/*.fastq.gz; do
        [[ -f "$f" ]] || continue
        local lines
        lines=$(gzip -dc "$f" | wc -l)
        total_reads=$((total_reads + lines / 4))
    done

    echo ""
    printf "  %-25s %-12s %s\n" "Sample" "Reads" "% of Total"
    echo "  -----------------------------------------------"
    for f in "${out_dir}"/*.fastq.gz; do
        [[ -f "$f" ]] || continue
        local sname
        sname=$(basename "$f" .fastq.gz)
        local lines
        lines=$(gzip -dc "$f" | wc -l)
        local count=$((lines / 4))
        local pct
        pct=$(awk "BEGIN {printf \"%.1f\", ($total_reads > 0) ? ($count / $total_reads) * 100 : 0}")
        printf "  %-25s %-12s %s%%\n" "$sname" "$count" "$pct"
    done
    echo "  -----------------------------------------------"
    echo "  Total: $total_reads reads"
}

# ── Chromosome filter ───────────────────────────────────────────────────────
# Removes contigs like GL000220.1, KI270733.1, etc.
# This reduces data size and improves downstream processing speed
filter_canonical_chromosomes() {
    local input_bam="$1"
    local output_bam="$2"
    
    log_info "Filtering to canonical chromosomes (chr1-22, X, Y, M)..."
    
    # Create list of canonical chromosomes
    local chr_list="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10"
    chr_list+=" chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
    chr_list+=" chr20 chr21 chr22 chrX chrY chrM"
    
    # Count reads before filtering
    local before_count=$(samtools view -c "$input_bam")
    
    # Filter BAM to canonical chromosomes
    samtools view -b "$input_bam" $chr_list > "$output_bam" 2>> "${LOG_FILE}"
    
    if [[ $? -eq 0 && -s "$output_bam" ]]; then
        samtools index "$output_bam"
        local after_count=$(samtools view -c "$output_bam")
        local filtered=$((before_count - after_count))
        log_info "Chromosome filtering complete: $after_count reads kept, $filtered reads on contigs removed"
    else
        log_warning "Chromosome filtering failed. Using unfiltered BAM."
        cp "$input_bam" "$output_bam"
        samtools index "$output_bam"
    fi
}

# ── Deduplication functions moved to lib/dedup.sh ───────────────────────────
# run_dedup(), _fastq_collapse_core(), detect_eclip_umi_length(),
# reformat_eclip_umi_to_sequence(), strip_eclip_barcode()
# are defined in lib/dedup.sh, sourced by CLIPittyClip.sh directly.
# ─────────────────────────────────────────────────────────────────────────────

# ── eCLIP Input Validation ────────────────────────────────────────────────────
# validate_eclip_input — checks that the input FASTQ matches the expected eCLIP mode
# Args: $1 = fastq_file (gzipped or plain)
#       $2 = expected_mode ("pe" or "se")
# Exits 1 with a clear error message on any mismatch.

# ── eCLIP input validation ──────────────────────────────────────────────────
validate_eclip_input() {
    local fastq_file="$1"
    local expected_mode="$2"

    log_info "Validating eCLIP input: $fastq_file (expected mode: $expected_mode)"

    # --- Step 1: Read first 1000 headers ---
    local headers
    if [[ "$fastq_file" == *.gz ]]; then
        headers=$(gzip -cd "$fastq_file" 2>/dev/null | awk 'NR%4==1' | head -1000)
    else
        headers=$(awk 'NR%4==1' "$fastq_file" | head -1000)
    fi

    if [[ -z "$headers" ]]; then
        log_error "eCLIP input validation failed: could not read headers from $fastq_file"
        exit 1
    fi

    # --- Step 2: Detect R1 vs R2 from Illumina comment field ---
    local r1_count r2_count total_count
    r1_count=$(echo "$headers" | awk '{n=split($0,a," "); if(n>1 && a[2]~/^1:N:/) count++} END{print count+0}')
    r2_count=$(echo "$headers" | awk '{n=split($0,a," "); if(n>1 && a[2]~/^2:N:/) count++} END{print count+0}')
    total_count=$(echo "$headers" | wc -l | tr -d ' ')

    local detected_read_end=""
    if [[ "$total_count" -gt 0 ]]; then
        local r1_pct r2_pct
        r1_pct=$(awk "BEGIN {printf \"%.0f\", ($r1_count / $total_count) * 100}")
        r2_pct=$(awk "BEGIN {printf \"%.0f\", ($r2_count / $total_count) * 100}")
        if [[ "$r1_pct" -ge 95 ]]; then
            detected_read_end="1"
        elif [[ "$r2_pct" -ge 95 ]]; then
            detected_read_end="2"
        fi
    fi

    if [[ -z "$detected_read_end" ]]; then
        log_error "eCLIP input validation failed: could not determine R1 or R2 from read headers (checked $total_count headers; ${r1_count} matched 1:N:, ${r2_count} matched 2:N:). Headers must contain a standard Illumina comment field (e.g. 1:N:0:... or 2:N:0:...). Check your FASTQ source."
        exit 1
    fi
    log_info "Detected read end: R${detected_read_end} (from ${total_count} headers)"

    # --- Step 3: Detect UMI location from first header ---
    local first_header
    first_header=$(echo "$headers" | head -1)
    local read_name="${first_header%% *}"   # everything before first space
    local read_name_no_at="${read_name#@}"  # strip leading @
    local token0="${read_name_no_at%%:*}"   # token before first colon

    local umi_location="sequence"
    local detected_umi_len=0
    if [[ "${#token0}" -ge 5 && "${#token0}" -le 10 ]] && [[ "$token0" =~ ^[ACGTN]+$ ]]; then
        umi_location="header_colon"
        detected_umi_len="${#token0}"
    elif [[ "$read_name_no_at" =~ _[ACGTN]{5,10}$ ]]; then
        umi_location="header_underscore"
    fi
    log_info "Detected UMI location: $umi_location"

    # --- Step 4: Validate against expected mode ---
    if [[ "$expected_mode" == "pe" ]]; then
        if [[ "$detected_read_end" == "1" ]]; then
            log_error "Input validation failed for --eclip pe: detected Read 1 input. PE eCLIP analysis requires Read 2, which contains the cross-link site at its 5' end. Please supply the corresponding R2 fastq file. See README section eCLIP-PE."
            exit 1
        fi
        # detected_read_end == "2" from here
        if [[ "$umi_location" == "sequence" ]]; then
            log_error "Input validation failed for --eclip pe: detected R2 read with UMI in sequence (raw/pre-eclipdemux format). CLIPittyClip --eclip pe requires post-eclipdemux files where the UMI has been moved to the read header by eclipdemux. Please obtain the post-eclipdemux R2 fastq from ENCODE. See README section eCLIP-PE."
            exit 1
        elif [[ "$umi_location" == "header_underscore" ]]; then
            log_error "Input validation failed for --eclip pe: detected umi_tools-style UMI in read header. CLIPittyClip --eclip pe expects post-eclipdemux format (UMI as colon-prefixed token, e.g. @NTACGTTGAT:...). Please use the post-eclipdemux file from ENCODE."
            exit 1
        fi
        # PASS: R2, header_colon
        log_info "Validation passed: PE eCLIP R2, post-eclipdemux format. UMI ${detected_umi_len}nt detected in read header."
        return 0

    elif [[ "$expected_mode" == "se" ]]; then
        if [[ "$detected_read_end" == "2" ]]; then
            log_error "Input validation failed for --eclip se: detected Read 2 input. --eclip se expects Read 1 from single-end eCLIP (seCLIP) sequencing. If you have paired-end eCLIP data, use --eclip pe with the R2 file instead."
            exit 1
        fi
        # detected_read_end == "1" from here
        if [[ "$umi_location" == "header_colon" ]]; then
            log_error "Input validation failed for --eclip se: detected UMI in read header (eclipdemux-style). --eclip se expects raw seCLIP fastq where the UMI is still in the read sequence (first 10nt). Please supply the unprocessed fastq from ENCODE."
            exit 1
        elif [[ "$umi_location" == "header_underscore" ]]; then
            log_error "Input validation failed for --eclip se: detected umi_tools-style UMI in read header. --eclip se expects raw seCLIP fastq where the UMI is still in the read sequence. Please supply the unprocessed fastq."
            exit 1
        fi
        # PASS: R1, sequence
        log_info "Validation passed: SE eCLIP R1, raw format. UMI assumed 10nt in sequence (seCLIP standard, Blue et al. 2022)."
        return 0
    fi
}

# ── eCLIP preprocessing ─────────────────────────────────────────────────────
# ── eCLIP PE Preprocessing ───────────────────────────────────────────────────
# run_eclip_pe_preprocessing — full PE eCLIP preprocessing chain
# Expected input: post-eclipdemux R2 fastq (UMI in read header as colon-prefix token)
# Flow: validate → UMI to seq → Deduplicate → Extract UMI → Adapter Trim
# Args: $1 = input_file, $2 = output_prefix, $3 = threads, $4 = sample_size, $5 = umi_len (hint)
run_eclip_pe_preprocessing() {
    local input_file="$1"
    local output_prefix="$2"
    local threads="$3"
    local sample_size="$4"
    local umi_len="${5:-0}"

    # Get path to eCLIP adapter FASTA (same dir as this script)
    local script_dir="$(dirname "${BASH_SOURCE[0]}")"
    local eclip_adapters_fasta="${script_dir}/eclip_adapters.fa"

    log_info "eCLIP PE mode: Expecting Read 2, post-eclipdemux format (UMI in read header)."
    validate_eclip_input "$input_file" "pe"

    log_info "eCLIP PE mode: Preprocessing workflow (validate → UMI to seq → Deduplicate → Extract UMI → Adapter Trim)"

    # Step 1: Detect UMI length from header
    local detected_umi_len
    detected_umi_len=$(detect_eclip_umi_length "$input_file" "$umi_len")
    umi_len="$detected_umi_len"

    # Step 2: Move UMI from header to sequence (required before collapse)
    update_status_first "eCLIP PE Preprocessing"
    echo -ne "(UMI to Sequence"
    local umi_seq_file="${output_prefix}_umi_in_seq.fastq.gz"
    reformat_eclip_umi_to_sequence "$input_file" "$umi_seq_file" "$umi_len"
    if [[ ! -s "$umi_seq_file" ]]; then
        log_error "Failed to reformat eCLIP UMI to sequence"
        exit 1
    fi

    # Step 3: Deduplicate (collapse exact duplicates with UMI in sequence)
    echo -ne " → Deduplicating"
    local collapsed_file="${output_prefix}_collapsed.fastq"
    local collapsed_file_gz="${output_prefix}_collapsed.fastq.gz"
    gzip -dc "$umi_seq_file" > "${output_prefix}_umi_temp.fastq"
    _fastq_collapse_core "${output_prefix}_umi_temp.fastq" "$collapsed_file"
    if [[ ! -s "$collapsed_file" ]]; then
        log_error "Deduplication (eCLIP PE collapse) failed"
        exit 1
    fi
    gzip -c "$collapsed_file" > "$collapsed_file_gz"
    rm -f "${output_prefix}_umi_temp.fastq" "$collapsed_file" "$umi_seq_file"

    # Step 4: Strip UMI from sequence, attach to header after count (CTK format: READ#count#UMI)
    echo -ne " → Extract UMI"
    local stripped_file="${output_prefix}_stripped.fastq.gz"
    strip_eclip_barcode "$collapsed_file_gz" "$stripped_file" "$umi_len"
    rm -f "$collapsed_file_gz"
    if [[ ! -s "$stripped_file" ]]; then
        log_error "stripBarcode.pl failed"
        exit 1
    fi

    # Step 5: Adapter trimming with fastp (using all eCLIP inline-barcode + TruSeq R2 adapters)
    echo -ne " → Adapter Trim) > "
    local final_file="${output_prefix}_cleaned.fastq"
    local fastp_cmd="fastp -i ${stripped_file} -o ${final_file} \
        --thread ${threads} \
        --adapter_fasta ${eclip_adapters_fasta} \
        --length_required 20 \
        --cut_tail --cut_tail_mean_quality 5 \
        --overlap_len_require 1 \
        --html ${output_prefix}_fastp.html \
        --json ${output_prefix}_fastp.json"
    if [ "$sample_size" -gt 0 ]; then fastp_cmd+=" --reads_to_process $sample_size"; fi
    log_info "Running: $fastp_cmd"
    execute_cmd "$fastp_cmd"
    rm -f "$stripped_file"

    if [[ ! -s "$final_file" ]]; then
        log_error "fastp failed to create cleaned file"
        exit 1
    fi

    log_info "eCLIP PE preprocessing complete: $final_file"
    log_info "Read ID format: READ#count#UMI (CTK-compatible)"

    # Export detected UMI length for downstream tag2collapse.pl
    ECLIP_UMI_LEN="$umi_len"
}

# ── eCLIP SE Preprocessing ───────────────────────────────────────────────────
# run_eclip_se_preprocessing — full SE eCLIP (seCLIP) preprocessing chain
# Expected input: raw Read 1 fastq (UMI in first 10nt of sequence, Blue et al. 2022)
# Flow: validate → Deduplicate → Extract UMI → Adapter Trim
# Args: $1 = input_file, $2 = output_prefix, $3 = threads, $4 = sample_size
run_eclip_se_preprocessing() {
    local input_file="$1"
    local output_prefix="$2"
    local threads="$3"
    local sample_size="$4"

    # Hardcoded SE parameters (Blue et al. 2022 — not user-configurable)
    local umi_len=10
    local adapter_seq="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # TruSeq Read 1 adapter

    log_info "eCLIP SE mode: Expecting raw Read 1, seCLIP format (UMI in read sequence)."
    log_info "SE eCLIP UMI length: ${umi_len}nt (seCLIP standard, Blue et al. 2022)"
    log_info "SE eCLIP adapter: TruSeq R1 (${adapter_seq})"
    validate_eclip_input "$input_file" "se"

    log_info "eCLIP SE mode: Preprocessing workflow (validate → Deduplicate → Extract UMI → Adapter Trim)"

    update_status_first "eCLIP SE Preprocessing"

    # Step 1: Deduplicate (UMI is in sequence — collapse on full read including UMI prefix)
    echo -ne "(Deduplicating"
    local temp_fastq="${output_prefix}_se_temp.fastq"
    local collapsed_plain="${output_prefix}_collapsed.fastq"
    local collapsed_gz="${output_prefix}_collapsed.fastq.gz"
    if [[ "$input_file" == *.gz ]]; then
        gzip -dc "$input_file" > "$temp_fastq"
    else
        cp "$input_file" "$temp_fastq"
    fi
    _fastq_collapse_core "$temp_fastq" "$collapsed_plain"
    if [[ ! -s "$collapsed_plain" ]]; then
        log_error "Deduplication (eCLIP SE collapse) failed"
        exit 1
    fi
    gzip -c "$collapsed_plain" > "$collapsed_gz"
    rm -f "$temp_fastq" "$collapsed_plain"

    # Step 2: Strip UMI from sequence, attach to header after count (CTK format: READ#count#UMI)
    echo -ne " → Extract UMI"
    local stripped_gz="${output_prefix}_stripped.fastq.gz"
    strip_eclip_barcode "$collapsed_gz" "$stripped_gz" "$umi_len"
    rm -f "$collapsed_gz"
    if [[ ! -s "$stripped_gz" ]]; then
        log_error "stripBarcode.pl failed"
        exit 1
    fi

    # Step 3: Adapter trimming with fastp (TruSeq R1 adapter, passed as sequence string)
    echo -ne " → Adapter Trim) > "
    local final_file="${output_prefix}_cleaned.fastq"
    local fastp_cmd="fastp -i ${stripped_gz} -o ${final_file} \
        --thread ${threads} \
        --adapter_sequence ${adapter_seq} \
        --length_required 20 \
        --cut_tail --cut_tail_mean_quality 5 \
        --overlap_len_require 1 \
        --html ${output_prefix}_fastp.html \
        --json ${output_prefix}_fastp.json"
    if [ "$sample_size" -gt 0 ]; then fastp_cmd+=" --reads_to_process $sample_size"; fi
    log_info "Running: $fastp_cmd"
    execute_cmd "$fastp_cmd"
    rm -f "$stripped_gz"

    if [[ ! -s "$final_file" ]]; then
        log_error "fastp failed to create cleaned file"
        exit 1
    fi

    log_info "eCLIP SE preprocessing complete: $final_file"
    log_info "Read ID format: READ#count#UMI (CTK-compatible)"

    # Export detected UMI length for downstream tag2collapse.pl
    ECLIP_UMI_LEN="$umi_len"
}

# ── Adapter trimming ────────────────────────────────────────────────────────
# 1b. Adapter trimming and quality filtering with fastp (standard mode only)
run_fastp() {
    local input_file="$1"
    local output_prefix="$2"
    local umi_len="$3"
    local adapter3="$4"
    local threads="$5"
    local sample_size="$6"
    local bc_len="${7:-0}"
    local spacer_len="${8:-0}"
    local bc_first="${9:-false}"   # --bc-first: layout is [BC][UMI][sp][READ] not [UMI][BC][sp][READ]

    local cleaned="${output_prefix}_cleaned.fastq"

    update_status_first "Adapter Trimming"

    # ── Common fastp quality / length flags ───────────────────────────────────
    local qc_flags="--thread ${threads} --length_required 16 --average_qual 30"
    local sample_flag=""
    if [ "$sample_size" -gt 0 ]; then sample_flag="--reads_to_process $sample_size"; fi
    local adapter_flag=""
    if [ -n "$adapter3" ]; then adapter_flag="--adapter_sequence ${adapter3}"; fi

    if [[ "$bc_first" == "true" ]]; then
        # ── BC-first layout: [BC][UMI][spacer][READ] ─────────────────────────
        # fastp --umi_loc=read1 always extracts from position 0, so if we ran
        # UMI extraction on the raw read we would grab part of the barcode.
        # Solution: two fastp passes.
        #   Pass 1 — strip BC only  →  [UMI][spacer][READ]
        #   Pass 2 — extract UMI, strip spacer, trim adapter  →  [READ]
        log_info "BC-first mode: [BC($bc_len)][UMI($umi_len)][sp($spacer_len)][READ]"
        log_info "  Pass 1: trim ${bc_len}nt barcode from 5' end"

        local tmp="${output_prefix}_bcfirst_tmp.fastq"

        local pass1="fastp -i ${input_file} -o ${tmp} \
            ${qc_flags} \
            --disable_adapter_trimming \
            --html /dev/null --json /dev/null"
        if [ "$bc_len" -gt 0 ]; then pass1+=" --trim_front1 ${bc_len}"; fi
        if [ "$sample_size" -gt 0 ]; then pass1+=" $sample_flag"; fi

        log_info "Running (pass 1): $pass1"
        execute_cmd "$pass1"
        if [ $? -ne 0 ] || [ ! -s "$tmp" ]; then
            log_error "fastp BC-first pass 1 failed."
            rm -f "$tmp"; exit 1
        fi

        log_info "  Pass 2: extract UMI(${umi_len}nt), trim spacer(${spacer_len}nt), trim adapter"
        local pass2="fastp -i ${tmp} -o ${cleaned} \
            ${qc_flags} \
            --html ${output_prefix}_fastp.html \
            --json ${output_prefix}_fastp.json \
            ${adapter_flag}"
        if [ "$umi_len" -gt 0 ]; then pass2+=" --umi --umi_loc=read1 --umi_len=${umi_len} --umi_delim=#"; fi
        if [ "$spacer_len" -gt 0 ]; then pass2+=" --trim_front1 ${spacer_len}"; fi

        log_info "Running (pass 2): $pass2"
        execute_cmd "$pass2"
        local exit_code=$?
        rm -f "$tmp"

    else
        # ── UMI-first layout (default): [UMI][BC][spacer][READ] ──────────────
        log_info "Standard mode: fastp adapter trimming"

        local fastp_cmd="fastp -i ${input_file} -o ${cleaned} \
            ${qc_flags} \
            --html ${output_prefix}_fastp.html \
            --json ${output_prefix}_fastp.json \
            ${adapter_flag}"
        if [ "$umi_len" -gt 0 ]; then fastp_cmd+=" --umi --umi_loc=read1 --umi_len=${umi_len} --umi_delim=#"; fi
        local front_trim=$(( bc_len + spacer_len ))
        if [ "$front_trim" -gt 0 ]; then fastp_cmd+=" --trim_front1 ${front_trim}"; fi
        if [ "$sample_size" -gt 0 ]; then fastp_cmd+=" $sample_flag"; fi

        log_info "Running: $fastp_cmd"
        execute_cmd "$fastp_cmd"
        local exit_code=$?
    fi

    if [ $exit_code -eq 0 ]; then
        log_info "Adapter trimming complete."
    else
        log_error "fastp failed."
        exit 1
    fi
}

# ── ncRNA pre-filter ────────────────────────────────────────────────────────
# 1c. ncRNA Pre-filtering with Bowtie2

# Filters out rRNA, tRNA, and other ncRNA reads before genome alignment
# Input: FASTQ from fastp
# Output: Unmapped reads (for genome alignment), Mapped reads (QC)
run_ncrna_filter() {
    local input_fastq="$1"
    local output_unmapped="$2"    # Reads that didn't map to ncRNA (continue to genome)
    local output_dir="$3"         # Directory for ncRNA mapping outputs
    local index_dir="$4"
    local threads="$5"
    local sample_name="$6"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    local ncrna_bam="${output_dir}/${sample_name}_ncrna.bam"
    local ncrna_stats="${output_dir}/${sample_name}_ncrna_stats.txt"
    
    update_status "ncRNA Filter"
    
    # Run Bowtie2 mapping to ncRNA index
    # --un: write unmapped reads to plain .fastq (these go to genome mapping)
    # Mapped reads are saved to BAM for QC
    local bt2_cmd="bowtie2 -x \"${index_dir}/ncrna\" \
        -U \"$input_fastq\" \
        --un \"$output_unmapped\" \
        -p $threads \
        2> \"$ncrna_stats\" \
        | samtools view -bS - > \"$ncrna_bam\""
    
    log_info "Running ncRNA filter: $bt2_cmd"
    
    if execute_cmd "$bt2_cmd"; then
        # Index the BAM for potential downstream use
        samtools index "$ncrna_bam" 2>/dev/null || true
        
        # Extract alignment rate from stats and display on console
        local align_rate=$(grep "overall alignment rate" "$ncrna_stats" | grep -oE "[0-9]+\.[0-9]+%" || echo "N/A")
        local total_reads=$(grep "reads; of these:" "$ncrna_stats" | grep -oE "^[0-9]+" || echo "N/A")
        local aligned_reads=$(grep "aligned exactly 1 time" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        local multi_aligned=$(grep "aligned >1 times" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        local ncrna_reads=$((aligned_reads + multi_aligned))
        
        log_info "ncRNA alignment rate: $align_rate"
        # Note: Per-sample stats logged to file; summary table shown after batch
        
        return 0
    else
        log_error "ncRNA filtering failed"
        return 1
    fi
}

# ── STAR alignment ──────────────────────────────────────────────────────────
# 2. Mapping with STAR
run_mapping_star() {
    local input_fastq="$1"
    local output_prefix="$2"
    local genome_dir="$3"
    local threads="$4"
    local mismatch_max="$5"

    # Resolve to absolute path — STAR's --readFilesCommand spawns cat/gzip in its
    # own temp dir, so relative paths in --readFilesIn will not be found.
    input_fastq="$(cd "$(dirname "$input_fastq")" && pwd)/$(basename "$input_fastq")"

    update_status "Mapping (STAR)"
    log_info "Starting mapping with STAR..."

    # For .fastq.gz: use --readFilesCommand gzip -dc (STAR spawns a decompressor).
    # For plain .fastq: omit --readFilesCommand entirely so STAR reads the file
    # directly without spawning a subprocess. Using --readFilesCommand cat forces
    # STAR into FIFO mode which fails on macOS even with absolute paths.
    local reads_command_flag=""
    if [[ "$input_fastq" == *.gz ]]; then
        reads_command_flag="--readFilesCommand gzip -dc"
    fi

    # Create temp directory for STAR (prevents FIFO errors on exFAT/NTFS drives)
    local star_tmp="${TMPDIR:-/tmp}/star_$(basename "$output_prefix")_$$"
    mkdir -p "$star_tmp"
    log_info "STAR temp directory: $star_tmp"

    local cmd="STAR --runThreadN ${threads} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${input_fastq} \
        ${reads_command_flag} \
        --outFileNamePrefix ${output_prefix}. \
        --outTmpDir ${star_tmp}/STARtmp \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverReadLmax 0.1 \
        --outFilterMismatchNmax 5 \
        --outSAMattributes NH HI AS nM NM MD \
        --alignEndsType EndToEnd \
        --scoreDelOpen -1 --scoreDelBase -1 \
        --scoreInsOpen -1 --scoreInsBase -1 \
        $ADV_ALIGNER_ARGS"
    # Mismatch filtering:
    #   --outFilterMismatchNoverReadLmax 0.1: fractional filter (allows ~3 mismatches in 30bp,
    #   ~2 in 20bp). Scales with post-trim read length like BWA's -n 0.06 philosophy.
    #   This is the primary filter; replaces the old absolute --outFilterMismatchNmax 2 which
    #   would discard reads carrying a genuine crosslink deletion + 2 sequencing errors (NM=3).
    #   --outFilterMismatchNmax 5: hard backstop only — blocks extreme cases.
    # Gap penalties:
    #   --scoreDelOpen/Base -1: lower deletion penalty to match mismatch cost.
    #   --scoreInsOpen/Base -1: lower insertion penalty symmetrically.
    #   STAR's default gap penalties (-2 open, -2 base) make 1-nt indels score worse than
    #   substitutions on short reads (~20-30bp post-trim iCLIP), causing CIMS-relevant
    #   deletion-containing reads to be suppressed or realigned as mismatches.

    log_info "Running: $cmd"
    execute_cmd "$cmd"
    local star_exit=$?
    
    # Cleanup temp directory (always, even on failure)
    rm -rf "$star_tmp"
    
    if [ $star_exit -ne 0 ]; then
        log_error "STAR mapping failed. Check the log file for details."
        exit 1
    fi
    
    samtools index "${output_prefix}.Aligned.sortedByCoord.out.bam"
    if [ $? -ne 0 ]; then
        log_error "samtools index failed. The BAM file might be empty or invalid."
        exit 1
    fi

    log_info "Mapping complete. Output: ${output_prefix}.Aligned.sortedByCoord.out.bam"
}

# ── Bowtie2 alignment ───────────────────────────────────────────────────────
# 2b. Mapping with Bowtie2
run_mapping_bowtie2() {
    local input_file="$1"
    local output_prefix="$2"
    local genome_index="$3"
    local threads="$4"
    
    update_status "Mapping (Bowtie2)"
    log_info "Starting mapping with Bowtie2..."
    log_info "Input: $input_file"
    log_info "Index: $genome_index"
    
    # Verify index (look for .1.bt2 or .1.bt2l, excluding ncRNA patterns)
    # First try top-level (maxdepth 1), then fall back to deeper search
    local found_idx=""
    
    # Try .1.bt2 at top level first, excluding ncRNA
    found_idx=$(find "$genome_index" -maxdepth 1 -name "*.1.bt2" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    
    # If not found, try .1.bt2l at top level
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -maxdepth 1 -name "*.1.bt2l" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    
    # Fall back to deeper search (excluding ncRNA subfolder)
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -name "*.1.bt2" ! -path "*/ncRNA/*" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -name "*.1.bt2l" ! -path "*/ncRNA/*" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    
    if [[ -z "$found_idx" ]]; then
        log_error "Bowtie2 index files (*.1.bt2) not found in $genome_index"
        return 1
    fi
    
    # Construct base name for index
    # Standard: /path/hg38.1.bt2 -> Base: /path/hg38
    # We strip the .1.bt2 suffix
    local idx_base="${found_idx%.1.bt2}"
    if [[ "$idx_base" == "$found_idx" ]]; then
         idx_base="${found_idx%.1.bt2l}"
    fi

    local sam_file="${output_prefix}.sam"
    local bam_file="${output_prefix}.Aligned.sortedByCoord.out.bam" # Match STAR naming for compatibility

    # CIMS-optimized Bowtie2 parameters (BWA aln -n 0.06 equivalent):
    #
    #   --end-to-end          No soft clipping; BWA aln is also end-to-end
    #   --mp 2,2              Uniform mismatch penalty (removes quality weighting);
    #                         BWA counts edits, not quality-weighted costs
    #   --rdg 1,1             Read gap open=1, extend=1 → 1-base deletion costs 2 pts
    #   --rfg 1,1             Ref gap (insertion) same → equal to 1 mismatch at --mp 2,2
    #                         This gives deletion/mismatch parity, matching BWA's edit budget
    #   --score-min L,0,-0.12 Linear floor: 0.06 edits/bp × 2 pts/edit = 0.12 pts/bp
    #                         Mirrors BWA -n 0.06 fractional edit distance, scaled by read length
    #                         e.g. 36bp → -4.3 (≈2 edits), 50bp → -6.0 (≈3 edits)
    #   -N 1                  Allow 1 mismatch in seed (BWA backtracking is natively more sensitive)
    #   -L 16                 Short seed for CLIP read lengths (25-50bp); more anchor positions
    #   -k 10                 Report up to 10 alignments (consistent with STAR multimapper cap)
    #   --no-unal             Suppress unaligned reads from SAM output
    local cmd="bowtie2 -p $threads \
  --end-to-end \
  --mp 2,2 \
  --rdg 1,1 \
  --rfg 1,1 \
  --score-min L,0,-0.12 \
  -N 1 \
  -L 16 \
  -k 10 \
  --no-unal \
  -x '$idx_base' -U '$input_file' -S '$sam_file' $ADV_ALIGNER_ARGS"

    log_info "Running Bowtie2 (CIMS-tuned, BWA aln -n 0.06 equivalent)..."
    execute_cmd "$cmd"
    
    if [ $? -ne 0 ]; then
        log_error "Bowtie2 alignment failed."
        return 1
    fi
    
    # Convert to BAM -> Sort -> Index
    # Note: "Processing Alignment" status is now in run_parse_alignment
    
    local sort_cmd="samtools view -bS '$sam_file' | samtools sort -@ $threads -o '$bam_file' -"
    execute_cmd "$sort_cmd"
    
    if [ -f "$bam_file" ]; then
        samtools index "$bam_file"
        rm -f "$sam_file" # Cleanup SAM
    else
        log_error "BAM conversion failed."
        return 1
    fi
    log_info "Mapping complete. Output: $bam_file"
}

# ── CTK alignment parser ────────────────────────────────────────────────────

# 2c. Parse Alignment for CIMS/CITS
# parseAlignment.pl requires SAM input (not BAM)
# This function converts BAM→SAM, optionally runs calmd for MD tags, then parses
run_parse_alignment() {
    local bam_file="$1"
    local output_bed="$2"
    local mutation_file="$3"
    local genome_index="$4"
    
    update_status "Processing Alignment"
    log_info "Parsing alignment for CIMS/CITS..."
    log_info "Input BAM: $bam_file"
    
    # Step 1: Convert BAM to SAM (parseAlignment.pl requires SAM input)
    local sam_file="${bam_file%.bam}.sam"
    log_info "Converting BAM to SAM for parseAlignment.pl..."
    samtools view -h "$bam_file" > "$sam_file"
    
    if [[ ! -s "$sam_file" ]]; then
        log_error "BAM to SAM conversion failed or empty output: $sam_file"
        return 1
    fi
    
    local work_sam="$sam_file"

    # Step 2: (Optional) Run samtools calmd for MD tag standardization
    # calmd recalculates MD tags from the reference FASTA authoritatively.
    # parseAlignment.pl uses both CIGAR and MD to classify mutations; inconsistent MD
    # (especially at deletion boundaries in homopolymer runs) = missed/misclassified deletions.
    # Priority: --genome-fasta flag > find in index dir (STAR index dirs never contain FASTA).
    local ref_fasta=""
    if [[ -n "${GENOME_FASTA:-}" ]] && [[ -f "$GENOME_FASTA" ]]; then
        ref_fasta="$GENOME_FASTA"
        log_info "Using provided genome FASTA for calmd: $ref_fasta"
    else
        # Fallback: search genome index directory (works for Bowtie2, not STAR)
        ref_fasta=$(find "$genome_index" -maxdepth 2 -name "*.fa" -o -name "*.fasta" 2>/dev/null | grep -v '^\._' | head -n 1)
        if [[ -n "$ref_fasta" ]]; then
            log_info "Reference FASTA found in index dir: $ref_fasta"
        fi
    fi

    if [[ -n "$ref_fasta" ]]; then
        log_info "Running 'samtools calmd' to standardize MD tags..."

        local calmd_sam="${sam_file%.sam}_calmd.sam"
        # calmd can take BAM input and output SAM (-S flag forces SAM output)
        samtools calmd -S "$bam_file" "$ref_fasta" > "$calmd_sam" 2>/dev/null

        if [[ -s "$calmd_sam" ]]; then
            work_sam="$calmd_sam"
            log_info "Using calmd-processed SAM: $calmd_sam"
        else
            log_warning "samtools calmd failed or empty. Using native tags."
            rm -f "$calmd_sam" 2>/dev/null
        fi
    else
        if [[ "${RUN_CIMS:-false}" == "true" ]]; then
            log_warning "CIMS WARNING: No reference FASTA available for samtools calmd."
            log_warning "  MD tags from STAR may be inconsistent at deletion boundaries,"
            log_warning "  which can cause parseAlignment.pl to miss or misclassify crosslink"
            log_warning "  deletions. Provide --genome-fasta for optimal CIMS detection."
        else
            log_info "No reference FASTA found. Relying on aligner's native MD tags."
        fi
    fi

    # Step 3: Run parseAlignment.pl
    # Options:
    #   -v: verbose
    #   --map-qual 1: Require unique mapping (MAPQ >= 1)
    #   --min-len 16: Minimum read length
    #   --mutation-file: Output mutation file for CIMS
    log_info "Running parseAlignment.pl..."
    local cmd="parseAlignment.pl -v --map-qual 1 --min-len 16 --mutation-file '$mutation_file' '$work_sam' '$output_bed'"
    
    execute_cmd "$cmd"
    local parse_exit=$?
    
    # Step 4: Cleanup temp SAM files
    rm -f "$sam_file" 2>/dev/null
    if [[ -f "${sam_file%.sam}_calmd.sam" ]]; then
        rm -f "${sam_file%.sam}_calmd.sam" 2>/dev/null
    fi
    
    if [[ $parse_exit -ne 0 ]]; then
        log_error "parseAlignment.pl failed with exit code $parse_exit"
        return $parse_exit
    fi
    
    # Verify outputs
    if [[ ! -s "$output_bed" ]]; then
        log_warning "parseAlignment.pl produced empty BED file: $output_bed"
    else
        local bed_count=$(wc -l < "$output_bed")
        log_info "parseAlignment.pl output: $bed_count tags in $output_bed"
    fi
    
    if [[ ! -s "$mutation_file" ]]; then
        log_warning "parseAlignment.pl produced empty mutation file: $mutation_file"
    else
        local mut_count=$(wc -l < "$mutation_file")
        log_info "parseAlignment.pl output: $mut_count mutations in $mutation_file"
    fi
    
    log_info "Alignment parsing complete."
}
