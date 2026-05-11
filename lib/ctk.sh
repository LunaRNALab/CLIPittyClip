#!/bin/bash
# lib/ctk.sh — CLIPittyClip v3.4
# Part of the CLIPittyClip pipeline. Source via lib/modules.sh or directly.
# Auto-split from modules.sh by build_v34_modules.sh.

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ── PCR duplicate collapse (tag2collapse) ───────────────────────────────────
# 3. PCR Duplicate Removal (using -c for chromosome-based processing to prevent OOM)
run_collapse_pcr() {
    local input_bed="$1"
    local output_bed="$2"
    local umi_len="${3:-0}"  # Optional: UMI length, defaults to 0
    local dedup_mode="${4:-true}"  # Was fastq2collapse.pl run? If not, no count in read names

    update_status "Collapsing"
    log_info "Collapsing PCR duplicates with CTK tag2collapse.pl..."

    # Only use --random-barcode when UMI is present in read names
    local barcode_flag=""
    local weight_flags=""
    if [ "${umi_len:-0}" -gt 0 ]; then
        barcode_flag="--random-barcode"
        if [[ "$dedup_mode" == "true" ]]; then
            # fastq2collapse.pl ran first: read names are READ#count#UMI
            # --weight: use tag count as weight
            # --weight-in-name: read count from #count# in read ID
            # -EM 30: EM iterations for barcode collapse
            # --seq-error-model alignment: estimate error from alignment
            weight_flags="--weight --weight-in-name -EM 30 --seq-error-model alignment"
            log_info "UMI mode (length=$umi_len): using random barcode collapse with EM weighting"
        else
            # No fastq2collapse.pl: read names are READ#UMI only, no count embedded
            # Use EM barcode collapse without weight-in-name
            weight_flags="-EM 30 --seq-error-model alignment"
            log_info "UMI mode (length=$umi_len): using random barcode collapse (no pre-collapse counts)"
        fi
    else
        log_info "No UMI: using position-only collapse"
    fi
    
    # Create temp cache directory path for -c option (tag2collapse.pl creates it)
    local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/collapse_cache.XXXXXX")
    
    # Use -big (memory-mapped BIG format) and -c (chromosome-based) for max memory efficiency
    local cmd="$CONDA_PREFIX/bin/perl $(which tag2collapse.pl) -big -c \"${cache_dir}\" --keep-tag-name --keep-max-score ${barcode_flag} ${weight_flags} \
        \"${input_bed}\" \"${output_bed}\""

    execute_cmd "$cmd"
    local exit_code=$?
    
    # Cleanup cache directory
    rm -rf "$cache_dir"

    if [ $exit_code -eq 0 ] && [[ -s "$output_bed" ]]; then
        local output_count=$(wc -l < "$output_bed")
        log_info "Collapsing complete. Output: $output_count tags"
    else
        log_error "PCR duplicate removal failed."
        exit 1
    fi
}

# ── CTK preprocessing ───────────────────────────────────────────────────────
# CTK Preprocessing: Filter mutations and extract mutation types
# Called after tag2collapse.pl, before CIMS/CITS
run_ctk_preprocessing() {
    local collapsed_bed="$1"
    local raw_mutation_file="$2"
    local output_dir="$3"
    
    # Status update removed - preprocessing is already done in main pipeline
    log_info "Preprocessing mutations for CTK analysis..."
    
    mkdir -p "$output_dir"
    
    # Step 1: Filter mutations to only those in collapsed tags
    # selectRow.pl uses zero-based column indexing: column 3 = read name
    log_info "Filtering mutations to collapsed tags (selectRow.pl -q 3 -f 3)..."
    local matched_file="${output_dir}/mutations_matched.txt"
    
    selectRow.pl -q 3 -f 3 "$raw_mutation_file" "$collapsed_bed" > "$matched_file"
    
    if [[ ! -s "$matched_file" ]]; then
        log_warning "No matching mutations found after filtering."
        return 1
    fi
    
    local matched_count=$(wc -l < "$matched_file")
    log_info "Matched mutations: $matched_count"
    
    # Step 2: Extract mutation types using getMutationType.pl
    log_info "Extracting deletion mutations..."
    local del_file="${output_dir}/deletions.bed"
    getMutationType.pl -t del "$matched_file" "$del_file"
    
    log_info "Extracting substitution mutations..."
    local sub_file="${output_dir}/substitutions.bed"
    getMutationType.pl -t sub "$matched_file" "$sub_file"
    
    # Report counts
    if [[ -s "$del_file" ]]; then
        local del_count=$(wc -l < "$del_file")
        log_info "Deletions extracted: $del_count"
    else
        log_warning "No deletions found."
    fi
    
    if [[ -s "$sub_file" ]]; then
        local sub_count=$(wc -l < "$sub_file")
        log_info "Substitutions extracted: $sub_count"
    else
        log_warning "No substitutions found."
    fi
    
    log_info "CTK preprocessing complete. Output: $output_dir"
}

# ── CIMS analysis ───────────────────────────────────────────────────────────
# 5. CIMS Analysis (CTK) - with parallel chromosome processing
# Detects crosslinking-induced mutation sites
# Input: collapsed BED + mutation file (deletions or substitutions)
# Output: CIMS.txt with significant mutation sites
run_cims() {
    local input_collapsed_bed="$1"
    local mutation_bed="$2"          # Already BED6 from getMutationType.pl
    local output_file="$3"
    local cims_iterations="${4:-10}"  # Default: 10 iterations
    local cims_fdr="${5:-0.001}"      # Default: FDR 0.001
    local threads="${THREADS:-4}"     # Use global THREADS or default
    
    log_info "Running CIMS analysis..."
    log_info "Input tags: $input_collapsed_bed"
    log_info "Input mutations: $mutation_bed"
    log_info "Iterations: $cims_iterations, FDR threshold: $cims_fdr"
    
    # Verify inputs exist
    if [[ ! -s "$input_collapsed_bed" ]]; then
        log_error "CIMS: Collapsed BED file is empty or missing: $input_collapsed_bed"
        return 1
    fi
    if [[ ! -s "$mutation_bed" ]]; then
        log_error "CIMS: Mutation BED file is empty or missing: $mutation_bed"
        return 1
    fi
    
    # Set up CTK environment
    export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
    
    # Check if parallel processing is available and beneficial
    local use_parallel="false"
    local chr_count=$(cut -f1 "$input_collapsed_bed" | sort -u | wc -l)
    
    if has_gnu_parallel && [[ $chr_count -gt 1 ]]; then
        use_parallel="true"
        log_info "Using parallel processing ($chr_count chromosomes, $threads threads)"
    fi
    
    if [[ "$use_parallel" == "true" ]]; then
        # Parallel mode: split BOTH collapsed.bed AND mutation file by chromosome
        local chunk_dir=$(mktemp -d "${TMPDIR:-/tmp}/cims_parallel.XXXXXX")
        
        # Split collapsed BED by chromosome (standard chromosomes only: chr1-22, X, Y, M)
        # This filters out contigs like KI270737.1, GL000220.1, etc.
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$input_collapsed_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".bed")}'
        
        # CRITICAL FIX: Also split mutation file by chromosome to match collapsed.bed chunks
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$mutation_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".mut.bed")}'
        
        # Create processing script (updated to use per-chromosome mutation file)
        local process_script="$chunk_dir/run_cims_chunk.sh"
        cat > "$process_script" << 'CIMS_SCRIPT'
#!/bin/bash
chunk_file="$1"
iterations="$2"
output_dir="$3"
chunk_name=$(basename "$chunk_file" .bed)
# Use matching mutation file for this chromosome
mutation_chunk="${output_dir}/${chunk_name}.mut.bed"
output_chunk="${output_dir}/${chunk_name}.cims.txt"
cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cims_cache.XXXXXX")
export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
# Only run CIMS if both chunk files exist and are non-empty
if [[ -s "$chunk_file" && -s "$mutation_chunk" ]]; then
    CIMS.pl -big -c "$cache_dir" -n "$iterations" "$chunk_file" "$mutation_chunk" "$output_chunk" >/dev/null 2>&1
fi
rm -rf "$cache_dir" 2>/dev/null
CIMS_SCRIPT
        chmod +x "$process_script"
        
        # Calculate optimal parallel jobs based on RAM and file size
        local optimal_jobs=$(calculate_optimal_parallel_jobs "$threads" "$input_collapsed_bed")
        local parallel_jobs=$(( optimal_jobs < chr_count ? optimal_jobs : chr_count ))
        log_info "CIMS: Using $parallel_jobs/$threads threads based on available RAM"
        
        # Pass only chunk_file, iterations, and chunk_dir (mutation file is derived from chunk name)
        ls "$chunk_dir"/chr_*.bed 2>/dev/null | grep -v '\.mut\.bed$' | \
            parallel --memfree 4G -j "$parallel_jobs" \
            "$process_script" {} "$cims_iterations" "$chunk_dir"
        
        # Merge results: CIMS.pl writes its header last (starting with #), so we can't
        # use tail -n +2 (that would strip data lines). Instead collect all data lines
        # first, then prepend the header from whichever chunk had it.
        local merged_header=""
        > "$output_file"
        for result in "$chunk_dir"/*.cims.txt; do
            if [[ -s "$result" ]]; then
                [[ -z "$merged_header" ]] && merged_header=$(grep "^#" "$result")
                grep -v "^#" "$result" >> "$output_file"
            fi
        done
        # Prepend header so it sits at the top of the merged file
        if [[ -n "$merged_header" ]]; then
            local tmp_header=$(mktemp)
            echo "$merged_header" > "$tmp_header"
            cat "$output_file" >> "$tmp_header"
            mv "$tmp_header" "$output_file"
        fi

        # Cleanup
        rm -rf "$chunk_dir"
    else
        # Sequential mode (fallback)
        if [[ "$use_parallel" == "false" ]] && has_gnu_parallel; then
            log_info "Single chromosome detected, using sequential processing"
        elif ! has_gnu_parallel; then
            log_info "GNU parallel not available, using sequential processing"
        fi
        # Use -big -c for memory efficiency in sequential mode
        local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cims_cache.XXXXXX")
        # No -v: per-chromosome verbose output suppressed on console; captured in LOG_FILE
        local cmd="CIMS.pl -big -c '$cache_dir' -n $cims_iterations '$input_collapsed_bed' '$mutation_bed' '$output_file'"
        execute_cmd "$cmd" 2>>"${LOG_FILE:-/dev/null}"
        rm -rf "$cache_dir" 2>/dev/null

        # Normalize: CIMS.pl writes header last; move it to top so downstream tools see it first
        if [[ -s "$output_file" ]] && grep -q "^#" "$output_file"; then
            local tmp_norm=$(mktemp)
            grep "^#" "$output_file" > "$tmp_norm"
            grep -v "^#" "$output_file" >> "$tmp_norm"
            mv "$tmp_norm" "$output_file"
        fi
    fi

    # Process results
    if [[ -s "$output_file" ]]; then
        local raw_count=$(grep -v "^#" "$output_file" | wc -l)
        log_info "CIMS complete: $raw_count total sites in $output_file"

        # Filter for significance only if FDR < 1
        if (( $(echo "$cims_fdr < 1" | bc -l) )); then
            # Preserve raw output when -k is passed (for threshold exploration)
            if [[ "${KEEP_INTERMEDIATE:-no}" == "yes" ]]; then
                cp "$output_file" "${output_file%.txt}_raw.txt"
                log_info "CIMS raw output preserved: ${output_file%.txt}_raw.txt"
            fi
            log_info "Filtering CIMS results by FDR < $cims_fdr..."
            local temp_file="${output_file}.tmp"
            # Header is now at top (^#); sort only data lines, then reattach header
            local cims_header=$(grep "^#" "$output_file")
            {
                [[ -n "$cims_header" ]] && echo "$cims_header"
                grep -v "^#" "$output_file" | \
                    awk -F'\t' -v fdr="$cims_fdr" '$9+0 < fdr' | \
                    sort -k9,9n -k8,8nr -k7,7n
            } > "$temp_file"
            mv "$temp_file" "$output_file"

            local filtered_count=$(grep -v "^#" "$output_file" | wc -l)
            log_info "CIMS filtered: $filtered_count sites (FDR < $cims_fdr)"
        fi
    else
        echo -e "[WARNING] CIMS produced empty output: $output_file" >> "${LOG_FILE}"
        echo -ne "${YELLOW}[WARNING] Empty Output${NC} > "
        return 1
    fi
    
    return 0
}

# ── CITS analysis ───────────────────────────────────────────────────────────
# 6. CITS Analysis (CTK) - with parallel chromosome processing
# Detects crosslinking-induced truncation sites
# Input: collapsed BED + deletion file (used to EXCLUDE read-through tags)
# Output: CITS.txt with significant truncation sites (single-nucleotide singletons)
run_cits() {
    local input_collapsed_bed="$1"
    local deletion_bed="$2"           # Used to exclude read-through tags
    local output_file="$3"            # Should be .txt extension
    local cits_pvalue="${4:-0.001}"   # Default: p-value 0.001
    local cits_gap="${5:-25}"         # Default: gap 25 for clustering
    local threads="${THREADS:-4}"     # Use global THREADS or default
    
    log_info "Running CITS analysis..."
    log_info "Input tags: $input_collapsed_bed"
    log_info "Deletion file (to exclude): $deletion_bed"
    log_info "P-value: $cits_pvalue, Gap: $cits_gap"
    
    # Verify inputs exist
    if [[ ! -s "$input_collapsed_bed" ]]; then
        log_error "CITS: Collapsed BED file is empty or missing: $input_collapsed_bed"
        return 1
    fi
    if [[ ! -s "$deletion_bed" ]]; then
        log_warning "CITS: Deletion file is empty. Running without read-through filtering."
    fi
    
    # Set up CTK environment
    export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
    
    # Check if parallel processing is available and beneficial
    local use_parallel="false"
    local chr_count=$(cut -f1 "$input_collapsed_bed" | sort -u | wc -l)
    
    if has_gnu_parallel && [[ $chr_count -gt 1 ]]; then
        use_parallel="true"
        log_info "Using parallel processing ($chr_count chromosomes, $threads threads)"
    fi
    
    local cits_raw="${output_file%.txt}_tmp.bed"
    
    if [[ "$use_parallel" == "true" ]]; then
        # Parallel mode: split BOTH collapsed.bed AND deletion file by chromosome
        local chunk_dir=$(mktemp -d "${TMPDIR:-/tmp}/cits_parallel.XXXXXX")
        
        # Split collapsed BED by chromosome (standard chromosomes only: chr1-22, X, Y, M)
        # This filters out contigs like KI270737.1, GL000220.1, etc.
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$input_collapsed_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".bed")}'
        
        # CRITICAL FIX: Also split deletion file by chromosome to match collapsed.bed chunks
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$deletion_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".del.bed")}'
        
        # Create processing script (updated to use per-chromosome deletion file)
        local process_script="$chunk_dir/run_cits_chunk.sh"
        cat > "$process_script" << 'CITS_SCRIPT'
#!/bin/bash
chunk_file="$1"
pvalue="$2"
gap="$3"
output_dir="$4"
chunk_name=$(basename "$chunk_file" .bed)
# Use matching deletion file for this chromosome
deletion_chunk="${output_dir}/${chunk_name}.del.bed"
output_chunk="${output_dir}/${chunk_name}.cits.bed"
cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cits_cache.XXXXXX")
export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
# Only run CITS if collapsed chunk exists (deletion chunk may be absent/empty)
if [[ -s "$chunk_file" ]]; then
    actual_del="$deletion_chunk"
    if [[ ! -s "$deletion_chunk" ]]; then
        # No deletions for this chromosome: pass an empty file so CITS.pl
        # retains all reads as truncation candidates rather than skipping entirely
        actual_del=$(mktemp "${TMPDIR:-/tmp}/empty_del.XXXXXX")
    fi
    CITS.pl -big -c "$cache_dir" -p "$pvalue" --gap "$gap" "$chunk_file" "$actual_del" "$output_chunk" >/dev/null 2>&1
    [[ "$actual_del" != "$deletion_chunk" ]] && rm -f "$actual_del"
fi
rm -rf "$cache_dir" 2>/dev/null
CITS_SCRIPT
        chmod +x "$process_script"
        
        # Calculate optimal parallel jobs based on RAM and file size
        local optimal_jobs=$(calculate_optimal_parallel_jobs "$threads" "$input_collapsed_bed")
        local parallel_jobs=$(( optimal_jobs < chr_count ? optimal_jobs : chr_count ))
        log_info "CITS: Using $parallel_jobs/$threads threads based on available RAM"
        
        # Pass only chunk_file, pvalue, gap, and chunk_dir (deletion file is derived from chunk name)
        ls "$chunk_dir"/chr_*.bed 2>/dev/null | grep -v '\.del\.bed$' | \
            parallel --memfree 4G -j "$parallel_jobs" \
            "$process_script" {} "$cits_pvalue" "$cits_gap" "$chunk_dir"
        
        # Merge results
        > "$cits_raw"
        for result in "$chunk_dir"/*.cits.bed; do
            if [[ -s "$result" ]]; then
                cat "$result" >> "$cits_raw"
            fi
        done
        
        # Cleanup
        rm -rf "$chunk_dir"
    else
        # Sequential mode (fallback)
        if [[ "$use_parallel" == "false" ]] && has_gnu_parallel; then
            log_info "Single chromosome detected, using sequential processing"
        elif ! has_gnu_parallel; then
            log_info "GNU parallel not available, using sequential processing"
        fi
        # Use -big -c for memory efficiency in sequential mode
        local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cits_cache.XXXXXX")
        # No -v: per-chromosome verbose output suppressed on console; captured in LOG_FILE
        local cmd="CITS.pl -big -c '$cache_dir' -p $cits_pvalue --gap $cits_gap '$input_collapsed_bed' '$deletion_bed' '$cits_raw'"
        execute_cmd "$cmd"
        rm -rf "$cache_dir" 2>/dev/null
    fi
    
    # Process results
    if [[ -s "$cits_raw" ]]; then
        local raw_count=$(wc -l < "$cits_raw")
        log_info "CITS raw output: $raw_count sites"
        
        # Preserve raw output when -k is passed (for threshold exploration)
        if [[ "${KEEP_INTERMEDIATE:-no}" == "yes" ]]; then
            cp "$cits_raw" "${output_file%.txt}_raw.bed"
            log_info "CITS raw output preserved: ${output_file%.txt}_raw.bed"
        fi
        
        # Filter to single-nucleotide sites (singleton) - this is the main output
        # Add header for consistency with CIMS output
        echo -e "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand" > "$output_file"
        awk '{if($3-$2==1) {print $0}}' "$cits_raw" >> "$output_file"
        
        local singleton_count=$(($(wc -l < "$output_file") - 1))  # Subtract header
        log_info "CITS complete: $singleton_count singleton sites in $output_file"
        
        # Remove intermediate raw file
        rm -f "$cits_raw"
    else
        echo -e "[WARNING] CITS produced empty output" >> "${LOG_FILE}"
        echo -ne "${YELLOW}[WARNING] Empty Output${NC} > "
        return 1
    fi
    
    return 0
}

# ── Flanked BED for motif analysis ──────────────────────────────────────────
# 7. Flanked BED Generation for Motif Analysis
# Generates ±10nt flanked regions around CIMS/CITS sites
# Creates flanked BED file alongside the input file (same directory)
# Users can run their own motif analysis tools on these files
generate_flanked_bed() {
    local input_bed="$1"
    local flank_nt="${2:-10}"         # Default: ±10 nucleotides
    
    if [[ ! -s "$input_bed" ]]; then
        return 1
    fi
    
    # Create flanked file alongside input: sample_CIMS_sub.txt → sample_CIMS_sub_flanked.bed
    local flanked_bed="${input_bed%.txt}_flanked.bed"
    
    awk -v n="$flank_nt" 'BEGIN{OFS="\t"} {
        start = $2 - n
        if (start < 0) start = 0
        print $1, start, $3 + n, $4, $5, $6
    }' "$input_bed" > "$flanked_bed"
    
    local site_count=$(wc -l < "$flanked_bed")
    log_info "Generated flanked BED (±${flank_nt}nt, $site_count sites): $(basename "$flanked_bed")"
}

# ── Full CTK pipeline (BAM → CIMS/CITS) ────────────────────────────────────
# 8. Full CTK Analysis Pipeline
# Orchestrates the complete CIMS/CITS workflow based on RUN_CIMS and RUN_CITS flags
run_ctk_full_analysis() {
    local bam_file="$1"
    local output_dir="$2"
    local genome_fasta="$3"
    local cims_iterations="${4:-10}"
    local cims_fdr="${5:-1}"
    local cits_pvalue="${6:-1}"
    local cits_gap="${7:-25}"
    local motif_flank="${8:-10}"
    local run_motif="${9:-yes}"
    local run_cims="${10:-true}"
    local run_cits="${11:-true}"
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK CIMS/CITS FULL ANALYSIS"
    log_info "  CIMS: $run_cims | CITS: $run_cits"
    log_info "═══════════════════════════════════════════════════════════════"
    
    # Create directories based on what's enabled
    mkdir -p "$output_dir/preprocessing"
    [[ "$run_cims" == "true" ]] && mkdir -p "$output_dir/CIMS"
    [[ "$run_cits" == "true" ]] && mkdir -p "$output_dir/CITS"
    [[ "$run_motif" == "yes" ]] && mkdir -p "$output_dir/motif_analysis"
    
    local sample_name=$(basename "${bam_file%.bam}" | sed 's/.Aligned.sortedByCoord.out//')
    
    # Phase 1: Preprocessing
    log_info "Phase 1: Parsing alignment..."
    local tags_bed="${output_dir}/preprocessing/${sample_name}_tags.bed"
    local mutation_file="${output_dir}/preprocessing/${sample_name}_mutations.txt"
    
    # run_parse_alignment handles calmd + parseAlignment.pl
    run_parse_alignment "$bam_file" "$tags_bed" "$mutation_file" "$(dirname "$genome_fasta")"
    
    if [[ ! -s "$tags_bed" ]]; then
        log_error "parseAlignment.pl failed to produce tags. Aborting CTK analysis."
        return 1
    fi
    
    # Phase 1b: Collapse tags
    log_info "Phase 1b: Collapsing PCR duplicates..."
    local collapsed_bed="${output_dir}/preprocessing/${sample_name}_collapsed.bed"
    tag2collapse.pl --keep-tag-name "$tags_bed" "$collapsed_bed"
    
    if [[ ! -s "$collapsed_bed" ]]; then
        log_error "tag2collapse.pl failed. Aborting CTK analysis."
        return 1
    fi
    
    # Phase 1c: CTK Preprocessing (selectRow + getMutationType)
    log_info "Phase 1c: Filtering and extracting mutation types..."
    run_ctk_preprocessing "$collapsed_bed" "$mutation_file" "${output_dir}/preprocessing"
    
    local del_bed="${output_dir}/preprocessing/deletions.bed"
    local sub_bed="${output_dir}/preprocessing/substitutions.bed"
    
    # Phase 2: CIMS Analysis (only if enabled)
    if [[ "$run_cims" == "true" ]]; then
        log_info "Phase 2: CIMS Analysis..."
        
        if [[ -s "$del_bed" ]]; then
            log_info "Running CIMS on deletions..."
            run_cims "$collapsed_bed" "$del_bed" \
                "${output_dir}/CIMS/${sample_name}_CIMS_del.txt" \
                "$cims_iterations" "$cims_fdr"
        fi
        
        if [[ -s "$sub_bed" ]]; then
            log_info "Running CIMS on substitutions..."
            run_cims "$collapsed_bed" "$sub_bed" \
                "${output_dir}/CIMS/${sample_name}_CIMS_sub.txt" \
                "$cims_iterations" "$cims_fdr"
        fi
    else
        log_info "Phase 2: CIMS Analysis... SKIPPED (not enabled)"
    fi
    
    # Phase 3: CITS Analysis (only if enabled)
    # run_cits handles missing/empty deletion file internally (runs without read-through filter).
    # For standard iCLIP, deletions at crosslink sites are rare — CITS must still run.
    # Do NOT gate CITS on [[ -s "$del_bed" ]]: that silently skips all iCLIP samples.
    if [[ "$run_cits" == "true" ]]; then
        log_info "Phase 3: CITS Analysis..."
        run_cits "$collapsed_bed" "$del_bed" \
            "${output_dir}/CITS/${sample_name}_CITS.bed" \
            "$cits_pvalue" "$cits_gap"
    else
        log_info "Phase 3: CITS Analysis... SKIPPED (not enabled)"
    fi
    
    # Phase 4: Flanked BED Generation (for user's motif analysis)
    if [[ "$run_motif" == "yes" ]]; then
        log_info "Phase 4: Generating flanked BED files..."
        
        if [[ "$run_cims" == "true" ]]; then
            local cims_del_sig="${output_dir}/CIMS/${sample_name}_CIMS_del_significant.bed"
            local cims_sub_sig="${output_dir}/CIMS/${sample_name}_CIMS_sub_significant.bed"
            [[ -s "$cims_del_sig" ]] && generate_flanked_bed "$cims_del_sig" "$motif_flank"
            [[ -s "$cims_sub_sig" ]] && generate_flanked_bed "$cims_sub_sig" "$motif_flank"
        fi
        
        if [[ "$run_cits" == "true" ]]; then
            local cits_singleton="${output_dir}/CITS/${sample_name}_CITS_singleton.bed"
            [[ -s "$cits_singleton" ]] && generate_flanked_bed "$cits_singleton" "$motif_flank"
        fi
    fi
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK ANALYSIS COMPLETE"
    log_info "  Output: $output_dir"
    log_info "═══════════════════════════════════════════════════════════════"
}

# ── Streamlined CTK analysis (reuses existing collapsed.bed) ────────────────
# 9. CTK Analysis Pipeline (Streamlined - uses pre-existing collapsed.bed and mutations.txt)
# This function reuses the standard pipeline's outputs to avoid duplicate preprocessing
run_ctk_analysis() {
    local collapsed_bed="$1"          # From standard pipeline
    local mutation_file="$2"          # From standard pipeline
    local output_dir="$3"
    local genome_fasta="$4"
    local sample_name="$5"
    local cims_iterations="${6:-10}"
    local cims_fdr="${7:-0.05}"
    local cits_pvalue="${8:-0.05}"
    local cits_gap="${9:-25}"
    local motif_flank="${10:-10}"
    local run_motif="${11:-yes}"
    local run_cims="${12:-true}"
    local run_cits="${13:-true}"
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK CIMS/CITS ANALYSIS"
    log_info "  Sample: $sample_name"
    log_info "  CIMS: $run_cims | CITS: $run_cits"
    log_info "═══════════════════════════════════════════════════════════════"
    
    # Create directories based on what's enabled
    [[ "$run_cims" == "true" ]] && mkdir -p "$output_dir/CIMS"
    [[ "$run_cits" == "true" ]] && mkdir -p "$output_dir/CITS"
    
    # Verify inputs exist
    if [[ ! -s "$collapsed_bed" ]]; then
        log_error "CTK: Collapsed BED file is empty or missing: $collapsed_bed"
        return 1
    fi
    if [[ ! -s "$mutation_file" ]]; then
        log_warning "CTK: Mutation file is empty or missing: $mutation_file"
        log_warning "CTK: CIMS/CITS requires mutation information. Skipping."
        return 1
    fi
    
    # Step 1: CTK Preprocessing (selectRow + getMutationType)
    log_info "Step 1: Filtering and extracting mutation types..."
    run_ctk_preprocessing "$collapsed_bed" "$mutation_file" "$output_dir"
    
    local del_bed="${output_dir}/deletions.bed"
    local sub_bed="${output_dir}/substitutions.bed"
    
    # Step 2: CIMS Analysis (only if enabled)
    if [[ "$run_cims" == "true" ]]; then
        update_status "CIMS"
        
        local cims_del_file="${output_dir}/CIMS/${sample_name}_CIMS_del.txt"
        local cims_sub_file="${output_dir}/CIMS/${sample_name}_CIMS_sub.txt"
        
        # CIMS on deletions: the primary signal for standard/iCLIP crosslink-induced mutations.
        # CTK is designed around deletion-based CIMS (BWA alignment produces deletions readily).
        # STAR EndToEnd tends to realign deletion reads as substitutions due to gap penalties,
        # so deletion counts may be very low with STAR. If del_bed is empty, CIMS is skipped.
        if [[ -s "$del_bed" ]]; then
            run_cims "$collapsed_bed" "$del_bed" "$cims_del_file" \
                "$cims_iterations" "$cims_fdr"
        else
            log_warning "CIMS: No deletions found — CIMS deletion analysis skipped."
            log_warning "  This is expected with STAR EndToEnd alignment on short reads."
            log_warning "  Consider BWA (--mapper bowtie2 or external BWA) for richer deletion signal."
        fi

        # CIMS on substitutions: secondary signal; useful for C→T transitions (iCLIP crosslink signature).
        # CIMS.pl is patched (tagNum==0 guard + count>0 q-value guard) so this is safe to run.
        if [[ -s "$sub_bed" ]]; then
            run_cims "$collapsed_bed" "$sub_bed" "$cims_sub_file" \
                "$cims_iterations" "$cims_fdr"
        fi

        # Generate flanked BED for CIMS results (for user's motif analysis)
        if [[ "$run_motif" == "yes" ]]; then
            [[ -s "$cims_del_file" ]] && generate_flanked_bed "$cims_del_file" "$motif_flank"
            [[ -s "$cims_sub_file" ]] && generate_flanked_bed "$cims_sub_file" "$motif_flank"
        fi
    fi
    
    # Step 3: CITS Analysis (only if enabled)
    if [[ "$run_cits" == "true" ]]; then
        update_status "CITS"

        local cits_file="${output_dir}/CITS/${sample_name}_CITS.txt"

        # run_cits handles missing/empty deletion file internally (runs without read-through filter).
        # Do NOT gate CITS on [[ -s "$del_bed" ]]: that silently skips all standard iCLIP samples.
        run_cits "$collapsed_bed" "$del_bed" "$cits_file" \
            "$cits_pvalue" "$cits_gap"

        # Generate flanked BED for CITS results (for user's motif analysis)
        if [[ "$run_motif" == "yes" ]]; then
            [[ -s "$cits_file" ]] && generate_flanked_bed "$cits_file" "$motif_flank"
        fi
    fi
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK ANALYSIS COMPLETE"
    log_info "  Output: $output_dir"
    log_info "═══════════════════════════════════════════════════════════════"
}

# ── Group CTK analysis ──────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════
# Group-Based CTK Analysis (Bash 3.x Compatible)
# Aggregates samples by group before running CIMS/CITS
# ═══════════════════════════════════════════════════════════════════════════

run_group_ctk_analysis() {
    local groups_file="$1"
    local output_root="$2"
    local genome_index="$3"
    local cims_iterations="$4"
    local cims_fdr="$5"
    local cits_pvalue="$6"
    local cits_gap="$7"
    local motif_flank="$8"
    local run_motif="$9"
    local run_cims="${10}"
    local run_cits="${11}"
    local work_dir="${12:-$(pwd)}"   # directory containing *_analysis subdirs
    
    console_msg "\n[GROUP CTK ANALYSIS]"
    
    # 1. Parse groups file into temp file (sample<TAB>group format)
    local groups_map=$(mktemp)
    parse_groups_file "$groups_file" "$groups_map"
    
    # 2. Create temp file for group→samples mapping
    local group_samples_file=$(mktemp)
    
    # 3. Find all sample analysis directories (in work_dir) and assign to groups
    for sample_dir in "${work_dir}"/*_analysis; do
        [[ ! -d "$sample_dir" ]] && continue
        
        local sample_name="$(basename "${sample_dir%_analysis}")"
        
        # Look up group for this sample using grep (bash 3.x compatible)
        local group=$(grep -w "^${sample_name}" "$groups_map" 2>/dev/null | cut -f2)
        
        if [[ -z "$group" ]]; then
            # Sample not in groups file → individual group (use sample name)
            group="$sample_name"
            log_info "Sample '$sample_name' not in groups file, treating as individual"
        fi
        
        # Append to group_samples_file: group<TAB>sample
        echo -e "${group}\t${sample_name}" >> "$group_samples_file"
    done
    
    # 4. Determine CTK output folder name
    local ctk_folder_name
    if [[ "$run_cims" == "true" ]] && [[ "$run_cits" == "true" ]]; then
        ctk_folder_name="5_CTK_Analysis"
    elif [[ "$run_cims" == "true" ]]; then
        ctk_folder_name="5_CIMS_Analysis"
    elif [[ "$run_cits" == "true" ]]; then
        ctk_folder_name="5_CITS_Analysis"
    fi
    
    # 5. Get unique groups and process alphabetically (skip unknown)
    local unique_groups=$(cut -f1 "$group_samples_file" | sort -u | grep -v "^unknown$")
    
    for group in $unique_groups; do
        # Get all samples for this group
        local samples=$(grep -w "^${group}" "$group_samples_file" | cut -f2 | tr '\n' ' ')
        local sample_count=$(echo $samples | wc -w | tr -d ' ')
        
        # Print group header on same line (no newline) so status updates appear after it
        if [[ "$sample_count" -eq 1 ]]; then
            printf "  > Processing %s: " "$group"
        else
            printf "  > Processing %s (%d samples): " "$group" "$sample_count"
        fi
        
        # 6. Create group CTK directory
        local group_ctk_dir="$output_root/$ctk_folder_name/$group"
        mkdir -p "$group_ctk_dir"
        
        # 7. Aggregate collapsed.bed files
        local group_collapsed="$group_ctk_dir/${group}_collapsed.bed"
        > "$group_collapsed"  # Clear file
        for sample in $samples; do
            local sample_dir="${work_dir}/${sample}_analysis"
            if [[ -d "$sample_dir" ]]; then
                local sample_collapsed=$(find "$sample_dir" -name "*_collapsed.bed" 2>/dev/null | head -n 1)
                if [[ -s "$sample_collapsed" ]]; then
                    cat "$sample_collapsed" >> "$group_collapsed"
                    log_info "Added $sample_collapsed to group $group"
                else
                    log_warning "Collapsed BED not found for sample: $sample"
                fi
            else
                log_warning "Analysis directory not found: $sample_dir"
            fi
        done
        
        # 8. Aggregate mutations.txt files
        local group_mutations="$group_ctk_dir/${group}_mutations.txt"
        > "$group_mutations"  # Clear file
        for sample in $samples; do
            local sample_dir="${work_dir}/${sample}_analysis"
            if [[ -d "$sample_dir" ]]; then
                local sample_mutations=$(find "$sample_dir" -name "*_mutations.txt" 2>/dev/null | head -n 1)
                if [[ -s "$sample_mutations" ]]; then
                    cat "$sample_mutations" >> "$group_mutations"
                fi
            fi
        done
        
        # 9. Get genome fasta for motif analysis (prioritize genome/primary, exclude ncrna)
        local genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*genome*.fa" -o -name "*genome*.fasta" \) 2>/dev/null | head -n 1)
        if [[ -z "$genome_fasta" ]]; then
            genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*primary*.fa" -o -name "*primary*.fasta" \) 2>/dev/null | head -n 1)
        fi
        if [[ -z "$genome_fasta" ]]; then
            genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) ! -name "*ncrna*" 2>/dev/null | head -n 1)
        fi
        
        # 10. Run CTK analysis on aggregated data
        if [[ -s "$group_collapsed" ]]; then
            run_ctk_analysis "$group_collapsed" "$group_mutations" \
                "$group_ctk_dir" "$genome_fasta" "$group" \
                "$cims_iterations" "$cims_fdr" "$cits_pvalue" "$cits_gap" \
                "$motif_flank" "$run_motif" "$run_cims" "$run_cits"
        else
            log_error "CTK: Collapsed BED file is empty or missing: $group_collapsed"
        fi
        
        update_status_done
    done
    
    # Cleanup temp files
    rm -f "$groups_map" "$group_samples_file"
    
    console_msg "  > Group CTK analysis complete"
}
