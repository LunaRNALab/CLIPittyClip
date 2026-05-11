#!/bin/bash
# lib/clink.sh — CLIPittyClip v3.4
# Part of the CLIPittyClip pipeline. Source via lib/modules.sh or directly.
# Auto-split from modules.sh by build_v34_modules.sh.

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ═══════════════════════════════════════════════════════════════════════════════
# CLINK PIPELINE  (UMI-tools + Python pileup engine)
# Replaces: tag2collapse.pl + parseAlignment.pl + CIMS.pl + CITS.pl
# ═══════════════════════════════════════════════════════════════════════════════

# Resolve the path to the Clink Python scripts bundled in lib/clink/
_clink_dir() {
    echo "$(dirname "${BASH_SOURCE[0]}")/clink"
}

# Execute a Clink Python command, respecting VERBOSE (mirrors execute_cmd behavior).
# All Clink Python scripts print progress to stderr; with VERBOSE=true those lines
# are teed to the console so the user can watch in real time.
_clink_exec() {
    local log="${LOG_FILE:-/dev/null}"
    if [[ "${VERBOSE:-false}" == "true" ]]; then
        eval "$*" 2>&1 | tee -a "$log"
        return "${PIPESTATUS[0]}"
    else
        eval "$*" >> "$log" 2>&1
        return $?
    fi
}

# ---------------------------------------------------------------------------
# check_clink_deps — hard dependency check, called once before any Clink work
# ---------------------------------------------------------------------------
check_clink_deps() {
    local ok=true

    # Python packages: pysam, numpy, scipy
    if ! python3 -c "import pysam, numpy, scipy" 2>/dev/null; then
        log_error "Clink requires Python packages: pysam, numpy, scipy"
        log_error "  Install with: conda install -n clipittyclip pysam umi_tools -c bioconda"
        ok=false
    fi

    # umi_tools binary: check PATH first, then common conda env locations
    local umi_bin
    umi_bin=$(command -v umi_tools 2>/dev/null)
    if [[ -z "$umi_bin" ]]; then
        local candidates=(
            "$CONDA_PREFIX/bin/umi_tools"
            "$HOME/anaconda3/envs/umi_tools/bin/umi_tools"
            "$HOME/miniconda3/envs/umi_tools/bin/umi_tools"
            "/opt/anaconda3/envs/umi_tools/bin/umi_tools"
            "/opt/miniconda3/envs/umi_tools/bin/umi_tools"
        )
        for p in "${candidates[@]}"; do
            if [[ -x "$p" ]]; then umi_bin="$p"; break; fi
        done
    fi
    if [[ -z "$umi_bin" ]]; then
        log_error "Clink requires umi_tools. Not found in PATH or common conda locations."
        log_error "  Install with: conda install -n clipittyclip umi_tools -c bioconda"
        log_error "  Or in a separate env: conda create -n umi_tools -c bioconda umi_tools python=3.11"
        ok=false
    else
        export CLINK_UMI_TOOLS="$umi_bin"
        log_info "umi_tools found: $umi_bin"
    fi

    [[ "$ok" == "true" ]]
}

# ---------------------------------------------------------------------------
# run_clink_collapse — BAM-level UMI deduplication via umi_tools dedup
#
# Args:
#   $1  input sorted BAM
#   $2  output deduplicated BAM path
#   $3  UMI length (0 = position-only; -1 = auto-detect)
#   $4  threads (default 4)
# ---------------------------------------------------------------------------
run_clink_collapse() {
    local bam_in="$1"
    local bam_out="$2"
    local umi_len="${3:--1}"
    local threads="${4:-4}"
    local clink_dir
    clink_dir=$(_clink_dir)

    log_info "Clink collapse: UMI-aware deduplication (umi_tools directional)"
    log_info "  Input:  $bam_in"
    log_info "  Output: $bam_out"

    local extra_args=""
    if [[ -n "${CLINK_UMI_TOOLS:-}" ]]; then
        extra_args="--umi-tools $CLINK_UMI_TOOLS"
    fi
    if [[ "$umi_len" -ge 0 ]] 2>/dev/null; then
        extra_args="$extra_args --umi-len $umi_len"
    fi

    _clink_exec python3 "$clink_dir/collapse.py" \
        --bam     "$bam_in" \
        --out     "$bam_out" \
        --threads "$threads" \
        $extra_args

    if [[ $? -ne 0 ]] || [[ ! -s "$bam_out" ]]; then
        log_error "Clink collapse failed. Check log for details."
        return 1
    fi
    log_info "Clink collapse complete: $bam_out"
}

# ---------------------------------------------------------------------------
# run_clink_pileup — single BAM scan → compressed NPZ
#
# Args:
#   $1  deduplicated BAM
#   $2  output .npz path
#   $3  threads (default 1; passed to pileup.py --threads for chromosome-level parallelism)
# ---------------------------------------------------------------------------
run_clink_pileup() {
    local bam_in="$1"
    local npz_out="$2"
    local threads="${3:-1}"
    local clink_dir
    clink_dir=$(_clink_dir)

    log_info "Clink pileup: scanning BAM → $npz_out (threads=$threads)"

    _clink_exec python3 "$clink_dir/pileup.py" \
        "$bam_in" \
        --out "$npz_out" \
        --threads "$threads"

    if [[ $? -ne 0 ]] || [[ ! -s "$npz_out" ]]; then
        log_error "Clink pileup failed. Check log for details."
        return 1
    fi
    log_info "Clink pileup saved: $npz_out"
}

# ---------------------------------------------------------------------------
# run_clink_cits — truncation site calling from NPZ
#
# Args:
#   $1  pileup .npz
#   $2  output prefix
#   $3  min coverage (default 5)
#   $4  min fraction (default 0.05)
#   $5  FDR threshold (default 0.05)
# ---------------------------------------------------------------------------
run_clink_cits() {
    local npz_in="$1"
    local prefix="$2"
    local min_cov="${3:-5}"
    local min_frac="${4:-0.05}"
    local fdr="${5:-0.05}"
    local clink_dir
    clink_dir=$(_clink_dir)

    log_info "Clink CITS: truncation site calling"

    _clink_exec python3 "$clink_dir/cits.py" \
        --pileup  "$npz_in" \
        --prefix  "$prefix" \
        --min-cov "$min_cov" \
        --min-frac "$min_frac" \
        --fdr     "$fdr"

    if [[ $? -ne 0 ]]; then
        log_error "Clink CITS failed. Check log for details."
        return 1
    fi

    local out_bed="${prefix}_truncations.bed"
    local n=0
    [[ -f "$out_bed" ]] && n=$(grep -c '' "$out_bed" 2>/dev/null || echo 0)
    log_info "Clink CITS complete: $(( n > 0 ? n - 1 : 0 )) significant sites → $out_bed"
}

# ---------------------------------------------------------------------------
# run_clink_cims — deletion + substitution site calling from NPZ
#
# Args:
#   $1  pileup .npz
#   $2  output prefix
#   $3  min coverage (default 5)
#   $4  min fraction (default 0.05)
#   $5  FDR threshold (default 0.05)
#   $6  sub-types to output (e.g. "TC" for PAR-CLIP; "" = all)
#   $7  no-subs flag: "true" = deletions only
# ---------------------------------------------------------------------------
run_clink_cims() {
    local npz_in="$1"
    local prefix="$2"
    local min_cov="${3:-5}"
    local min_frac="${4:-0.05}"
    local fdr="${5:-0.05}"
    local sub_types="${6:-}"
    local no_subs="${7:-false}"
    local clink_dir
    clink_dir=$(_clink_dir)

    log_info "Clink CIMS: deletion + substitution site calling"

    local extra_args=""
    [[ "$no_subs" == "true" ]]   && extra_args="$extra_args --no-subs"
    [[ -n "$sub_types" ]]        && extra_args="$extra_args --sub-types $sub_types"

    _clink_exec python3 "$clink_dir/cims.py" \
        --pileup   "$npz_in" \
        --prefix   "$prefix" \
        --min-cov  "$min_cov" \
        --min-frac "$min_frac" \
        --fdr      "$fdr" \
        $extra_args

    if [[ $? -ne 0 ]]; then
        log_error "Clink CIMS failed. Check log for details."
        return 1
    fi
    log_info "Clink CIMS complete → ${prefix}_deletions.bed + substitution beds"
}

# ---------------------------------------------------------------------------
# run_clink_full — orchestrates the full Clink sub-pipeline for one sample
#
# Args:
#   $1  input sorted BAM (from STAR/Bowtie2)
#   $2  output directory  (e.g. 6_Clink/{sample}/)
#   $3  sample name
#   $4  UMI length (-1 = auto)
#   $5  threads
#   $6  run_cits  "true"|"false"
#   $7  run_cims  "true"|"false"
#   $8  min coverage
#   $9  min fraction
#   $10 FDR threshold
# ---------------------------------------------------------------------------
run_clink_full() {
    local bam_in="$1"
    local out_dir="$2"
    local sample_name="$3"
    local umi_len="${4:--1}"
    local threads="${5:-4}"
    local run_cits="${6:-true}"
    local run_cims="${7:-true}"
    local min_cov="${8:-5}"
    local min_frac="${9:-0.05}"
    local fdr="${10:-0.05}"
    local prebuilt_dedup="${11:-}"  # optional: pre-existing dedup BAM from early Clink-only path

    mkdir -p "$out_dir"

    local dedup_bam="${out_dir}/${sample_name}_dedup.bam"
    local npz="${out_dir}/${sample_name}_pileup.npz"
    local prefix="${out_dir}/${sample_name}"

    # Collapse step: skip if dedup BAM was already built in the early Clink-only path
    if [[ -n "$prebuilt_dedup" ]] && [[ -f "$prebuilt_dedup" ]]; then
        log_info "Clink collapse: reusing pre-built dedup BAM: $prebuilt_dedup"
        dedup_bam="$prebuilt_dedup"
    else
        update_status "Clink collapse"
        run_clink_collapse "$bam_in" "$dedup_bam" "$umi_len" "$threads" || return 1
    fi

    update_status "Clink pileup"
    run_clink_pileup "$dedup_bam" "$npz" "$threads" || return 1

    if [[ "$run_cits" == "true" ]]; then
        update_status "Clink CITS"
        run_clink_cits "$npz" "$prefix" "$min_cov" "$min_frac" "$fdr" || true
    fi

    if [[ "$run_cims" == "true" ]]; then
        update_status "Clink CIMS"
        run_clink_cims "$npz" "$prefix" "$min_cov" "$min_frac" "$fdr" || true
    fi

    # Generate collapsed BED for peak calling (bedtools bamtobed -split handles spliced reads)
    update_status "Clink BED"
    local collapsed_bed="${out_dir}/${sample_name}_collapsed.bed"
    bedtools bamtobed -i "$dedup_bam" -split 2>/dev/null \
        | sort -k1,1 -k2,2n > "$collapsed_bed"
    if [[ -s "$collapsed_bed" ]]; then
        log_info "Clink collapsed BED: $collapsed_bed"
    else
        log_warning "Clink: bedtools bamtobed produced empty BED — peak calling may be affected"
    fi

    update_status_done
    log_info "Clink full pipeline complete for $sample_name"
    echo "$collapsed_bed"   # Return path for caller to use in peak calling
}

# ---------------------------------------------------------------------------
# run_group_clink_analysis — grouped Clink CITS/CIMS (--group-xlsite)
#
# Merges per-sample dedup BAMs by group, then runs pileup → CITS/CIMS on the
# merged BAM.  Mirrors run_group_ctk_analysis but for the Clink sub-pipeline.
#
# Args:
#   $1  groups_file         — same groups.txt used by CTK/bedgraph
#   $2  clink_output_root   — e.g. OUTPUT_ROOT/5_Clink  or  OUTPUT_ROOT/6_Clink
#                             (per-sample sub-dirs live here after aggregation)
#   $3  threads             (default 4)
#   $4  run_cits            "true"|"false"
#   $5  run_cims            "true"|"false"
#   $6  min_coverage        (default 5)
#   $7  min_fraction        (default 0.05)
#   $8  fdr                 (default 0.05)
# ---------------------------------------------------------------------------
run_group_clink_analysis() {
    local groups_file="$1"
    local clink_output_root="$2"
    local threads="${3:-4}"
    local run_cits="${4:-true}"
    local run_cims="${5:-true}"
    local min_cov="${6:-5}"
    local min_frac="${7:-0.05}"
    local fdr="${8:-0.05}"

    console_msg "\n[GROUP CLINK ANALYSIS]"

    # --- Parse groups file → sample→group map ---
    local groups_map
    groups_map=$(mktemp)
    parse_groups_file "$groups_file" "$groups_map"

    # --- Scan clink_output_root for per-sample dirs; skip GROUP_* dirs ---
    local group_samples_file
    group_samples_file=$(mktemp)

    for sample_dir in "$clink_output_root"/*/; do
        [[ ! -d "$sample_dir" ]] && continue
        local sample_name
        sample_name="$(basename "$sample_dir")"
        # Skip group output dirs from this (or a previous) run
        [[ "$sample_name" == GROUP_* ]] && continue   # unquoted RHS = glob pattern in [[ ]]

        local group
        group=$(grep -w "^${sample_name}" "$groups_map" 2>/dev/null | cut -f2)
        if [[ -z "$group" ]]; then
            group="$sample_name"
            log_info "Clink group: '$sample_name' not in groups file → treated as individual"
        fi
        printf "%s\t%s\n" "$group" "$sample_name" >> "$group_samples_file"
    done

    if [[ ! -s "$group_samples_file" ]]; then
        log_warning "Clink group analysis: no sample directories found in $clink_output_root"
        rm -f "$groups_map" "$group_samples_file"
        return 0
    fi

    local unique_groups
    unique_groups=$(cut -f1 "$group_samples_file" | sort -u)

    for group in $unique_groups; do
        local samples
        samples=$(grep -w "^${group}" "$group_samples_file" | cut -f2 | tr '\n' ' ')
        local sample_count
        sample_count=$(echo $samples | wc -w | tr -d ' ')

        printf "  > %s (%d sample%s): " \
            "$group" "$sample_count" "$([[ $sample_count -eq 1 ]] && echo '' || echo 's')"

        local group_dir="$clink_output_root/GROUP_${group}"
        mkdir -p "$group_dir"

        # --- Collect per-sample dedup BAMs ---
        local bam_inputs=()
        for sample in $samples; do
            local bam="$clink_output_root/${sample}/${sample}_dedup.bam"
            if [[ -f "$bam" ]]; then
                bam_inputs+=("$bam")
                log_info "  Clink group $group: adding $bam"
            else
                log_warning "  Clink group $group: dedup BAM not found for $sample ($bam)"
            fi
        done

        if [[ ${#bam_inputs[@]} -eq 0 ]]; then
            log_warning "Clink group $group: no dedup BAMs found — skipping"
            printf "SKIPPED (no BAMs)\n"
            continue
        fi

        # --- Merge BAMs (or symlink if only one sample) ---
        local merged_bam="$group_dir/${group}_merged_dedup.bam"
        if [[ ${#bam_inputs[@]} -eq 1 ]]; then
            cp "${bam_inputs[0]}" "$merged_bam"
        else
            update_status "Clink group $group merge"
            samtools merge -f -@ "$threads" "$merged_bam" "${bam_inputs[@]}" \
                2>>"${LOG_FILE:-/dev/null}"
        fi

        if [[ ! -s "$merged_bam" ]]; then
            log_error "Clink group $group: BAM merge failed → $merged_bam"
            printf "FAILED\n"
            continue
        fi

        # Index merged BAM (required by pileup.py)
        samtools index "$merged_bam" 2>>"${LOG_FILE:-/dev/null}"

        # --- Pileup ---
        local npz="$group_dir/${group}_pileup.npz"
        update_status "Clink group $group pileup"
        if ! run_clink_pileup "$merged_bam" "$npz" "$threads"; then
            log_error "Clink group $group: pileup failed"
            printf "FAILED\n"
            continue
        fi

        local prefix="$group_dir/${group}"

        # --- CITS ---
        if [[ "$run_cits" == "true" ]]; then
            update_status "Clink group $group CITS"
            run_clink_cits "$npz" "$prefix" "$min_cov" "$min_frac" "$fdr" || true
        fi

        # --- CIMS ---
        if [[ "$run_cims" == "true" ]]; then
            update_status "Clink group $group CIMS"
            run_clink_cims "$npz" "$prefix" "$min_cov" "$min_frac" "$fdr" || true
        fi

        update_status_done
        log_info "Clink group $group complete → $group_dir"
    done

    rm -f "$groups_map" "$group_samples_file"
}
