#!/bin/bash
# lib/peaks.sh — CLIPittyClip v3.4
# Part of the CLIPittyClip pipeline. Source via lib/modules.sh or directly.
# Auto-split from modules.sh by build_v34_modules.sh.

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ── HOMER peak calling ──────────────────────────────────────────────────────
# 4. Peak Calling (HOMER)
run_peak_calling_homer() {
    local input_bed="$1"
    local out_dir="$2"
    local peak_dist="$3"
    local peak_size="$4"
    local frag_len="$5"
    local log_file="${out_dir}_homer.log"

    update_status "Peaks"
    log_info "Calling peaks with HOMER..."

    echo "Running HOMER makeTagDirectory..." > "$log_file"
    makeTagDirectory "${out_dir}" "${input_bed}" -single -format bed >> "$log_file" 2>&1

    echo "Running HOMER findPeaks..." >> "$log_file"
    findPeaks "${out_dir}" -o auto -style factor -L 2 -localSize 10000 -strand separate \
        -minDist "${peak_dist}" -size "${peak_size}" -fragLength "${frag_len}" $ADV_PEAK_CALLER_ARGS >> "$log_file" 2>&1

    if [[ -f "${out_dir}/peaks.txt" ]]; then
        echo "Converting peaks.txt to BED format..." >> "$log_file"
        sed '/^[[:blank:]]*#/d;s/#.*//' "${out_dir}/peaks.txt" > "${out_dir}/peaksTemp.bed"
        awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' "${out_dir}/peaksTemp.bed" > "${out_dir}/peaks.bed"
        rm "${out_dir}/peaksTemp.bed"
        sort -k1,1 -k2,2n "${out_dir}/peaks.bed" > "${out_dir}/peaks_Sorted.bed"
    fi

    log_info "Peak calling complete. Log saved to $log_file"
}

# ── CTK peak calling ────────────────────────────────────────────────────────
# 4. Peak Calling (CTK tag2peak.pl)
run_peak_calling_ctk() {
    local input_bed="$1"
    local out_dir="$2"
    local peak_dist="$3"
    local log_file="${out_dir}_ctk.log"

    update_status "Peaks"
    log_info "Calling peaks with CTK tag2peak.pl..."

    local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/tag2peak_cache.XXXXXX")
    local raw_peaks="${out_dir}_raw.bed"

    echo "Running CTK tag2peak.pl..." > "$log_file"
    $CONDA_PREFIX/bin/perl $(which tag2peak.pl) -big -ss --valley-seeking -minPH 2 -gap "${peak_dist}" \
        ${ADV_PEAK_CALLER_ARGS} -c "${cache_dir}" "${input_bed}" "${raw_peaks}" >> "$log_file" 2>&1
    local exit_code=$?
    rm -rf "$cache_dir"

    if [[ $exit_code -eq 0 && -s "$raw_peaks" ]]; then
        log_info "Peak calling complete. Log saved to $log_file"
    else
        log_error "CTK tag2peak.pl failed."
        rm -f "$raw_peaks"
        exit 1
    fi
}

# ── Peak calling dispatcher ─────────────────────────────────────────────────
# 4. Peak Calling - dispatcher
run_peak_calling() {
    if [[ "${PEAK_CALLER:-homer}" == "ctk" ]]; then
        run_peak_calling_ctk "$@"
    else
        run_peak_calling_homer "$@"
    fi
}

# ── Peak coverage matrix ────────────────────────────────────────────────────
# 4b. Add Enhanced Columns to Peak Coverage Matrix
# Adds: BC (groups), Raw Group Counts, Normalized Counts, BedGraph Stats
# Column Order: BC -> Raw Counts -> Normalized Counts -> BG Stats
add_matrix_columns() {
    local peak_matrix="$1"      # Path to peakCoverage.txt
    local peaks_bed="$2"        # Path to peaks_Sorted.bed
    local bg_dir="$3"           # BedGraph directory
    local scale_file="$4"       # Scale factors TSV
    local groups_file="$5"      # Optional groups file
    
    log_info "Adding enhanced columns to peak matrix..."
    
    # Validate inputs
    if [[ ! -f "$peak_matrix" ]]; then
        log_error "Peak matrix not found: $peak_matrix"
        return 1
    fi
    
    # Extract sample names from scale_factors.tsv (reliable source of actual samples)
    # NOT from matrix header (which may include CTK columns)
    if [[ ! -f "$scale_file" ]]; then
        log_warning "Scale factors file not found: $scale_file. Enhanced columns will be limited."
        local samples=()
    else
        local samples=($(cut -f1 "$scale_file"))
    fi
    
    log_info "  Samples detected: ${samples[*]}"
    
    # Prepare temp files for new columns
    local new_cols_file=$(mktemp)
    local new_header=""
    
    # -------------------------------------------
    # STEP 1: Biological Complexity (BC) - Groups Only
    # -------------------------------------------
    if [[ -n "$groups_file" && -f "$groups_file" ]]; then
        log_info "  Calculating Biological Complexity (BC)..."
        local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
        
        for group in $unique_groups; do
            log_info "    Group: $group"
            # Get sample column indices for this group
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file" | tr '\n' ' ')
            
            # For each peak (row), count samples with count > 0
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                # Build map of sample name -> column index (1-based, offset by 6)
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
            }
            NR==1 { print "BC_'"$group"'"; next }
            {
                bc=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map) {
                        c = col_map[s]
                        if($c + 0 > 0) bc++
                    }
                }
                print bc
            }
            ' "$peak_matrix" > "bc_${group}.col"
            
            # Add to new columns
            if [[ ! -s "$new_cols_file" ]]; then
                cat "bc_${group}.col" > "$new_cols_file"
            else
                paste "$new_cols_file" "bc_${group}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
            fi
            rm -f "bc_${group}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 2: Normalized Read Counts (Per Sample)
    # -------------------------------------------
    if [[ -f "$scale_file" ]]; then
        log_info "  Adding normalized read counts..."
        
        for sample in "${samples[@]}"; do
            # Get scale factor for this sample using awk for reliable tab-delimited matching
            local sf=$(awk -F'\t' -v s="$sample" '$1==s {print $3; exit}' "$scale_file")
            if [[ -z "$sf" ]]; then
                log_warning "    Scale factor not found for $sample, using 1.0"
                sf="1.0"
            fi
            
            # Get column index for this sample (1-based)
            local col_idx=0
            for i in "${!samples[@]}"; do
                if [[ "${samples[$i]}" == "$sample" ]]; then
                    col_idx=$((i + 7))  # Offset by 6 base columns + 1 for 1-indexing
                    break
                fi
            done
            
            # Calculate normalized value: raw_count * scale_factor
            awk -F'\t' -v col="$col_idx" -v sf="$sf" '
            NR==1 { print "NormedTC_'"${sample}"'"; next }
            { printf "%.4f\n", $col * sf }
            ' "$peak_matrix" > "normed_${sample}.col"
            
            # Add to new columns
            if [[ ! -s "$new_cols_file" ]]; then
                cat "normed_${sample}.col" > "$new_cols_file"
            else
                paste "$new_cols_file" "normed_${sample}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
            fi
            rm -f "normed_${sample}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 3: Group Columns (Raw Sum + Normalized Sum) - Groups Only
    # -------------------------------------------
    if [[ -n "$groups_file" && -f "$groups_file" ]]; then
        log_info "  Adding group aggregate columns..."
        local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
        
        for group in $unique_groups; do
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file" | tr '\n' ' ')
            
            # Sum raw counts for group
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
            }
            NR==1 { print "TC_'"$group"'"; next }
            {
                sum=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map) sum += $col_map[s]
                }
                print sum
            }
            ' "$peak_matrix" > "grp_raw_${group}.col"
            
            paste "$new_cols_file" "grp_raw_${group}.col" > "${new_cols_file}.tmp"
            mv "${new_cols_file}.tmp" "$new_cols_file"
            rm -f "grp_raw_${group}.col"
            
            # Sum normalized counts for group
            # This needs to reference the normalized columns we just added
            # For simplicity, recalculate using scale factors
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" -v sf_file="$scale_file" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
                # Load scale factors
                while((getline line < sf_file) > 0) {
                    split(line, sf_parts, "\t")
                    # Extract sample name from path
                    n = split(sf_parts[1], path_parts, "/")
                    sname = path_parts[n]
                    scale[sname] = sf_parts[3]
                }
            }
            NR==1 { print "NormedTC_'"$group"'"; next }
            {
                sum=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map && s in scale) {
                        sum += $col_map[s] * scale[s]
                    }
                }
                printf "%.4f\n", sum
            }
            ' "$peak_matrix" > "grp_normed_${group}.col"
            
            paste "$new_cols_file" "grp_normed_${group}.col" > "${new_cols_file}.tmp"
            mv "${new_cols_file}.tmp" "$new_cols_file"
            rm -f "grp_normed_${group}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 4: BedGraph Stats (Sum/Avg/Max) Per Sample
    # -------------------------------------------
    if [[ -d "$bg_dir" ]] && [[ ${#samples[@]} -gt 0 ]]; then
        log_info "  Adding BedGraph statistics..."
        
        # Split peaks by strand and sort for bedtools compatibility (same as STEP 5)
        awk -F'\t' '$6=="+"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_pos.tmp.bed
        awk -F'\t' '$6=="-"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_neg.tmp.bed
        
        for sample in "${samples[@]}"; do
            local bg_pos="${bg_dir}/${sample}_pos.bedgraph"
            local bg_neg="${bg_dir}/${sample}_neg.bedgraph"
            
            if [[ ! -f "$bg_pos" || ! -f "$bg_neg" ]]; then
                log_warning "    BedGraph not found for $sample, skipping."
                continue
            fi
            
            # Sort bedgraphs for bedtools compatibility (strip track header if present)
            grep -v "^track" "$bg_pos" | sort -k1,1 -k2,2n > "${sample}_pos_sorted.bg.tmp"
            grep -v "^track" "$bg_neg" | sort -k1,1 -k2,2n > "${sample}_neg_sorted.bg.tmp"
            
            for stat in sum mean max; do
                # Run bedtools map on sorted files (same as STEP 5)
                bedtools map -a peaks_pos.tmp.bed -b "${sample}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_pos_${stat}.tmp" 2>/dev/null
                bedtools map -a peaks_neg.tmp.bed -b "${sample}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_neg_${stat}.tmp" 2>/dev/null
                
                # Build lookup and match to original peak order (same logic as STEP 5)
                awk -F'\t' '
                BEGIN { 
                    while((getline < "bg_pos_'"$stat"'.tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                    while((getline < "bg_neg_'"$stat"'.tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                }
                NR==1 { print "Cov" toupper(substr("'"$stat"'",1,1)) substr("'"$stat"'",2) "_'"$sample"'"; next }
                {
                    key = $1"\t"$2"\t"$3
                    if($6 == "+") print (key in pos) ? pos[key] : 0
                    else print (key in neg) ? neg[key] : 0
                }
                ' "$peaks_bed" > "bg_${sample}_${stat}.col"
                
                paste "$new_cols_file" "bg_${sample}_${stat}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
                rm -f "bg_${sample}_${stat}.col"
            done
            
            rm -f "${sample}_pos_sorted.bg.tmp" "${sample}_neg_sorted.bg.tmp"
        done
        
        rm -f peaks_pos.tmp.bed peaks_neg.tmp.bed "bg_pos_"*.tmp "bg_neg_"*.tmp
        
        # -------------------------------------------
        # STEP 5: Group BedGraph Stats (from combined bedgraph) - Groups Only
        # -------------------------------------------
        if [[ -n "$groups_file" && -f "$groups_file" ]]; then
            log_info "  Adding group BedGraph statistics..."
            local combined_bg_dir="${bg_dir}/COMBINED_BEDGRAPH"
            local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
            
            # Re-split and sort peaks for bedtools compatibility
            awk -F'\t' '$6=="+"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_pos.tmp.bed
            awk -F'\t' '$6=="-"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_neg.tmp.bed
            
            for group in $unique_groups; do
                local grp_bg_pos="${combined_bg_dir}/${group}_combined_pos.bedgraph"
                local grp_bg_neg="${combined_bg_dir}/${group}_combined_neg.bedgraph"
                
                if [[ ! -f "$grp_bg_pos" || ! -f "$grp_bg_neg" ]]; then
                    log_warning "    Combined BedGraph not found for $group, skipping."
                    continue
                fi
                
                # Sort group bedgraphs for bedtools compatibility (strip track header, ensure TABs)
                grep -v "^track" "$grp_bg_pos" | tr ' ' '\t' | sort -k1,1 -k2,2n > "${group}_pos_sorted.bg.tmp"
                grep -v "^track" "$grp_bg_neg" | tr ' ' '\t' | sort -k1,1 -k2,2n > "${group}_neg_sorted.bg.tmp"
                
                for stat in sum mean max; do
                    bedtools map -a peaks_pos.tmp.bed -b "${group}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_pos_${stat}.tmp"
                    bedtools map -a peaks_neg.tmp.bed -b "${group}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_neg_${stat}.tmp"
                    
                    awk -F'\t' '
                    BEGIN { 
                        while((getline < "bg_pos_'"$stat"'.tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                        while((getline < "bg_neg_'"$stat"'.tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                    }
                    NR==1 { print "Cov" toupper(substr("'$stat'",1,1)) substr("'$stat'",2) "_'$group'"; next }
                    {
                        key = $1"\t"$2"\t"$3
                        if($6 == "+") print (key in pos) ? pos[key] : 0
                        else print (key in neg) ? neg[key] : 0
                    }
                    ' "$peaks_bed" > "bg_${group}_${stat}.col"
                    
                    paste "$new_cols_file" "bg_${group}_${stat}.col" > "${new_cols_file}.tmp"
                    mv "${new_cols_file}.tmp" "$new_cols_file"
                    rm -f "bg_${group}_${stat}.col"
                done
                
                rm -f "${group}_pos_sorted.bg.tmp" "${group}_neg_sorted.bg.tmp"
            done
            
            rm -f peaks_pos.tmp.bed peaks_neg.tmp.bed "bg_pos_"*.tmp "bg_neg_"*.tmp
        fi
    fi
    
    # -------------------------------------------
    # FINAL: Paste all new columns to matrix
    # -------------------------------------------
    if [[ -s "$new_cols_file" ]]; then
        paste "$peak_matrix" "$new_cols_file" > "${peak_matrix}.enhanced"
        mv "${peak_matrix}.enhanced" "$peak_matrix"
        log_info "Enhanced columns added to $peak_matrix"
    fi
    
    # -------------------------------------------
    # REORDER: Group columns by type (prefix)
    # Order: base -> TC_ -> NormedTC_ -> BC_ -> DEL_ -> SUB_ -> TRUNC_ -> BG*
    # -------------------------------------------
    log_info "Reordering columns by type..."
    awk -F'\t' '
    BEGIN { OFS="\t" }
    NR==1 {
        # Parse header and categorize columns by prefix
        for(i=1; i<=NF; i++) {
            h = $i
            if(h ~ /^(chr|start|end|name|score|strand)$/) { order[i] = 1; base[i] = h }
            else if(h ~ /^TC_/)      { order[i] = 2; tc[i] = h }
            else if(h ~ /^NormedTC_/) { order[i] = 3; ntc[i] = h }
            else if(h ~ /^BC_/)      { order[i] = 4; bc[i] = h }
            else if(h ~ /^DEL_/)     { order[i] = 5; del[i] = h }
            else if(h ~ /^SUB_/)     { order[i] = 6; sub_[i] = h }
            else if(h ~ /^TRUNC_/)   { order[i] = 7; trunc[i] = h }
            else if(h ~ /^BG/)       { order[i] = 8; bg[i] = h }
            else                     { order[i] = 9; other[i] = h }
            headers[i] = h
        }
        # Build column order
        n = 0
        for(o=1; o<=9; o++) {
            for(i=1; i<=NF; i++) {
                if(order[i] == o) { col_order[++n] = i }
            }
        }
        total_cols = n
        # Print reordered header
        for(j=1; j<=total_cols; j++) {
            printf "%s%s", headers[col_order[j]], (j<total_cols ? OFS : ORS)
        }
        next
    }
    {
        # Print reordered data
        for(j=1; j<=total_cols; j++) {
            printf "%s%s", $col_order[j], (j<total_cols ? OFS : ORS)
        }
    }
    ' "$peak_matrix" > "${peak_matrix}.reordered"
    mv "${peak_matrix}.reordered" "$peak_matrix"
    log_info "Columns reordered by type prefix"
    
    rm -f "$new_cols_file"
}

# ── CTK columns for peak matrix ─────────────────────────────────────────────
# 4b. Add CTK columns to peak coverage matrix
# Adds _del, _sub, _trunc columns for each sample/group with CTK outputs
add_ctk_columns_to_peak_matrix() {
    local peak_matrix="$1"         # Input/output: peak coverage matrix file
    local peaks_bed="$2"           # Sorted peaks BED file for bedtools
    local ctk_dir="$3"             # Directory containing CTK outputs
    local cims_fdr="${4:-0.05}"    # FDR threshold for CIMS filtering

    local cits_pval="${5:-0.05}"   # P-value threshold for CITS filtering
    local groups_file="$6"         # Optional: Groups file for aggregation
    
    log_info "Adding CTK site counts to peak matrix..."
    log_info "CTK directory: $ctk_dir"
    

    
    # Helper function to add column from CTK file
    add_ctk_column() {
        local ctk_file="$1"
        local column_name="$2"
        local ctk_type="$3"  # "cims" or "cits"
        local threshold="$4"
        
        if [[ ! -s "$ctk_file" ]]; then
            log_warning "CTK file not found or empty: $ctk_file"
            return 1
        fi
        
        local temp_filtered="${ctk_file}.filtered.bed"
        local temp_coverage="${ctk_file}.coverage.txt"
        
        # Filter by threshold and convert to BED
        if [[ "$ctk_type" == "cims" ]]; then
            # CIMS: Column 9 is FDR, skip header line
            grep -v "^#" "$ctk_file" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$temp_filtered"
        else
            # CITS: P-value is embedded in name as [P=value]
            grep -v "^#" "$ctk_file" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$temp_filtered"
        fi
        
        local site_count=$(wc -l < "$temp_filtered")
        log_info "  $column_name: $site_count significant sites"
        
        if [[ "$site_count" -eq 0 ]]; then
            # Add column of zeros
            local num_rows=$(($(wc -l < "$peak_matrix") - 1))  # Subtract header
            echo "$column_name" > temp_col.txt
            yes "0" | head -n "$num_rows" >> temp_col.txt
        else
            # Count sites per peak using bedtools
            bedtools coverage -s -a "$peaks_bed" -b "$temp_filtered" > "$temp_coverage"
            
            # Extract count column (7th) and combine with header
            echo "$column_name" > temp_col.txt
            awk '{print $7}' "$temp_coverage" >> temp_col.txt
        fi
        
        # Paste to matrix
        paste "$peak_matrix" temp_col.txt > temp_matrix.txt
        mv temp_matrix.txt "$peak_matrix"
        
        # Cleanup
        rm -f "$temp_filtered" "$temp_coverage" temp_col.txt
    }
    

    
    # -------------------------------------------------------------------------
    # Scenario A: Group Aggregation Mode (if groups_file provided)
    # -------------------------------------------------------------------------
    if [[ -n "$groups_file" ]]; then
        log_info "Using Groups File for CTK Aggregation..."
        
        # Parse groups
        local groups_map=$(mktemp)
        parse_groups_file "$groups_file" "$groups_map"
        local unique_groups=$(cut -f2 "$groups_map" | sort -u)
        
        for group in $unique_groups; do
            # Get samples for this group
            local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
            
            # --- Aggregate CIMS (Deletions) ---
            local group_del_bed="${ctk_dir}/${group}_aggregated_CIMS_del.txt"
            local precalc_del="${ctk_dir}/${group}/CIMS/${group}_CIMS_del.txt"
            if [[ -s "$precalc_del" ]]; then
                cp "$precalc_del" "$group_del_bed"
            else
                > "$group_del_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    # Look for sample CIMS file in CTK dir or subdirs
                    local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_del.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_del_bed"
                done
            fi
            if [[ -s "$group_del_bed" ]]; then
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            else
                # still add empty column for consistency?
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CIMS (Substitutions) ---
            local group_sub_bed="${ctk_dir}/${group}_aggregated_CIMS_sub.txt"
            local precalc_sub="${ctk_dir}/${group}/CIMS/${group}_CIMS_sub.txt"
             if [[ -s "$precalc_sub" ]]; then
                cp "$precalc_sub" "$group_sub_bed"
            else
                > "$group_sub_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_sub.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_sub_bed"
                done
            fi
            if [[ -s "$group_sub_bed" ]]; then
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            else
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CITS (Truncations) ---
            local group_cits_bed="${ctk_dir}/${group}_aggregated_CITS.txt"
            local precalc_cits="${ctk_dir}/${group}/CITS/${group}_CITS.bed"
            # Try .bed first, then .txt
            if [[ ! -s "$precalc_cits" ]]; then precalc_cits="${ctk_dir}/${group}/CITS/${group}_CITS.txt"; fi
            
             if [[ -s "$precalc_cits" ]]; then
                cp "$precalc_cits" "$group_cits_bed"
            else
                > "$group_cits_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    local s_file=$(find "$ctk_dir" -name "${sample}_CITS.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_cits_bed"
                done
            fi
            if [[ -s "$group_cits_bed" ]]; then
                add_ctk_column "$group_cits_bed" "TRUNC_${group}" "cits" "$cits_pval"
            else
                add_ctk_column "$group_cits_bed" "TRUNC_${group}" "cits" "$cits_pval"
            fi
            
            # Cleanup aggregated files
            rm -f "$group_del_bed" "$group_sub_bed" "$group_cits_bed"
        done
        
        rm "$groups_map"
        return 0
    fi
    
    # -------------------------------------------------------------------------
    # Scenario B: Default / Legacy Mode (All found in dir)
    # -------------------------------------------------------------------------
    # Find and process CIMS deletion files
    # Check both: 1) Direct structure: ctk_dir/CIMS/  2) Group structure: ctk_dir/*/CIMS/
    local cims_found=false
    
    # First try direct structure
    if [[ -d "${ctk_dir}/CIMS" ]]; then
        cims_found=true
        for cims_del_file in "${ctk_dir}/CIMS/"*_CIMS_del.txt; do
            if [[ -f "$cims_del_file" ]]; then
                local name=$(basename "$cims_del_file" _CIMS_del.txt)
                add_ctk_column "$cims_del_file" "DEL_${name}" "cims" "$cims_fdr"
            fi
        done
        
        for cims_sub_file in "${ctk_dir}/CIMS/"*_CIMS_sub.txt; do
            if [[ -f "$cims_sub_file" ]]; then
                local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                add_ctk_column "$cims_sub_file" "SUB_${name}" "cims" "$cims_fdr"
            fi
        done
    fi
    
    # If no direct CIMS folder, try group subfolders (ctk_dir/*/CIMS/)
    if [[ "$cims_found" == "false" ]]; then
        for group_dir in "${ctk_dir}"/*/; do
            if [[ -d "${group_dir}CIMS" ]]; then
                cims_found=true
                for cims_del_file in "${group_dir}CIMS/"*_CIMS_del.txt; do
                    if [[ -f "$cims_del_file" ]]; then
                        local name=$(basename "$cims_del_file" _CIMS_del.txt)
                        add_ctk_column "$cims_del_file" "DEL_${name}" "cims" "$cims_fdr"
                    fi
                done
                
                for cims_sub_file in "${group_dir}CIMS/"*_CIMS_sub.txt; do
                    if [[ -f "$cims_sub_file" ]]; then
                        local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                        add_ctk_column "$cims_sub_file" "SUB_${name}" "cims" "$cims_fdr"
                    fi
                done
            fi
        done
    fi
    
    # Find and process CITS files
    # Check both: 1) Direct structure: ctk_dir/CITS/  2) Group structure: ctk_dir/*/CITS/
    local cits_found=false
    
    # First try direct structure
    if [[ -d "${ctk_dir}/CITS" ]]; then
        cits_found=true
        for cits_file in "${ctk_dir}/CITS/"*_CITS.txt "${ctk_dir}/CITS/"*_CITS.bed; do
            if [[ -f "$cits_file" ]]; then
                local name=$(basename "$cits_file" | sed 's/_CITS\.\(txt\|bed\)$//')
                add_ctk_column "$cits_file" "TRUNC_${name}" "cits" "$cits_pval"
            fi
        done
    fi
    
    # If no direct CITS folder, try group subfolders
    if [[ "$cits_found" == "false" ]]; then
        for group_dir in "${ctk_dir}"/*/; do
            if [[ -d "${group_dir}CITS" ]]; then
                cits_found=true
                for cits_file in "${group_dir}CITS/"*_CITS.txt "${group_dir}CITS/"*_CITS.bed; do
                    if [[ -f "$cits_file" ]]; then
                        local name=$(basename "$cits_file" | sed 's/_CITS\.\(txt\|bed\)$//')
                        add_ctk_column "$cits_file" "TRUNC_${name}" "cits" "$cits_pval"
                    fi
                done
            fi
        done
    fi
    
    log_info "CTK columns added to peak matrix."
}

# ── Coverage / BedGraph ─────────────────────────────────────────────────────
# 7. Coverage Analysis (Bedgraph)
# 7. Coverage Analysis (Bedgraph)
run_coverage() {
    local input_bed="$1"      
    local output_prefix="$2"
    local genome_file="$3"     
    local bam_file="$4"       # [NEW] Explicit BAM path

    update_status "Bedgraph"
    log_info "generating normalized, filtered bedGraphs..."
    
    # Validation
    if [[ ! -f "$bam_file" ]]; then
         log_error "Combined/Normalized BedGraph requires BAM input, but could not locate: $bam_file"
         return 1
    fi
    
    # 2. Calculate Scale Factor for Normalization
    # Count mapped reads (primary alignments only? or all mapped?)
    # usually -F 4 (mapped). User script used -F 4.
    local mapped=$(samtools view -c -F 4 "$bam_file")
    
    if [[ "$mapped" -eq 0 ]]; then
        log_warning "No mapped reads found in $bam_file. Skipping BedGraph."
        return 0
    fi
    
    local scale=$(echo "scale=6; 1000000 / $mapped" | bc)
    log_info "  Mapped Reads: $mapped | Scale Factor: $scale"
    
    # Store scale factor for later use (normalized counts)
    # Use same directory as bedgraph output (output_prefix parent) to avoid OUTPUT_ROOT scoping issues
    local bg_dir=$(dirname "$output_prefix")
    local scale_file="${bg_dir}/scale_factors.tsv"
    echo -e "$(basename "$output_prefix")\t${mapped}\t${scale}" >> "$scale_file"
    
    # 3. Generate BedGraph (Filtered & Normalized)
    # Using 'cigar !~ "N"' to remove junction reads
    # Streaming samtools -> bedtools genomecov
    
    # Extract sample name from output_prefix for track header
    local sample_name
    sample_name=$(basename "$output_prefix")

    # Positive Strand: generate sorted bedgraph then prepend track header
    local pos_bg="${output_prefix}_pos.bedgraph"
    local pos_tmp="${output_prefix}_pos.bedgraph.tmp"
    samtools view -h -e 'cigar !~ "N"' "$bam_file" | \
    bedtools genomecov -ibam stdin -bg -strand + -scale "$scale" | \
    sort -k1,1 -k2,2n > "$pos_tmp"
    echo "track type=bedGraph name=\"${sample_name}\" description=\"Positive Strand\"" > "$pos_bg"
    cat "$pos_tmp" >> "$pos_bg"
    rm -f "$pos_tmp"

    # Negative Strand: generate sorted bedgraph then prepend track header
    local neg_bg="${output_prefix}_neg.bedgraph"
    local neg_tmp="${output_prefix}_neg.bedgraph.tmp"
    samtools view -h -e 'cigar !~ "N"' "$bam_file" | \
    bedtools genomecov -ibam stdin -bg -strand - -scale "$scale" | \
    sort -k1,1 -k2,2n > "$neg_tmp"
    echo "track type=bedGraph name=\"${sample_name}\" description=\"Negative Strand\"" > "$neg_bg"
    cat "$neg_tmp" >> "$neg_bg"
    rm -f "$neg_tmp"

    log_info "Bedgraphs generated: ${output_prefix}_pos.bedgraph, ${output_prefix}_neg.bedgraph"
}

# ── Combined BedGraph (group averaging) ─────────────────────────────────────
# 8. Combined BedGraph Generation (Group Averaging)
run_combined_bedgraph() {
    local output_dir="$1"
    local groups_file="$2"
    local bedgraph_dir="$3"
    
    log_info "Generating combined average bedgraphs..."
    
    mkdir -p "$bedgraph_dir/COMBINED_BEDGRAPH"
    
    # Identify Groups
    # If groups file provided, use it. Else, maybe regex?
    # For robust implementation, we'll assume GROUPS_FILE structure: SampleName<TAB>GroupName
    
    if [[ -z "$groups_file" || ! -f "$groups_file" ]]; then
        log_warning "No valid groups file provided. Skipping combined bedgraph generation."
        return 0
    fi
    
    # Extract unique groups
    local groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort | uniq)
    
    for group in $groups; do
        log_info "Processing Group: $group"
        
        # Get samples for this group
        local samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file")
        
        # Process Positive and Negative strands separately
        for strand in "pos" "neg"; do
            local bg_files=""
            local count=0
            
            for sample in $samples; do
                # Construct expected filename from run_coverage output
                # {sample}_pos.bedgraph
                local f="$bedgraph_dir/${sample}_${strand}.bedgraph"
                if [[ -f "$f" ]]; then
                    bg_files="$bg_files $f"
                    ((count++))
                fi
            done
            
            if [[ $count -eq 1 ]]; then
                # Single sample in group — copy directly (no union needed)
                local output_file="$bedgraph_dir/COMBINED_BEDGRAPH/${group}_combined_${strand}.bedgraph"
                local strand_desc
                if [[ "$strand" == "pos" ]]; then strand_desc="Combined Positive Strand"; else strand_desc="Combined Negative Strand"; fi
                echo "track type=bedGraph name=\"${group}\" description=\"${strand_desc}\"" > "$output_file"
                grep -v "^track" $bg_files >> "$output_file"
                log_info "  Generated: $(basename "$output_file") (1 sample, copied directly)"
            elif [[ $count -gt 1 ]]; then
                # Union and Average
                # Unionbedg produces: chrom start end val1 val2 ... valN
                # Column 1,2,3 are coords. Columns 4 to 3+N are values.
                # Average = sum(col 4..NF) / (NF-3)
                
                local output_file="$bedgraph_dir/COMBINED_BEDGRAPH/${group}_combined_${strand}.bedgraph"
                local strand_desc
                if [[ "$strand" == "pos" ]]; then
                    strand_desc="Combined Positive Strand"
                else
                    strand_desc="Combined Negative Strand"
                fi
                local combined_tmp="${output_file}.tmp"

                # unionbedg cannot read multiple files from stdin — strip track
                # headers to individual temp files and pass as explicit args
                local tmp_dir
                tmp_dir=$(mktemp -d)
                local tmp_files=()
                for bg in $bg_files; do
                    local tmp_bg="${tmp_dir}/$(basename "$bg")"
                    grep -v "^track" "$bg" > "$tmp_bg"
                    tmp_files+=("$tmp_bg")
                done

                bedtools unionbedg -i "${tmp_files[@]}" | \
                awk -v N="$count" 'BEGIN{OFS="\t"} {sum=0; for(i=4;i<=NF;i++) sum+=$i; print $1,$2,$3,sum/N}' | \
                sort -k1,1 -k2,2n > "$combined_tmp"
                echo "track type=bedGraph name=\"${group}\" description=\"${strand_desc}\"" > "$output_file"
                cat "$combined_tmp" >> "$output_file"
                rm -f "$combined_tmp"
                rm -rf "$tmp_dir"

                log_info "  Generated: $(basename "$output_file") ($count replicates)"
            else
                log_warning "  No bedgraph files found for group $group ($strand)"
            fi
        done
    done
}
