#!/bin/bash

# lib/utils.sh - Utility functions for CLIPittyClip

# define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get terminal width for separator
if [ -t 1 ]; then
    terminal_width=$(tput cols)
else
    terminal_width=80
fi
separator_line=$(printf "%${terminal_width}s" | tr ' ' '*')

log_info() {
    # Stick to FILE ONLY for info logs to keep console clean
    echo -e "[INFO] $1" >> "${LOG_FILE:-/dev/null}"
}

log_warning() {
    # Warnings go to file and Console (in Yellow)
    echo -e "[WARNING] $1" >> "${LOG_FILE}"
    echo -e "${YELLOW}[WARNING] $1${NC}" >&2
}

log_error() {
    # Errors go to file and Console (in Red)
    echo -e "[ERROR] $1" >> "${LOG_FILE}"
    echo -e "${RED}[ERROR] $1${NC}" >&2
}

# New Function: explicit console message (always visible)
console_msg() {
    echo -e "$1"
    # Also log it simply
    echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g' >> "${LOG_FILE}"
}

# New Function: Update status on the same line (for progress bars)
# Used for middle steps: prints "Status > " format
update_status() {
    local msg="$1"
    # Print status trail (e.g. "Mapping > ") without newline
    echo -ne "${msg} > " 
    
    # Log it as an event
    echo -e "[STATUS] $msg" >> "${LOG_FILE}"
}

# First status in a chain (no leading >)
update_status_first() {
    local msg="$1"
    echo -ne "${msg} > "
    echo -e "[STATUS] $msg" >> "${LOG_FILE}"
}

# Final status with newline
update_status_done() {
    echo -e "Done!"
    echo -e "[STATUS] Done" >> "${LOG_FILE}"
}

# Section header for dedup/demux style output (with indented items)
print_section_item() {
    local msg="$1"
    echo -e "  > ${msg}"
    echo -e "[SECTION] $msg" >> "${LOG_FILE}"
}

check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Dependency '$1' not found. Please install it or activate the environment."
        exit 1
    fi
}

check_file() {
    if [ ! -f "$1" ]; then
        log_error "File not found: $1"
        return 1
    fi
    return 0
}

check_star_index() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        log_error "Genome index directory not found: $dir"
        return 1
    fi
    # Check for critical STAR index files
    local missing=0
    for file in Genome SA SAindex genomeParameters.txt; do
        if [ ! -f "$dir/$file" ]; then
            log_error "Missing STAR index file: $dir/$file"
            missing=1
        fi
    done
    
    if [ "$missing" -eq 1 ]; then
        log_error "The specified directory does not appear to contain a valid STAR index."
        log_error "Please generate one using 'STAR --runMode genomeGenerate ...'"
        return 1
    fi
    return 0
}

check_bowtie_index() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        log_error "Genome index directory not found: $dir"
        return 1
    fi
    # Check for Bowtie2 index files (.1.bt2 or .1.bt2l)
    local found=$(find "$dir" -name "*.1.bt2" -o -name "*.1.bt2l" 2>/dev/null | head -n 1)
    
    if [[ -z "$found" ]]; then
        log_error "No Bowtie2 index files found in $dir."
        return 1
    fi
    return 0
}

# Check for ncRNA Bowtie2 index in annotation directory
# Returns 0 if found, 1 if not found (silent - caller handles messaging)
check_ncrna_index() {
    local dir="$1"
    # Check ncRNA/ subfolder first (recommended structure)
    if [[ -f "${dir}/ncRNA/ncrna.1.bt2" ]] || [[ -f "${dir}/ncRNA/ncrna.1.bt2l" ]]; then
        echo "${dir}/ncRNA"  # Return the path where index was found
        return 0
    fi
    # Fallback: check top-level directory (backwards compatibility)
    if [[ -f "${dir}/ncrna.1.bt2" ]] || [[ -f "${dir}/ncrna.1.bt2l" ]]; then
        echo "${dir}"  # Return the path where index was found
        return 0
    fi
    return 1
}

show_header() {
    echo "$separator_line"
    echo -e "${BLUE}CLIPittyClip: Modern CLIP-seq Analysis Pipeline${NC}"
    echo "$separator_line"
    echo "Version 3.1.0"
    echo "Author: Soon Yi (Updated by Antigravity)"
    echo "Last updated: $(date +'%Y-%m-%d')"
    echo "$separator_line"
    # Log the header too
    {
        echo "$separator_line"
        echo "CLIPittyClip: Modern CLIP-seq Analysis Pipeline"
        echo "Version 3.1.0"
        echo "$separator_line"
    } >> "${LOG_FILE}"
}

execute_cmd() {
    local cmd="$1"
    echo "[EXEC] $cmd" >> "${LOG_FILE}"
    
    local exit_code=0
    if [[ "$VERBOSE" == "true" ]]; then
        # Verbose: Print to stdout and log file
        eval "$cmd" 2>&1 | tee -a "${LOG_FILE}"
        exit_code=${PIPESTATUS[0]}
    else
        # Quiet: Redirect everything to log file
        eval "$cmd" >> "${LOG_FILE}" 2>&1
        exit_code=$?
    fi
    return $exit_code
}

send_notification() {
    local title="$1"
    local message="$2"
    
    if [[ "$NOTIFY_MODE" == "true" ]]; then
        # 1. Sound (Terminal Bell - Universal)
        # We print to stderr to ensure it passes through pipelines if needed, but stdout is fine here
        echo -e "\a"
        
        # 2. Pop-up (macOS specific)
        if command -v osascript &> /dev/null; then
            # Escape quotes in message/title to avoid breaking osascript
            # Simple approach: just use single quotes for osascript command
            osascript -e "display notification \"$message\" with title \"$title\" sound name \"Glass\"" 2>/dev/null
        fi
        
        # 3. Linux (notify-send) - Future proofing
        if command -v notify-send &> /dev/null; then
             notify-send "$title" "$message" 2>/dev/null
        fi
    fi
}

# ── Parallel job helpers (moved from modules.sh v3.4) ────────────────────────

# Check if GNU parallel is installed
has_gnu_parallel() {
    command -v parallel &>/dev/null && parallel --version 2>&1 | head -1 | grep -q "GNU parallel"
}
calculate_optimal_parallel_jobs() {
    local user_threads="$1"
    local input_file="$2"
    local min_ram_per_job="${3:-2}"  # Default: 2GB per job minimum
    
    # Get available RAM in GB (use 'available' column from free)
    local available_ram_gb
    available_ram_gb=$(free -g 2>/dev/null | awk '/^Mem:/ {print $7}')
    
    # Fallback if free -g fails
    if [[ -z "$available_ram_gb" ]] || [[ "$available_ram_gb" -eq 0 ]]; then
        available_ram_gb=8  # Conservative default
    fi
    
    # Get file size and estimate RAM per job
    # Rule of thumb: jobs processing large files need more RAM
    local file_lines=0
    if [[ -f "$input_file" ]]; then
        file_lines=$(wc -l < "$input_file" 2>/dev/null || echo 0)
    fi
    
    # Estimate: for files with >1M lines, use more RAM per job
    local ram_per_job=$min_ram_per_job
    if [[ "$file_lines" -gt 1000000 ]]; then
        ram_per_job=4  # Large files need 4GB per job
    elif [[ "$file_lines" -gt 5000000 ]]; then
        ram_per_job=6  # Very large files need 6GB per job
    fi
    
    # Calculate RAM-based job limit
    local ram_based_jobs=$((available_ram_gb / ram_per_job))
    ram_based_jobs=$((ram_based_jobs > 0 ? ram_based_jobs : 1))  # At least 1
    
    # Final = minimum of user threads and RAM-based limit
    local optimal_jobs=$((user_threads < ram_based_jobs ? user_threads : ram_based_jobs))
    optimal_jobs=$((optimal_jobs > 0 ? optimal_jobs : 1))  # At least 1
    
    echo "$optimal_jobs"
}

# parse_groups_file - Parse groups.txt and output sample→group mapping
# Input: groups_file path, output_map temp file path
# Format: sample_name<TAB>group_name (lines starting with # are comments)
parse_groups_file() {
    local groups_file="$1"
    local output_map="$2"
    
    log_info "Parsing groups file: $groups_file"
    
    > "$output_map"  # Clear output file
    
    while IFS=$'\t' read -r sample group || [[ -n "$sample" ]]; do
        # Skip comments and empty lines
        [[ "$sample" =~ ^#.*$ ]] && continue
        [[ -z "$sample" ]] && continue
        
        # Strip common extensions
        sample="${sample%.fastq.gz}"
        sample="${sample%.fq.gz}"
        sample="${sample%.fastq}"
        sample="${sample%.fq}"
        
        # Write mapping
        echo -e "${sample}\t${group}" >> "$output_map"
    done < "$groups_file"
    
    local group_count=$(cut -f2 "$output_map" | sort -u | wc -l | tr -d ' ')
    local sample_count=$(wc -l < "$output_map" | tr -d ' ')
    log_info "Parsed $sample_count samples into $group_count groups"
}
