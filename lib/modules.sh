#!/bin/bash
# lib/modules.sh — compatibility shim (CLIPittyClip v3.4)
#
# In v3.4 the monolithic modules.sh was split into four domain files:
#   align.sh   — demux, eCLIP, fastp, STAR, Bowtie2, parseAlignment
#   peaks.sh   — HOMER/CTK peak calling, coverage, bedgraph, matrix
#   ctk.sh     — tag2collapse, CIMS, CITS, group CTK
#   clink.sh   — collapse.py, pileup, Clink CITS/CIMS, group Clink
#
# Sourcing this file continues to work unchanged for all existing callers.

_lib="$(dirname "${BASH_SOURCE[0]}")"
source "$_lib/align.sh"
source "$_lib/peaks.sh"
source "$_lib/ctk.sh"
source "$_lib/clink.sh"
