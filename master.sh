#!/bin/bash
ortholog_name=og_3922  #no file extensions

python3 ./scripts/orthologs.py ortho_subset.tsv seq_files.tsv ${ortholog_name} 4

#arguments for above; ortho_subset, seq_files, ortholog name, minimum species

./scripts/auto_align.sh ./results/unaligned/${ortholog_name}_unaligned.fna


