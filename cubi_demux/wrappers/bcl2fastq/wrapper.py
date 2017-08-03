# -*- coding: utf-8 -*-
"""Snakemake wrapper for bcl2fastq2.
"""

from snakemake import shell

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

shell.executable('/bin/bash')

# Get number of barcode mismatches, defaults to 0 for bcl2fastq v1.
barcode_mismatches = snakemake.config['cubi_demux'].get('barcode_mismatches')
if barcode_mismatches is None:
    barcode_mismatches = 0

shell(r"""
set -euo pipefail
IFS="\n\t"

# Print sample sheet back to the user.
head -n 10000 {snakemake.input.sheet}

# Prepare the output directory with Makefile etc.
configureBclToFastq.pl \
    --mismatches {barcode_mismatches} \
    --sample-sheet {sample_sheet} \
    --input-dir {args.input_dir}/Data/Intensities/BaseCalls \
    --output-dir {snakemake.params.output_dir} \
    --fastq-cluster-count 0 \
    {tiles_arg}

# Actually perform the demultiplexing using Make
make \
    -C {snakemake.params.output_dir} \
    -j {snakemake.config[cubi_demux][num_threads]}
""")
