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

# More than 8 threads will not work for bcl2fastq.
bcl2fastq_threads = min(8, snakemake.config['cubi_demux']['cores'])


shell(r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

#export TMPDIR=$(mktemp -d)
#trap "rm -rf $TMPDIR" EXIT
export TMPDIR=$TMPDIR/v1yayaya
mkdir -p $TMPDIR

# -------------------------------------------------------------------------------------------------
# Print Sample Sheet

head -n 10000 {snakemake.input.sheet}

# -------------------------------------------------------------------------------------------------
# Run blc2fastq v1

# Prepare output directory with Makefile etc.

configureBclToFastq.pl \
    --mismatches {barcode_mismatches} \
    --sample-sheet {snakemake.input.sheet} \
    --input-dir {snakemake.params.input_dir}/Data/Intensities/BaseCalls \
    --output-dir $TMPDIR/demux_out \
    --fastq-cluster-count 0 \
    --force \
    {snakemake.params.tiles_arg}

# Actually perform the demultiplexing using Make
make \
    -C $TMPDIR/demux_out \
    -j {bcl2fastq_threads}

# -------------------------------------------------------------------------------------------------
# Move Files to Destination.

# Move sample FASTQ files.
flowcell={snakemake.params.flowcell_token}
srcdir=$TMPDIR/demux_out/Project_Project
for path in $srcdir/*/*; do
    sample=$(basename $(dirname $path) | cut -d _ -f 2-)
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/$sample/$flowcell/$lane/$(basename $path)

    cp $path $dest
    pushd $(dirname $dest) \
    && md5sum $(basename $dest) >($basename $dest).md5 \
    && popd
done

# Move undetermined FASTQ files.
srcdir=$TMPDIR/demux_out
for path in $srcdir/Undetermined_indices/*/*; do
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/Undetermined/$flowcell/$lane/$(basename $path)

    cp $path $dest
    pushd $(dirname $dest) \
    && md5sum $(basename $dest) >($basename $dest).md5 \
    && popd
done
""")
