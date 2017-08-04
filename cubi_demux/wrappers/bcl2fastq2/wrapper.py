# -*- coding: utf-8 -*-
"""Snakemake wrapper for bcl2fastq2.
"""

from snakemake import shell

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

shell.executable('/bin/bash')

# Get number of barcode mismatches, defaults to 1 for bcl2fastq v2.
barcode_mismatches = snakemake.config['cubi_demux'].get('barcode_mismatches')
if barcode_mismatches is None:
    barcode_mismatches = 1

# More than 8 threads will not work for bcl2fastq.
bcl2fastq_threads = min(8, snakemake.config['cubi_demux']['cores'])

shell(r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# -------------------------------------------------------------------------------------------------
# Print Sample Sheet

head -n 10000 {snakemake.input.sheet}

# -------------------------------------------------------------------------------------------------
# Execute bcl2fastq v2

bcl2fastq \
    --barcode-mismatches {barcode_mismatches} \
    --sample-sheet {snakemake.input.sheet} \
    --runfolder-dir {snakemake.params.input_dir} \
    --output-dir $TMPDIR/demux_out \
    --demultiplexing-threads {bcl2fastq_threads} \
    --processing-threads {bcl2fastq_threads} \
    {snakemake.params.tiles_arg}

# -------------------------------------------------------------------------------------------------
# Move Files to Destination.

# Move sample FASTQ files.
flowcell={snakemake.params.flowcell_token}
srcdir=$TMPDIR/demux_out/Project
for path in $srcdir/*; do
    sample=$(basename $path | rev | cut -d _ -f 5- | rev)
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/$sample/$flowcell/$lane/$(basename $path)

    mv $path $dest
    pushd $(dirname $dest) \
    && md5sum $(basename $dest) >$(basename $dest).md5 \
    && popd
done

# Move undetermined FASTQ files.
srcdir=$TMPDIR/demux_out
for path in $srcdir/Undetermined_*; do
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/Undetermined/$flowcell/$lane/$(basename $path)

    mv $path $dest
    pushd $(dirname $dest) \
    && md5sum $(basename $dest) >$(basename $dest).md5 \
    && popd
done
""")
