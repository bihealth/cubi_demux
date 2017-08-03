# -*- coding: utf-8 -*-
"""Snakemake wrapper for HTS Screen utility.
"""

import os

from snakemake import shell

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

shell.executable('/bin/bash')

match_vector_py = os.path.join(os.path.dirname(__file__), 'match_vector_to_report.py')

introduction = r"""
set -euo pipefail

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
"""

functions = r"""
# -------------------------------------------------------------------------------------------------
# BWA Alignment Helper.

bwa_alignment()
{{
    num=$1
    name=$2
    index=$3

    set -euo pipefail

    # Align using BWA-MEM and remove headers from SAM.  Then,
    # filter to primary alignments and print out TSV format with
    # the following columns, sorted by read name.
    #
    # 1: read name
    # 2: number of alignments found
    # 3: sequence (if $num==1)
    # 4: PHRED qualities (if $num==1)
    bwa mem \
        -t 2 \
        $index \
        - \
    | tee $TMPDIR/bwa_$num.sam \
    | grep -v '^@' \
    | awk -v num=$num -F $'\t' '
        BEGIN {{ OFS=FS;}}
        (!and($2, 256) && !and($2, 1024) && !and($2, 2048)) {{
            if (and($2, 4)) {{
                count = 0;
            }} else {{
                count = 1;
                for (i = 12; i <= NF; ++i) {{
                    if ($i ~ "XA:Z") {{
                        s = $i;
                        gsub(/[^;]/, "", s);
                        count = count + length(s);
                    }}
                }}
                $2 = count;
            }}

            if (num == 1) {{
                print $1, count, $10, $11;
            }} else {{
                print $1, count;
            }}
        }}' \
    | sort -k1,1 \
    > $TMPDIR/bwa_$num.tsv
    touch $TMPDIR/bwa_$num.tsv.done
}}

# -------------------------------------------------------------------------------------------------
# BWA Postprocessing Helper.

process_bwa_results()
{{
    set -euo pipefail

    header="#name"
    paste_cmd="paste"
    i=1
    for ref in $*; do
        num=$i
        name=$(echo $ref | cut -d : -f 1)
        index=$(echo $ref | cut -d : -f 2)

        while [[ ! -f $TMPDIR/bwa_$num.tsv.done ]]; do
            sleep 2
        done

        header+="\t$name"

        if [[ $num -eq 1 ]]; then
            paste_cmd+=" $TMPDIR/bwa_$num.tsv"
        else
            paste_cmd+=" <(cut -f 2 $TMPDIR/bwa_$num.tsv)"
        fi

        let "i=$i+1"
    done
    header+="\tUnaligned"

    echo -e "$header" > $TMPDIR/match_vector.tsv
    eval "$paste_cmd" \
    | tee >(
        awk -F $'\t' '{{
                x = 1;
                for (i = 5; i <= NF; ++i) {{
                    if ($i >= 1) {{ x = 0; }}
                }}
                if (x) {{
                    printf("@%s\n%s\n+\n%s\n", $1, $3, $4);
                }}
            }}' \
        | gzip -c \
        > $TMPDIR/unaligned_reads.fq.gz
        ) \
    | cut -f 1,2,5- \
    | awk -F $'\t' '
        BEGIN {{ OFS=FS }}
        {{
            x = 1;
            for (i = 2; i <= NF; ++i) {{
                if ($i >= 1) {{ x = 0; }}
            }}
            print $0, x;
        }}' \
    >> $TMPDIR/match_vector.tsv
}}

# -------------------------------------------------------------------------------------------------
# Match Vector Postprocessing.

postprocess_match_vector()
{{
    set -euo pipefail
    fname=$1

    num_reads=$(( $(cat $fname | wc -l) - 1))

    hit_no_genomes=$(rev $fname | cut -f 1 | grep -v '^0$' | wc -l)
    perc_hit_no_genomes=$(printf "%0.2f" $(echo $hit_no_genomes / $num_reads | bc -l))

    genomes=$(head -n 1 $fname | rev | cut -f 2- | rev | cut -f 2-)

    python {match_vector_py} --in-file $fname
}}

# -------------------------------------------------------------------------------------------------
# HTS Screen Processing.
"""

# Build main main() code.
tee_cmd = 'tee'
refs = []
for i, item in enumerate(snakemake.config['hts_screen']['references']):
    tee_cmd += ' >(bwa_alignment {i} {ref_name} {path_index})'.format(
        i=(i + 1),
        ref_name=item['name'].replace(' ', '_'),
        path_index=item['bwa_index'])
    refs.append('{ref_name}:{path_index}'.format(
        ref_name=item['name'].replace(' ', '_'),
        path_index=item['bwa_index']))

generated_main = r"""
# Perform BWA Alignment -----------------------------------------
(seqtk sample {{snakemake.input.fastq}} {sample_rate} || true) \
| eval {tee_cmd} >/dev/null

# Process BWA Results -------------------------------------------
process_bwa_results {refs}

# Create Match Vector -------------------------------------------
prefix=$(dirname {{snakemake.output.txt}})

mkdir -p $prefix
gzip -c $TMPDIR/match_vector.tsv \
> $prefix/{fastq_basename}_match_vector.tsv.gz

postprocess_match_vector $TMPDIR/match_vector.tsv \
> $prefix/{fastq_basename}_hts_screen.txt
""".format(
    sample_rate=snakemake.config['hts_screen']['sample_rate'],
    tee_cmd=tee_cmd,
    refs=' '.join(refs),
    fastq_basename=os.path.basename(
        str(snakemake.input.fastq)[:-len('.fastq.gz')]))

kraken_postprocessing = r"""
# -------------------------------------------------------------------------------------------------
# Kraken Postprocessing

prefix=$(dirname {{snakemake.output.kraken}})
mkdir -p $prefix

if [[ $((zcat $TMPDIR/unaligned_reads.fq.gz || true) | head | wc -l) -eq 0 ]]; then
    >&2 echo "Not executing Kraken; no unaligned reads"
    touch $prefix/NO_UNALIGNED_READS
    test -f $prefix/{fastq_basename}_kraken.gz \
    || gzip -c /dev/null > $prefix/{fastq_basename}_kraken.gz
    test -f $prefix/{fastq_basename}_kraken_report.gz \
    || gzip -c /dev/null > $prefix/{fastq_basename}_kraken_report.gz
else
    kraken \
        --threads 2 \
        --db {{snakemake.config[hts_screen][kraken_db]}} \
        --fastq-input \
        --gzip-compressed \
        $TMPDIR/unaligned_reads.fq.gz \
    | gzip -c \
    > $prefix/{fastq_basename}_kraken.gz

    kraken-report \
        --db {{snakemake.config[hts_screen][kraken_db]}} \
        <(zcat $prefix/{fastq_basename}_kraken.gz) \
    | gzip -c \
    > $prefix/{fastq_basename}_kraken_report.gz
fi
""".format(fastq_basename=os.path.basename(
    str(snakemake.input.fastq)[:-len('.fastq.gz')]))

shell(
    introduction +
    functions +
    generated_main +
    kraken_postprocessing
)
