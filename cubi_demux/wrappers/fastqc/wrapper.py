# -*- coding: utf-8 -*-
"""Snakemake wrapper for FastQC.
"""

from snakemake import shell

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

shell.executable('/bin/bash')

shell(r"""
fastqc \
    -o $(dirname {snakemake.output.html}) \
    {snakemake.input.fastq}
""")
