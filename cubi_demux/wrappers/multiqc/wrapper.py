# -*- coding: utf-8 -*-
"""Snakemake wrapper for FastQC.
"""

from snakemake import shell

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

shell.executable('/bin/bash')

shell(r"""
multiqc \
    --outdir $(dirname {snakemake.output.html}) \
    --interactive \
    {snakemake.input}
""")
