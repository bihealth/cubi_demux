# -*- coding: utf-8 -*-
"""Main entry point of the ``cubi-demux`` application."""

import os

from ruamel.yaml import YAML

from .utils import listify
from .bcl_to_fastq import Bcl2FastqV1Wrapper, Bcl2FastqV2Wrapper


def wrapper_path(path):
    """Generate path to wrapper"""
    return 'file://' + os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'wrappers', path))


class MissingConfiguration(Exception):
    """Raised on missing configuration."""


class InvalidConfiguration(Exception):
    """Raised on invalid configuration."""


class DemuxApp:
    """Application class."""

    def __init__(self, config):
        #: Configuration as loaded by Snakemake, possibly override with
        #: settings from command line.
        self.config = config
        # Check configuration.
        self._check_config()
        #: The wrapper for branching between RTA/bcl2fastq v1 and v2
        #: demultiplexing.
        self.demux_wrapper = self._build_wrapper()

    def _check_config(self):
        """Check configuration."""
        if not self.config['sample_sheet']:
            raise MissingConfiguration(
                'Configuration setting "sample_sheet" is missing from the '
                'YAML configuration file and/or --sample-sheet is missing '
                'from the command line.')
        if not self.config['cubi_demux']['input_dir']:
            raise MissingConfiguration(
                'Configuration setting "input_dir" is missing from the '
                'YAML configuration file and/or --input-dir is missing '
                'from the command line.')
        if not self.config['cubi_demux']['output_dir']:
            raise MissingConfiguration(
                'Configuration setting "output_dir" is missing from the '
                'YAML configuration and/or --output-dir is missing '
                'from the command line.')
        if not os.path.exists(self.config['cubi_demux']['input_dir']):
            raise InvalidConfiguration(
                'The input directory {} does not exist.'.format(
                    self.config['cubi_demux']['input_dir']))
        if (os.path.exists(self.config['cubi_demux']['output_dir']) and not
                self.config['cubi_demux']['continue']):
            raise InvalidConfiguration(
                'The output directory {} already exists!'.format(
                    self.config['cubi_demux']['output_dir']))

    def bcl2fastq_wrapper(self):
        """Return name of bcl2fastq wrapper to use."""
        if self.config['sample_sheet'][0]['rta_version'] == 2:
            return 'bcl2fastq2'
        else:
            return 'bcl2fastq'

    def out_prefix(self, fname):
        """Prefix ``fname`` with path to output file."""
        return os.path.join(self.config['cubi_demux']['output_dir'], fname)

    def _build_wrapper(self):
        """Build Bcl2FastqWrapper object."""
        if self.config['sample_sheet'][0]['rta_version'] == 2:
            return Bcl2FastqV2Wrapper(self.config)
        else:
            return Bcl2FastqV1Wrapper(self.config)

    def get_result_files_marker(self):
        """Returns list of result files for "rule all"."""
        return os.path.join(self.config['cubi_demux']['output_dir'],
                            'CUBI_Demux_Complete.txt')

    def get_result_files_demux(self):
        """Returns list of result files for "rule demux"."""
        return self.demux_wrapper.get_output_files()

    @listify
    def get_result_files_fastqc(self):
        """Returns list of result files for "rule fastqc"."""
        # We construct the FastQC output paths from the demux output paths.
        for path in self.demux_wrapper.get_output_files():
            ext = '.fastq.gz'
            if path.endswith(ext):
                folder = os.path.dirname(path)
                base = os.path.basename(path)[:-len(ext)]
                filename = base + '_fastqc.html'
                yield os.path.join(folder, 'qc', 'fastqc', filename)

    @listify
    def get_result_files_hts_screen(self):
        """Returns list of result files for "rule hts_screen"."""
        # We construct the hts_screen output paths from the demux output paths.
        for path in self.demux_wrapper.get_output_files():
            ext = '.fastq.gz'
            if path.endswith(ext):
                folder = os.path.dirname(path)
                base = os.path.basename(path)[:-len(ext)]
                filename = base + '_hts_screen.txt'
                yield os.path.join(folder, 'qc', 'hts_screen', filename)
                filename = base + '_kraken_report.gz'
                yield os.path.join(folder, 'qc', 'kraken', filename)

    def write_sample_sheet(self, output_path):
        """Write out sample sheet to ``output_path``."""
        with open(output_path, 'wt') as f:
            self.demux_wrapper.print_sheet(f)
