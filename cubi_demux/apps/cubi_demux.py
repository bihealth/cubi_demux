# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off demultiplexing

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import os
import sys
import tempfile

import snakemake

from ..__main__ import MissingConfiguration
from .. import __version__

#: Path to default configuration, from package (the reason why the package is
#: not zip_safe.
DEFAULT_CONFIG_YAML = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'config.yaml'))

#: Path to the Snakefile.
PATH_SNAKEFILE = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'Snakefile'))


def run(args):
    """Main program entry point after parsing the command line"""
    with tempfile.TemporaryDirectory(prefix='cubi-demux') as tmpdir:
        print('Temporary directory: {}'.format(tmpdir), file=sys.stderr)
        args = dict((k, v) for k, v in vars(args).items() if v is not None)
        snakemake.snakemake(
            PATH_SNAKEFILE, configfile=args['config'], config=args,
            verbose=args['verbose'])


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    # Setup the parser.
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--version', action='version', version='%%(prog)s %s' % __version__)
    parser.add_argument(
        '--verbose', action='store_true', default=False)
    parser.add_argument(
        '--config', type=str, default=DEFAULT_CONFIG_YAML,
        help='Path to configuration YAML file.')
    parser.add_argument(
        '--sample-sheet', type=str, default=None,
        help='Path to sample sheet YAML file, if not in configuration.')
    parser.add_argument(
        '--num-threads', type=int, default=None,
        help='Number of threads to run with.')
    parser.add_argument(
        '--input-dir', type=str, default=None,
        help='Path to input sequencer output folder.')

    # Parse command line arguments.
    args = parser.parse_args(argv)

    # Start the processing.
    run(args)


if __name__ == '__main__':
    sys.exit(main())
