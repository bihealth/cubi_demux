# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off demultiplexing

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse
import json
import os
import sys
import tempfile

import snakemake

from ruamel.yaml import YAML

from .. import __version__

#: File name of the configuration.
CONFIG_FILE = 'CUBI_Demux_Config.json'

#: Path to default configuration, from package (the reason why the package is
#: not zip_safe.
DEFAULT_CONFIG_YAML = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'config.yaml'))

#: Path to the Snakefile.
PATH_SNAKEFILE = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'Snakefile'))


def write_config_file(args):
    """Write configuration file into output directory.

    Returns path to config file and config.
    """
    # Create output directory.
    os.makedirs(args.output_dir, exist_ok=True)
    # Load configuration.
    with open(args.config, 'rt') as f:
        config = YAML(typ='safe').load(f)
    # Override configuration with command line.
    args_dict = dict((k, v) for k, v in vars(args).items() if v is not None)
    keys = ('input_dir', 'output_dir', 'lanes', 'tiles', 'barcode_mismatches',
            'cores', 'continue')
    for key in keys:
        if args_dict.get(key):
            config['cubi_demux'][key] = args_dict[key]
    if args.sample_sheet:
        config['sample_sheet'] = args.sample_sheet
    # Load the sheet information and put into configuration.
    if isinstance(config['sample_sheet'], str):
        with open(config['sample_sheet'], 'rt') as f:
            config['sample_sheet'] = YAML(typ='safe', pure=True).load(f)
    # Write out configuration as JSON.
    out_path = os.path.join(args.output_dir, CONFIG_FILE)
    print('Writing configuration+sheet to {}'.format(out_path),
          file=sys.stderr)
    with open(out_path, 'wt') as f:
        json.dump(config, f, sort_keys=True, indent=4)
    return out_path, config


def work(args, workdir, config_path, config):
    """Perform the actual work."""
    # Build arguments to pass to ``snakemake``.
    argv = [
        '--snakefile', PATH_SNAKEFILE,
        '--directory', workdir,
        '--configfile', config_path,
        '--cores', config['cubi_demux']['cores'],
        '--use-conda'
    ]
    if args.verbose:
        argv.append('--verbose')
    for k, v in vars(args).items():
        argv += ['--config', '{}={}'.format(k, v)]
    argv = list(map(str, argv))
    print('Executing "snakemake {}"...'.format(' '.join(argv)),
          file=sys.stderr)
    snakemake.main(argv)


def run(args):
    """Main program entry point after parsing the command line"""
    # Read in configuration file, override with command line settings, write
    # to configuration file.
    config_path, config = write_config_file(args)
    # Run Snakemake in temporary directory.
    if args.work_in_output:
        print('Working in output directory: {}'.format(args.output_dir),
              file=sys.stderr)
        work(args, args.output_dir, config_path, config)
    else:
        with tempfile.TemporaryDirectory(prefix='cubi-demux') as tmpdir:
            print('Working in tmp directory: {}'.format(tmpdir),
                  file=sys.stderr)
            work(args, tmpdir, config_path, config)


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    # Setup the parser.
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--version', action='version', version='%%(prog)s %s' % __version__)
    parser.add_argument(
        '--verbose', action='store_true', default=False)

    group = parser.add_argument_group('Input / Output Related')
    group.add_argument(
        '--input-dir', type=str, default=None,
        help=('Path to input sequencer output folder, overrides setting in '
              'config YAML.'))
    group.add_argument(
        '--output-dir', type=str, default=None,
        help='Path to output folder, overrides setting in config YAML.')
    group.add_argument(
        '--sample-sheet', type=str, default=None,
        help=('Path to sample sheet YAML file, overrides setting in '
              'config YAML.'))
    group.add_argument(
        '--config', type=str, default=DEFAULT_CONFIG_YAML,
        help='Path to configuration YAML file. Default: {}'.format(
            DEFAULT_CONFIG_YAML))

    group = parser.add_argument_group('Misc. Configuration')
    group.add_argument(
        '--barcode-mismatches', type=int, default=None,
        help=('Mismatches to allow in barcode, default is 0 for v1 and '
              '1 for v2'))
    group.add_argument(
        '--cores', type=int, default=None,
        help='Number of cores to use, overrides setting in config YAML.')

    # Configuration related to lane and tile selection.
    group = parser.add_argument_group('Lane/Tile Selection')
    mutex = group.add_mutually_exclusive_group()
    mutex.add_argument(
        '--lane', type=int, default=[], action='append', dest='lanes',
        help=('Select individual lanes for demultiplexing; default is to '
              'use all for which the sample sheet provides information; '
              'provide multiple times for selecting multiple lanes.'))
    mutex.add_argument(
        '--tiles', type=str, default=[], action='append', dest='tiles',
        help=('Select tile regex; provide multiple times for multiple '
              'regexes; conflicts with --lane'))

    group = parser.add_argument_group('Snakemake Execution Related')
    group.add_argument(
        '--continue', action='store_true', default=None,
        help='Do not exit if output dir exists but continue.')
    group.add_argument(
        '--work-in-output', action='store_true', default=False,
        help='Work output directory instead of temporary directory.')

    # Parse command line arguments.
    args = parser.parse_args(argv)

    # Start the processing.
    run(args)


if __name__ == '__main__':
    sys.exit(main())
