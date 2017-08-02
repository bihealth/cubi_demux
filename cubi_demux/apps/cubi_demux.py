# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off demultiplexing

As the entry point Snakefile would always be the same anyway, it is much
more convenient to wrap the call to snakemake itself.
"""

import argparse

import snakemake

from .. import __version__


def main(argv=None):
    """Main program entry point, starts parsing command line arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument('--version', action='version', version='%%(prog)s %s' % __version__)

    group = parser.add_argument_group(
        'Snakemake Basics',
        'Basic arguments from Snakemake')
    group.add_argument('targets', metavar='FILE', nargs='*', help='The files to create')
    group.add_argument(
        '-S', '--summary', action='store_true', default=False,
        help='Print paths that are to be created')
    group.add_argument(
        '-d', '--directory', default=os.getcwd(),
        help='Path to directory to run in, default is cwd')
    group.add_argument(
        '--cores', type=int,
        help='Number of cores to use for local processing')
    group.add_argument(
        '--unlock', action='store_true', default=False,
        help='Unlock working directory')
    group.add_argument(
        '--rerun-incomplete', action='store_true', default=False,
        help='Rerun incomplete jobs')

    group = parser.add_argument_group(
        'Snakemake Verbosity / Debugging',
        'Arguments from Snakemake that are useful for debugging, such as '
        'increasing verbosity')
    group.add_argument(
        '-p', '--printshellcmds', action='store_true', default=False,
        help='Print shell commands')
    group.add_argument(
        '--verbose', action='store_true', default=False,
        help='Enable verbose mode')
    group.add_argument(
        '--debug', action='store_true', default=False,
        help='Enable debugging')
    group.add_argument(
        '-n', '--dryrun', action='store_true', default=False,
        help='Simulate/run workflow without executing anything')
    group.add_argument(
        '--print-compilation', action='store_true', default=False,
        help='Print compilation and stop')

    group = parser.add_argument_group(
        'DRMAA/Cluster Configuration',
        'Arguments for enabling and controllingn basic DRMAA execution')
    group.add_argument(
        '--cubi-pipeline-use-drmaa', action='store_true', default=False,
        help='Enables running the pipeline with DRMA via Snakemake')
    group.add_argument(
        '--cubi-pipeline-drmaa-jobs', type=int, default=100,
        help='Number of DRMAA jobs to run')

    group = parser.add_argument_group(
        'DRMAA Robustness',
        'Snakemake settings to increase robustness of executing the pipeline '
        'during execution in DRMAA mode')
    group.add_argument(
        '--restart-times', type=int, default=5,
        help='Number of times to restart jobs automatically')
    group.add_argument(
        '--max-jobs-per-second', type=int, default=10,
        help='Maximal number of jobs to launch per second in DRMAA mode')

    group = parser.add_argument_group(
        'CUBI Pipeline Miscellaneous',
        'Various options specific to the CUBI pipeline')
    group.add_argument(
        '--cubi-pipeline-self-test', action='store_true', default=False,
        help='Perform self-test and exit')
    group.add_argument(
        '--step', type=str, metavar='STEP',
        choices=sorted(STEP_TO_MODULE.keys()),
        help='The type of the step to run')
    group.add_argument(
        '--no-use-conda', dest='use_conda', action='store_false', default=True,
        help='Disable usage of conda')

    args = parser.parse_args(argv)

    if args.cubi_pipeline_self_test:
        run_self_test()
    else:
        if not args.step:
            for cfg in CONFIG_FILES:
                path = os.path.join(args.directory, cfg)
                if not os.path.exists(path):
                    continue
                with open(path, 'rt') as f:
                    data = yaml.round_trip_load(f.read())
                try:
                    args.step = data['pipeline_step']['name']
                    break
                except KeyError:
                    print(('INFO: could not pick up pipelinestep/name '
                           'from {}').format(path), file=sys.stderr)
        if not args.step:
            parser.error('the following arguments are required: --step')
        try:
            return run(args)
        except DRMAANotAvailable:
            parser.error(
                'DRMAA not available but given --cubi-pipeline-use-drmaa')


if __name__ == '__main__':
    sys.exit(main())
