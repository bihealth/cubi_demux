==============
CUBI Demuxtool
==============

Automated demultiplexing + QC.
Companion to Flowcelltool.

--------
Overview
--------

``cubi-demuxtool`` is a single command line tool that will

- demultiplex Illumina data (both RTAv1 and RTAv2)
- perform QC with FastQC,
- screen reads for model organisms with HTS Screen (similar to FastqScreen),
- screen reads that do not align to model organisms with ``kraken`` for viral/bacterial DNA,
- aggregate the FastQC and HTS Screen reports using MultiQC.

------------
Installation
------------

The installation of ``cubi_demux`` itself is very simple but because of its nature, it has a dependency on the open source but not free ``bcl2fastq`` by Illumina.
We cannot distribute binary packages of that software so please bear with us through the following steps.

Building ``bcl2fastq`` Conda Packages
=====================================

First, setup Bioconda build installation (to ``~/miniconda3``, you might want to use a different path).

.. code-block:: shell

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    $ export PATH=$HOME/miniconda3/bin:$PATH

    $ conda config --add channels defaults
    $ conda config --add channels conda-forge
    $ conda config --add channels bioconda

Next, clone the ``cubi_demux`` Git repository

.. code-block:: shell

    $ git clone https://github.com/bihealth/cubi_demux.git
    $ cd cubi_demux
    $ git checkout v0.1.0

Next, we download the ``bcl2fastq`` source packages.

.. code-block:: shell

    $ mkdir -p downloads
    $ cd downloads
    $ wget \
        ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4.tar.bz2 \
        ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-tar.zip
    $ cd ..

Simply type (on the BIH cluster only at the moment):

.. code-block:: shell

    $ conda config --add channels http://cubi-conda.bihealth.org
    $ conda install cubi_demux

----------------
Overall Workflow
----------------

The overall workflow is as follows:

1. Generate a sample sheet YAML file using `CUBI Flowcelltool <https://github.com/bihealth/flowcelltool>`_.
2. Call ``cubi-demux`` pointing it to the sample sheet YAML, as well as the input and output folder.
3. Wait until the demultiplexing is complete.

And you're done.

-------------
Configuration
-------------

The default configuration is shown below together with documentation.

::

    # Configuration for cubi-demux.
    #
    # This file is a YAML configuration file.  The default configuration is
    # preconfigured for the BIH cluster and has to be adjusted accordingly.

    # Configuration for the demultiplexing.
    cubi_demux:
      input_dir: null  # path to input, override with `--input-dir`
      output_dir: null  # path to input, override with `--output-dir`
      cores: 8   # number of threads, override with `--cores`
      barcode_mismatches: null  # default is RTA version specific
      # Selecting lanes and tiles are mutually exclusive.
      lanes: null  # null or list of integers
      tiles: null  # tile specifications for bcl2fastq executable
      continue: false  # continue (do not break if output dir exists)

    # Configuration for the screening after demultiplexing.  You should provide
    # a list of BWA-indexed references and the path to a Kraken DB.  The data
    # will first be subsampled and screened versus the given model organisms'
    # genomes.  Unaligned reads will then be screened by Kraken.
    hts_screen:
      sample_rate: 0.001  # sample this rate of reads for screening
      kraken_db: '/fast/projects/cubit/current/static_data/app_support/kraken/0.10.5-cubi20160426/minikraken_20141208'
      references:
      - name: 'H. sapiens'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/GRCh37/hs37/hs37.fa'
      - name: 'M. musculus'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/NCBIM37/sanger/NCBIM37_um.fa'
      - name: 'D. rerio'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/danRer10/ucsc/danRer10.fa'
      - name: 'D. melanogaster'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/dm6/ucsc/dm6.fa'
      - name: 'S. cerevisiae'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/sacCer3/ucsc/sacCer3.fa'
      - name: 'E. coli'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/ecoli/GCA_000005845.2_ASM584v2/ecoli.fa'
      - name: 'Phi X 174'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/phix/illumina/phix.fa'
      - name: 'Univec 9'
        bwa_index: '/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.12/UniVec/9/UniVec.fa'

    # The sample sheet.  Either a path to the sample sheet or a dict with the
    # sample sheet.  The path can can also be set with `--sample-sheet`.
    sample_sheet: null

----------------------
Command Line Interface
----------------------

You can override certain settings from the configuration file directly on the command line.

::

    usage: cubi-demux [-h] [--version] [--verbose] [--work-in-output]
                    [--config CONFIG] [--sample-sheet SAMPLE_SHEET]
                    [--num-threads NUM_THREADS] [--input-dir INPUT_DIR]
                    [--output-dir OUTPUT_DIR]
                    [--barcode-mismatches BARCODE_MISMATCHES] [--cores CORES]
                    [--continue] [--lane LANES | --tiles TILES]

    optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit
    --verbose
    --work-in-output      Work output directory instead of temporary directory.
    --config CONFIG       Path to configuration YAML file. Default: /fast/users/
                            mholtgr/Development/demuxtool/cubi_demux/config.yaml
    --sample-sheet SAMPLE_SHEET
                            Path to sample sheet YAML file, overrides setting in
                            config YAML.
    --num-threads NUM_THREADS
                            Number of threads to run with, overrides setting in
                            config YAML.
    --input-dir INPUT_DIR
                            Path to input sequencer output folder, overrides
                            setting in config YAML.
    --output-dir OUTPUT_DIR
                            Path to output folder, overrides setting in config
                            YAML.
    --barcode-mismatches BARCODE_MISMATCHES
                            Mismatches to allow in barcode, default is 0 for v1
                            and 1 for v2
    --cores CORES         Number of cores to use, overrides setting in config
                            YAML.
    --continue            Do not exit if output dir exists but continue.

    Lane/Tile Selection:
    --lane LANES          Select individual lanes for demultiplexing; default is
                            to use all for which the sample sheet provides
                            information; provide multiple times for selecting
                            multiple lanes.
    --tiles TILES         Select tile regex; provide multiple times for multiple
                            regexes; conflicts with --lane
