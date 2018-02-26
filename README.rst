==============
CUBI Demuxtool
==============

Automated demultiplexing + quality control.
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

----------------
Overall Workflow
----------------

After installation (see below), the overall workflow is as follows:

1. Generate a global configuration file (you can reuse this file for future uses).
2. Generate a sample sheet YAML file using `CUBI Flowcelltool <https://github.com/bihealth/flowcelltool>`_.
3. Call ``cubi-demux`` pointing it to the sample sheet YAML, as well as the input and output folder.
4. Wait until the demultiplexing is complete.

And you're done.

------------
Installation
------------

The installation of ``cubi_demux`` itself is very simple but because of its nature, it has a dependency on the open source but not free ``bcl2fastq`` by Illumina.
We cannot distribute binary packages of that software so please bear with us through the following steps.

Prerequisites
=============

- Install Docker (e.g., `following the instructions from Docker.com <https://docs.docker.com/install/>`_.
- Get the ``bioconda-utils-build-env`` container:

    .. code-block:: shell

    $ docker pull bioconda/bioconda-utils-build-env

Building ``bcl2fastq`` Conda Packages
=====================================

First, setup Bioconda build installation (to ``~/miniconda3``, you might want to use a different path).

You can do this on a different server from the one that you will execute ``cubi_demux`` on.

.. code-block:: shell

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    $ export PATH=$HOME/miniconda3/bin:$PATH

    $ conda config --add channels defaults
    $ conda config --add channels conda-forge
    $ conda config --add channels bioconda

    $ conda install conda-build

Next, clone the ``cubi_demux`` Git repository

.. code-block:: shell

    $ git clone https://github.com/bihealth/cubi_demux.git
    $ cd cubi_demux
    $ git checkout master

Next, we download the ``bcl2fastq`` source packages.

.. code-block:: shell

    $ mkdir -p downloads
    $ cd downloads
    $ wget \
        ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4.tar.bz2 \
        ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-tar.zip
    $ cd ..

Then, we build the conda packages but inside the ``bioconda-utils-build-env`` container:

.. code-block:: shell

    host $ mkdir packages
    host $ docker run -v $PWD:/cubi_demux -i -t bioconda/bioconda-utils-build-env /bin/bash

    container $ cd /cubi_demux

    container $ conda build conda/bcl2fastq-v1.8.4
    [...]
    container $ cp /opt/conda/conda-bld/linux-64/bcl2fastq-1.8.4-pl5.20.3_4.tar.bz2 packages

    container $ conda build conda/bcl2fastq2-v2.17.1.14
    [...]
    container $ cp /opt/conda/conda-bld/linux-64/bcl2fastq2-2.17.1.14-2.tar.bz2 packages

We now have to create a local conda repository containing these packages somewhere on the file system **where you want to run demultiplexing**.
For example, this would be on the demultiplexing server in the case of working on one server or the shared cluster file system in the case of working with HPC.
For the sake of simplicity, we assume this is the same as the build machine and create the repository in your home folder:

.. code-block:: shell

    $ mkdir -p $HOME/local_channel/linux-64
    $ cp packages/* $HOME/local_channel/linux-64
    $ conda index $HOME/local_channel/linux-64

Building ``cubi_demux`` Package
===============================

If we were able to redistribute Illumina ``bcl2fastq`` packages via Bioconda, this would be much simpler.
However, here is how to build a conda package ``cubi_conda``.

Note that this will have the path to your ``local_channel`` baked into the package.
This means that it cannot be easily ported to another machine that does not have the ``local_channel`` Conda channel as well **in the same location**.

You can build the package on a different one from that you use for demultiplexing, but you have to specify the path to ``local_channel`` on the machine that you will perform demultiplexing one (the deployment machine).
As ``cubi_demux`` relies on Snakemake's conda integration, the instructions above are complicated but probably the best ones available.

.. code-block::

    $ export BCL2FASTQ_CHANNEL=file://$HOME/local_channel
    $ conda build conda/cubi_demux
    [...]
    anaconda upload /bioconda/2018-02/miniconda3/conda-bld/linux-64/cubi_demux-0.1.1-py36_1.tar.bz2
    [...]
    $ cp \
        /bioconda/2018-02/miniconda3/conda-bld/linux-64/cubi_demux-0.1.1-py36_1.tar.bz2 \
        $HOME/local_channel/linux-64
    $ conda index $HOME/local_channel/linux-64

Installing ``cubi_demux``
=========================

First, make your ``local_channel`` Conda channel known to conda

.. code-block::

    $ conda config --add channels file://$HOME/local_channel

Then, you can install ``cubi_demux``:

.. code-block::

    $ conda install cubi_demux

Create Data for Read Screening and Kraken
=========================================

-------------
Configuration
-------------

The default configuration is shown below together with documentation.

.. code-block:: yaml

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
