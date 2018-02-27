#!/bin/bash

# Helper script for downloading Kraken DB, reference sequence and building BWA indices.

if [[ -z "$BCL2FASTQ_CHANNEL" ]]; then
    >&2 echo "Environment variable BCL2FASTQ_CHANNEL is empty"
    exit 1
fi

if [[ -z "$1" ]]; then
    >&2 echo "Usage: $0 output_dir"
    exit 1
fi

if [[ -e "$1" ]]; then
    >&2 echo "Output directory $1 already exists!"
    exit 1
fi

set -x

mkdir -p $1
cd $1

# Create Conda Environment -------------------------------------------------------------------------

cat <<EOF >environment.yml
name: build_references
channels:
- $BCL2FASTQ_CHANNEL
- bioconda
- conda-forge
- defaults
packages:
- wget
- bwa
EOF

conda create -n build_references

source activate build_references

set -euo pipefail

# Download Data ------------------------------------------------------------------------------------

mkdir -p hg19/tmp
cd hg19/tmp
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr{{1..22},X,Y,M}.fa.gz
cd ..
zcat tmp/*.gz >hg19.fa
bwa index hg19.fa
cd ..

mkdir NCBIM37
cd NCBIM37
wget ftp://ftp-mouse.sanger.ac.uk/ref/NCBIM37_um.fa
bwa index NCBIM37_um.fa
cd ..

mkdir danRer10
cd danRer10
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz
gzip -d danRer10.fa.gz
bwa index danRer10.fa
cd ..

mkdir dm6
cd dm6
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gzip -d dm6.fa.gz
bwa index dm6.fa
cd ..

mkdir -p sacCer3/tmp
cd sacCer3/tmp
test -e chromFa.tar.gz || \
    wget http://hgdownload-test.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
tar xf chromFa.tar.gz
cd ..
cat tmp/*.fa >sacCer3.fa
bwa index sacCer3.fa
cd ..

mkdir -p ecoli/tmp
cd ecoli/tmp
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/all_assembly_versions/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz
cd ..
gzip -cd tmp/GCA_000005845.2_ASM584v2_genomic.fna.gz >ecoli.fa
bwa index ecoli.fa
cd ..

mkdir -p phix/tmp
cd phix/tmp
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
tar xf PhiX_Illumina_RTA.tar.gz
cd ..
cp tmp/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa phix.fa
bwa index phix.fa
cd ..

mkdir -p UniVec/tmp
cd UniVec/tmp
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
cd ..
mv tmp/UniVec UniVec.fa
bwa index UniVec.fa
cd ..

mkdir -p kraken/tmp
cd kraken/tmp
wget http://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
cd ..
tar tf tmp/minikraken.tgz
cd ..

# Print Result -------------------------------------------------------------------------------------

cat <<EOF

!! Data Download & Index Build Complete !!

== Recommended config.yaml =========================================================================

# BEGIN: config.yaml
#
# Configuration for cubi-demux.
#
# This file is a YAML configuration file.  The default configuration is
# preconfigured for the BIH cluster and has to be adjusted accordingly.

# Configuration for the demultiplexing.
cubi_demux:
  input_dir: null  # path to input, override with "--input-dir"
  output_dir: null  # path to input, override with "--output-dir"
  cores: 8   # number of threads, override with "--cores"
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
  kraken_db: '$(readlink -f kraken/minikraken_20141208)'
  references:
  - name: 'H. sapiens'
    bwa_index: '$(readlink -f hs37/hs37.fa)'
  - name: 'M. musculus'
    bwa_index: '$(readlink -f NCBIM37/NCBIM37_um.fa)'
  - name: 'D. rerio'
    bwa_index: '$(readlink -f danRer10/danRer10.fa)'
  - name: 'D. melanogaster'
    bwa_index: '$(readlink -f dm6/dm6.fa)'
  - name: 'S. cerevisiae'
    bwa_index: '$(readlink -f sacCer3/sacCer3.fa)'
  - name: 'E. coli'
    bwa_index: '$(readlink -f ecoli/ecoli.fa)'
  - name: 'Phi X'
    bwa_index: '$(readlink -f phix/phix.fa)'
  - name: 'Univec'
    bwa_index: '$(readlink -f UniVec/UniVec.fa)'

# The sample sheet.  Either a path to the sample sheet or a dict with the
# sample sheet.  The path can can also be set with "--sample-sheet".
sample_sheet: null

# END: config.yaml

EOF
