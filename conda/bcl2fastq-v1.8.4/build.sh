#!/bin/bash

set -e

# Setup path variables
BINARY_HOME=$PREFIX/bin
PACKAGE_HOME=$PREFIX/opt/$PKG_NAME-$PKG_VERSION

# Create destinatio directories
mkdir -p $BINARY_HOME
mkdir -p $PACKAGE_HOME

# Export explicit paths to compilers
export CC=$PREFIX/bin/gcc
export CXX=$PREFIX/bin/g++

# Extract source tarball
tar -xf $RECIPE_DIR/../../downloads/bcl2fastq-1.8.4.tar.bz2

# Patch Boost (*sigh*)
mkdir -p tmp && \
    tar -C tmp -xf bcl2fastq/redist/boost_1_44_0.tar.gz && \
    pushd tmp/boost_1_44_0 && \
    patch -p1 <$RECIPE_DIR/boost.1.patch && \
    patch -p1 <$RECIPE_DIR/boost.2.patch && \
    popd && \
    tar -C tmp -czf bcl2fastq/redist/boost_1_44_0.tar.gz boost_1_44_0

# Patch bcl2fastq
pushd bcl2fastq && \
    patch -p1 <$RECIPE_DIR/bcl2fastq.patch && \
    popd

# Launch building
mkdir bcl2fastq-1.8.4-build
cd bcl2fastq-1.8.4-build
../bcl2fastq/src/configure --prefix=$PREFIX --with-cmake=$(which cmake)
make
make install
