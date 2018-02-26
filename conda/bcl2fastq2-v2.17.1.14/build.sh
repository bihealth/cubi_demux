#!/bin/bash

set -e

# Setup path variables
BINARY_HOME=$PREFIX/bin
PACKAGE_HOME=$PREFIX/opt/$PKG_NAME-$PKG_VERSION

# Create destinatio directories
mkdir -p $BINARY_HOME
mkdir -p $PACKAGE_HOME

# Extract source tarball
tar --strip-components 1 \
    -xf $RECIPE_DIR/../../downloads/$PKG_NAME-v$PKG_VERSION.tar.gz

# Patch boost (*sigh*)
BDIR=$(mktemp -d)
trap "rm -rf $BDIR" EXIT
tar -C $BDIR -xf redist/boost_1_54_0.tar.bz2
pushd $BDIR/boost_1_54_0 &&
    patch -p2 <$RECIPE_DIR/boost.patch &&
    popd &&
    tar -C $BDIR -cjf redist/boost_1_54_0.tar.bz2 boost_1_54_0

# Launch building
export CPLUS_INCLUDE_PATH=$PREFIX/include:/usr/include/x86_64-linux-gnu
export LIBRARY_PATH=$PREFIX/lib:/usr/lib/x86_64-linux-gnu
cd src
./configure --with-cmake=$(which cmake) --force-builddir --prefix=$PREFIX
make
make install
