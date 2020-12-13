#!/bin/sh

# Create all the necessary directories and so on.
mkdir -p bin share include lib src
basedir=$(pwd)

cd src

# Build BLAS CBLAS LAPACK LAPACKE from git
git clone https://github.com/xianyi/openblas
cd openblas
make -j4
make install PREFIX=$basedir
cd ..

# Build arpack-ng from git
git clone https://github.com/opencollab/arpack-ng
cd arpack-ng
./bootstrap
export FCFLAGS=-O3
export FFLAGS=-O3
export CFLAGS=-O3
export CPPFLAGS=-O3
export CXXFLAGS=-O3
./configure --prefix=$basedir --with-blas=$basedir/lib/libopenblas.a --with-lapack=$basedir/lib/openblas.a --enable-icb
make -j4
make check
make install
cd ..
