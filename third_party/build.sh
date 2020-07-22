#!/bin/sh

mkdir -p bin share include lib src
basedir=$(pwd)

cd src

 Build fftw
tar xvzf $basedir/tars/fftw-3.3.8.tar.gz
cd fftw-3.3.8
./configure --prefix=$basedir
make
make install
cd ..

# Build BLAS CBLAS LAPACK LAPACKE
git clone https://github.com/xianyi/openblas
cd openblas
make -j4
make install PREFIX=$basedir
cd ..

# Build arpack
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
