---
layout: page
title: Installation
description: ~
---

`mtPGS` is implemented as a C++ package, which can be installed from GitHub by:

### Dependencies 
* C++ library: [Armadillo](https://arma.sourceforge.net/)
* Statistical genetics software: [PLINK](https://www.cog-genomics.org/plink/)


#### 1. Install `Armadillo` if necessary
Armadillo is a C++ library for linear algebra and scientific computing. Before installing Armadillo package, please ensure CMake tool, LAPACK and BLAS (or preferably OpenBLAS) are installed on your system. If these tools are not installed, you can download them from [http://www.cmake.org](http://www.cmake.org) (CMake tool) and [http://www.openblas.net/](http://www.openblas.net/) (OpenBLAS library). If you prefer to install the library and headers for Armadillo package in a user’s own directory, please use the option CMAKE_INSTALL_PREFIX as follows.

```
wget http://sourceforge.net/projects/arma/files/armadillo-9.200.7.tar.xz
tar xf armadillo-9.200.7.tar.xz
cd armadillo-9.200.7
cmake . -DCMAKE_INSTALL_PREFIX:PATH=/xxxx/armadillo-9.200.7
make && make install
```

#### 2. Install `mtPGS`
```
git clone https://github.com/xuchang0201/mtPGS.git
cd mtPGS/src
make
```
#### 3. Check the options included in the `mtPGS`
```
./mtPGS -h
```
