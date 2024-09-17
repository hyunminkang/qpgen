# qpgen
Libraries and Command Line Tools for PLINK files

## Overview

`qpgen` is a collection of C++ libraries and tools to help analysis involving PLINK-formatted files. This repository is under development and is not fully documented yet.

## Installing qpgen

You can install `qpgen` by following the instructions below.

```bash
## clone the repository
git clone --recursive https://github.com/hyunminkang/qpgen.git
cd qpgen

## build the submodules
cd submodules
sh -x build.sh
cd ..

## build qpgen
mkdir build
cd build
cmake ..
make

## list available package
../bin/qpgentools --help
```

To see the usage of individual commands, type:

```bash
../bin/qpgentools [command] --help
```

To see the compiled library, type:

```bash
ls -l ../lib/libqpgen.a
```

## Advanced Options for Installation

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  -DHTS_LIBRARIES=/htslib_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

Finally, to build the binary, run

<pre>
$ make
</pre>
