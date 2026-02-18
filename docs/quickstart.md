# Quickstart for qpgen

## Installing qpgen

Please follow the instruction below to install `qpgentools`

```sh
git clone --recursive https://github.com/hyunminkang/qpgen.git
cd qpgen
cd submodules
sh -x build.sh
cd ..
mkdir build
cd build
cmake ..
make
```

## List available tools

To list the available tools, run the following command:

```sh
../bin/qpgentools --help
```

If you encounter any difficulties, see [Install](install.md) for more details.