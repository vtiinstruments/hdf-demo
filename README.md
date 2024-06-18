# hdf-demo
## Required packages
 - C/C++ compiler
 - CMake 3.14 or later

## Installation

### HDF Windows support
Install latest HDF[https://github.com/HDFGroup/hdf5] release. 
Refer to USING_HDF5_CMake.txt. In particular, set HDF5_PLUGIN_PATH and HDF5_ROOT

### HDF Linux support
Many distributions also contains HDF5 package. For example in ubuntu, HDF5 can be installed with
```bash
sudo apt-get install libhdf5-dev
```

## Building example
### Windows
 - Open CMake and configure the source  to this example
 - Create a build directory and configure  build to build binaries
### Linux
```bash
sudo apt install cmake libhdf5-serial-dev build-essential ninja-build tar curl zip unzip
cd path.to.this.example.directory
mkdir _build
pushd _build
cmake ..
make
popd
```
