This is a repackaging of CADABRA code from CALTECH [http://www.vision.caltech.edu/cadabra]


# Build FFMPEG libraries

## Prepare the build (all platforms)
```
git submodule init
git submodule update
mkdir build
cd build
# configure libav - see the platform specific instructions
make install
```

## Configure LIBAV build instructions

### MacOSX
```
../libav/configure --prefix=.
```

### Win64 with mingw
```
sh ..\libav\configure --prefix=. --target-os=mingw32 --arch=x86_64 --enable-pic
```

### Linux
```
../libav/configure --prefix=. --enable-pic
```

# Build CADABRA tools

On linux make sure the gcc-6 binary and libraries are in the PATH and LD_LIBRARY_PATH respectively.
gcc-6 is the current version supported by mex. 

```
source env.linux.sh
```
```
make -f Makefile.<platform_extension>
```


# Running CADABRA tools.

In order to run qtrak you will need MATLAB-2018b runtime which can be installed following the directions provided [here|https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html]