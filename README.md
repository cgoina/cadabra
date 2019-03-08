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

make -f Makefile.<platform_extension>

