This is a repackaging of CADABRA code from CALTECH [http://www.vision.caltech.edu/cadabra]


# Build FFMPEG libraries

```
git submodule init
git submodule update
mkdir build
cd build
../libav/configure --prefix=.
make install
```

On windows to configure with mingw run:
```
sh ..\libav\configure --prefix=. --target-os=mingw32 --arch=x86_64
```

# Build CADABRA tools

## OSX
make -f Makefile.osx qtrak_preprocess qtrak_cluster analysis gVision
