g++ FFGrab.cpp -m64 -O3 -fPIC -shared -DMATLAB_MEX_FILE -DMX_COMPAT_32 -Iavbin/include -Iavbin/ffmpeg -I/usr/local/matlab/extern/include libavbin.so.64 -L/usr/local/matlab/bin/glnxa64 -lmx -lmex -o FFGrab.mexa64 -Wl,-R/usr/local/lib -Wl,-R.
#ld -shared -soname FFGrab.mexa64 -o FFGrab.o FFGrab.mexa64 -R /usr/local/lib -R . 

