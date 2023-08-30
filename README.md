#this is a realization of FFT algorithm in C++ programming language.

#construction:
FFT.h, IFFT.h: containing the solutions to FFT and IFFT algorithm.
complex_num.h: the definition of complex number.
series.h: generate a series of double value so that you can pass it to fft/ifft method.
main.cpp: a test for fft/ifft.

#NOTE that FFT method ONLY applied to series with length of pow(2,n).
e.g. 2,4,8,16,32,64,128,256,...,65536, 131072, ...

#how to compile it:
run the following command:
mingw32-make -f this.mk