main: main.o IFFT.o FFT.o complex_num.o
	g++ -o main main.o IFFT.o FFT.o complex_num.o

complex_num.o: complex_num.h
	g++ -c complex_num.cpp
FFT.o: FFT.h
	g++ -c FFT.cpp
IFFT.o: IFFT.h
	g++ -c IFFT.cpp
main.o: main.cpp
	g++ -c main.cpp	

clear:
	rm *.o