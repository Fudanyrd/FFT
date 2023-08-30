#pragma once
#ifndef _IFFT_
#define _IFFT_

#include "FFT.h"

double Max_length(const std::vector<complex_num>& Data);//find the maximum length of a complex_num vector.
void Erase_little(std::vector<complex_num>& Data, double err = 1.0e-8);
//erase the elements that is below err(err is set 1.0e-8)
void Erase_little(std::vector<double>& Data, double err = 1.0e-8);

std::vector<double> IDFT(unsigned int N,const std::vector<complex_num>& Data);
//IDFT algorithm.
double IDFT_j(unsigned int N, const std::vector<complex_num>& Data, unsigned int j);
//calculate the j element of the transformed DFT series.
complex_num IFFT_j(unsigned int _N_,unsigned int N, const std::vector<complex_num>& Data, unsigned int j,
unsigned int Gap, unsigned int start);
//calculate the jth element of the transformed FFT series.
std::vector<double> IFFT(unsigned int N, const std::vector<complex_num>& Data);
complex_num Wk(unsigned int N,unsigned int j,unsigned int k);
//calculate exp(2*pi*i/N*j*k)

#endif//_IFFT_
