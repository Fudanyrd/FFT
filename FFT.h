#pragma once
#ifndef _FFT_
#define _FFT_

#include <vector>
#include "complex_num.h"
#include <stdexcept>
#include <cmath>

complex_num Wn(unsigned int N,unsigned  int n,unsigned int j);
//calculate exp(-2*pi*sqrt(-1)/N*n*j)

std::vector<complex_num> gen_series(unsigned int N,unsigned  int j);
//generate a vector containing exp(-2*pi*n*j/N) with n(-[0,N-1]

complex_num DFT_j(unsigned int N,const std::vector<double>& Data,unsigned int j);
//calculate the j (- [0,N-1] element in the series.

std::vector<complex_num> DFT(unsigned int N,const std::vector<double>& Data);
//calculate the whole transformed series using slower algorithm.

complex_num FFT_j(unsigned int N,const std::vector<double>& Data,unsigned int j,
unsigned int Gap, unsigned int start);
//calculate the j(-[0,N-1] element in the series.

std::vector<complex_num> FFT(unsigned int N,const std::vector<double>& Data);
//calculate the whole series using FFT algorithm.

complex_num DFT_j(unsigned int N,const std::vector<double>& Data,unsigned int j);
std::vector<complex_num> DFT(unsigned int N,const std::vector<double>& Data);

complex_num FFT_j(unsigned int N, const std::vector<double>& Data, unsigned int j,
unsigned int Gap, unsigned int start);

void lengthenSeries(std::vector<double>& original);

complex_num** gen_array(unsigned int N, unsigned int m,const std::vector<double>& Data);

bool check(unsigned int N);

std::vector<complex_num> FFT(unsigned int N, const std::vector<double>& Data);

#endif//_FFT_
