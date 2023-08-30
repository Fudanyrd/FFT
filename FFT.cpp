#include "FFT.h"

complex_num Wn(unsigned int N,unsigned  int n,unsigned int j){
	double e = -2*_PI_*n*j/N;
	return complex_num(cos(e), sin(e));
}//calculate exp(-2*pi*sqrt(-1)/N*n*j)

std::vector<complex_num> gen_series(int N, int j){
	std::vector<complex_num> res;
	for(int n = 0; n != N; ++n){
		res.push_back(Wn(N,n,j));
	}
	return res;
}

complex_num DFT_j(unsigned int N,const std::vector<double>& Data,unsigned int j){
	complex_num res(0,0);
	for(int n=0; n!=N; ++n){
		res += Wn(N,n,j)*Data[n];
	}
	return res;
}
std::vector<complex_num> DFT(unsigned int N,const std::vector<double>& Data){
	std::vector<complex_num> res;
	for(int j=0;j!=N;++j){
		res.push_back(DFT_j(N,Data,j));
	}
	return res;
}

complex_num FFT_j(unsigned int N, const std::vector<double>& Data, unsigned int j,
unsigned int Gap, unsigned int start){
	if(N == 1){
		return complex_num(Data[start],0);
	}
	
	unsigned int m = N/2;
	unsigned int j0 = j/m;//in range [0,1]
	unsigned int j1 = j%m;// in range [0,m-1]
	if(j0%2){//j0 is odd
		return FFT_j(m, Data, j1, Gap*2,start) - Wn(N,1,j1)*FFT_j(m,Data,j1, Gap*2,start + Gap);
	}
	else{
		return FFT_j(m, Data, j1, Gap*2,start) + Wn(N,1,j1)*FFT_j(m,Data,j1, Gap*2,start + Gap);
	}
}

void lengthenSeries(std::vector<double>& original){
	size_t target_size = 1;
	while(target_size < original.size()) target_size *= 2;
	for(unsigned int i = 0;original.size() != target_size;i++){
		original.push_back(0.0);//to inicate the end of a series.
	}
	return;
}

complex_num** gen_array(unsigned int N, unsigned int m,const std::vector<double>& Data){
	//pt[start][j]
	unsigned int row = N/m, col = m;//Gap = row.
	complex_num** pt = new complex_num*[row];
	for(int i = 0;i!=row;++i){
		pt[i] = new complex_num[col];
	}
	
	if(m == 1){
		for(int r = 0;r != row; ++r){
			pt[r][0] = complex_num(Data[r],0.0);
		}
		return pt;
	}
	
	else{
		complex_num** result = gen_array(N,m/2,Data);
		if(result == 0) throw std::domain_error("invalid pointer");
		unsigned int half = col / 2;
		for(int r = 0;r!=row; ++r){
			for(int c = 0; c!=col;++c){
				if(c/half) pt[r][c] = result[r][c%half] - Wn(col,1,c%half)*result[r+row][c%half];
				else pt[r][c] = result[r][c%half] + Wn(col,1,c%half)*result[r+row][c%half];
			}
		}
		for(int r = 0;r!=row;++r){
			delete[] result[r];
		}
		delete[] result;
		
		return pt;
	}
}

bool check(unsigned int N){
	unsigned int i = 1;
	while(i < N) i*=2;
	return i == N;
}

std::vector<complex_num> FFT(unsigned int N, const std::vector<double>& Data){
	complex_num** pt;
	if(check(N))
		pt = gen_array(N,N,Data);
	else{
		std::vector<double> copy = Data;
		lengthenSeries(copy);
		unsigned int N_ = copy.size();
		pt = gen_array(N_,N_,copy);
		
		std::vector<complex_num> res(pt[0],pt[0] + N_);
		delete[] pt[0]; delete[] pt;
		return res;
	}
	
	std::vector<complex_num> res(pt[0],pt[0] + N);
	delete[] pt[0]; delete[] pt;
	return res;
}
