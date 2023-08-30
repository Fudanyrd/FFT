#include "FFT.h"
#include "IFFT.h"

double Max_length(const std::vector<complex_num>& Data){//find the maximum length of a complex_num vector.
	double res = 0.0;
	std::vector<complex_num>::size_type dat_sz = Data.size();
	for(std::vector<complex_num>::size_type i = 0;i!=dat_sz;++i){
		if(res < Data[i]._length())
			res = Data[i]._length();
	}
	return res;
}
void Erase_little(std::vector<complex_num>& Data, double err){
	std::vector<complex_num>::size_type dat_sz = Data.size();
	for(std::vector<complex_num>::size_type i = 0;i!=dat_sz;++i)
		if(Data[i]._length() < err)
			Data[i] = complex_num(0.0,0.0);
	return;
}//erase the elements that is below err(err is set 1.0e-8)
void Erase_little(std::vector<double>& Data, double err){
	std::vector<double>::size_type dat_sz = Data.size();
	for(std::vector<double>::size_type i = 0;i!=dat_sz;++i){
		if(fabs(Data[i]) < err)
			Data[i] = 0.0;
	}
	return;
}

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

//Now I will realize the funcions.
complex_num Wk(unsigned int N,unsigned int j,unsigned int k){
	double base = 2*_PI_/N;
	return complex_num(cos(base*j*k), sin(base*j*k));
}

double IDFT_j(unsigned int N, const std::vector<complex_num>& Data, unsigned int k){
	complex_num temp, res(0.0,0.0);
	for(unsigned int j=0;j!=N;++j){
		temp = Data[j]*Wk(N,j,k);
		res.real()+=temp.real()/N;
		res.imagery()+=temp.imagery()/N;
	}
	
	return res.real();
}

std::vector<double> IDFT(unsigned int N, const std::vector<complex_num>& Data){
	std::vector<double> res;
	for(unsigned int j =0u; j!=N;++j){
		res.push_back(IDFT_j(N,Data,j));
	}
	return res;
}

std::vector<complex_num> odd_series(const std::vector<complex_num>& Data){
	std::vector<complex_num> res;
	std::vector<complex_num>::size_type dat_sz = Data.size(); 
	for(std::vector<complex_num>::size_type i = 1; i < dat_sz; i += 2){
		res.push_back(Data[i]);
	}
	return res;
}
std::vector<complex_num> even_series(const std::vector<complex_num>& Data){
	std::vector<complex_num> res;
	std::vector<complex_num>::size_type dat_sz = Data.size(); 
	for(std::vector<complex_num>::size_type i = 0; i < dat_sz; i += 2){
		res.push_back(Data[i]);
	}
	return res;
}
complex_num IFFT_j(unsigned int _N_,unsigned int N, const std::vector<complex_num>& Data, unsigned int j,
unsigned int Gap, unsigned int start){
	if(N==1) return complex_num(Data[start].real()/_N_,Data[start].imagery()/_N_);
	unsigned int m = N/2;
	unsigned int j0 = j/m;//in range [0,1]
	unsigned int j1 = j%m;//in range [0,m-1]
	if(j0%2){//j0 is odd
		return IFFT_j(_N_,m, Data, j1,Gap*2,start) - Wk(N,1,j1)*IFFT_j(_N_,m,Data,j1,Gap*2,start+Gap);
	}
	else{
		return IFFT_j(_N_,m, Data, j1,Gap*2,start) + Wk(N,1,j1)*IFFT_j(_N_,m,Data,j1,Gap*2,start+Gap);
	}
}// Do NOT forget to divide N!

complex_num** gen_IFFT(unsigned int N,unsigned int m,const std::vector<complex_num>& Data){
	unsigned int row = N/m, col = m;//Gap = row.
	complex_num** pt = new complex_num*[row];
	for(int i = 0;i!=row;++i){
		pt[i] = new complex_num[col];
	}
	
	if(m == 1){
		for(int r = 0;r != row; ++r){
			pt[r][0] = complex_num(Data[r].real()/N,Data[r].imagery()/N);
		}
		return pt;
	}
	
	else{
		complex_num** result = gen_IFFT(N,m/2,Data);
		if(result == 0) throw std::domain_error("invalid pointer");
		unsigned int half = col / 2;
		for(int r = 0;r!=row; ++r){
			for(int c = 0; c!=col;++c){
				if(c/half) pt[r][c] = result[r][c%half] - Wk(col,1,c%half)*result[r+row][c%half];
				else pt[r][c] = result[r][c%half] + Wk(col,1,c%half)*result[r+row][c%half];
			}
		}
		for(int r = 0;r!=row;++r){
			delete[] result[r];
		}
		delete[] result;
		return pt;
	}
}

void lengthenSeries(std::vector<complex_num>& original){
	size_t target_size = 1;
	while(target_size < original.size()) target_size *= 2;
	for(unsigned int i = 0;original.size() != target_size;i++){
		original.push_back(0.0);
	}
	return;
}

std::vector<double> IFFT(unsigned int N, const std::vector<complex_num>& Data){
	std::vector<double> res;
	//lengthenSeries(Data);
	//N = Data.size();
	complex_num** pt;
	if(check(N)){
		pt = gen_IFFT(N,N,Data);
	}
	else{
		std::vector<complex_num> copy = Data;
		lengthenSeries(copy);
		unsigned int N_ = copy.size();
		pt = gen_IFFT(N_,N_, copy);
		
		for(unsigned int i = 0; i!=N_;++i){
			res.push_back(pt[0][i].real());
		}
		delete[] pt[0]; delete[] pt;
		return res;
	}
	
	for(unsigned int i = 0;i!=N;++i){
		res.push_back(pt[0][i].real());
	}
	delete[] pt[0]; delete[] pt;
	return res;
}