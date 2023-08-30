#include <iostream>
#include <ios>
#include <iomanip>
#include <string>
#include <vector>
#include <iomanip>
#include <list>
#include <cstdlib>
#include <ctime>

using std::setprecision;	using std::setw;
using std::setiosflags;		using std::cin;
using std::string;			using std::vector;
using std::cout;			using std::endl;
using std::setw;			using std::list;

#include <fstream>
using std::ifstream;		using std::ofstream;
using std::ostream;

#include "FFT.h"
#include "IFFT.h"
#include "Series.hpp"

#define N 131072u

template <typename _Tp>
void write_vector(std::ostream& os, const vector<_Tp>& v1,const vector<complex_num>& v, const vector<_Tp>& v2){
	os << "original,, real part, imagery part, ,transformed\n";
	for(int i = 0;i!=v1.size(); ++i){
		os << v1[i] << ", , "<< v[i].real() << ", " << v[i].imagery() << ", , " << v2[i] << endl;
	}
	return;
}

int main(int argc,char** argv){
	srand((unsigned)time(0));
	vector<double> res = seriesCreater(N,0.0,4.0);
	time_t begin, end1, end2;
	system("pause");
	
	ofstream out("hello.csv");
	
	begin = time(0);
	vector<complex_num> transformed = FFT(res.size(),res); 
	end1 = time(0);
	
	Erase_little(transformed);
	vector<double> fft = IFFT(N,transformed);
	end2 = time(0);
	write_vector(out,res,transformed,fft);
	Erase_little(fft);
	
	cout << "time 1: " << difftime(end1,begin) << endl;
	cout << "time 2: " << difftime(end2,end1) << endl;
	
	//system("pause");
	return 0;
} 
