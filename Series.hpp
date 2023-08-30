#ifndef _SERIES_
#define _SERIES_

#include <vector>
#include <cmath>
#define pi acos(-1)

double func(double x){
	return (x-(int)x)+1.0;
}
//2^15 = 32768
//2^16 = 65536
//2^17 = 131072
std::vector<double> seriesCreater(unsigned int N,double begin = 0.0,double end = 1.0){
	std::vector<double> result;
	unsigned int count = 0u;
	double gap = end - begin;
	for(count = 0;count != N; ++count){
		result.push_back(func(begin + gap/N*count));
	}
	return result;
}

#endif//_SERIES_

