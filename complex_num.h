#pragma once
#ifndef COMPLEX_NUM_H
#define COMPLEX_NUM_H

#include <cmath>
#define _PI_ acos(-1)
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

#ifndef _GRAMMAR_
bool not_minus(char c);
#endif

class complex_num{
private:
	double r;// real part;
	double i;//imagery part;
public:
	complex_num():r(0.0),i(0.0){
	}		
	complex_num(double dat):r(dat),i(0.0){
	}
	complex_num(double real,double imagery):
	r(real),i(imagery){}
	
	double& real(){
		return r;
	}
	const double& real()const{ return r;}
	double& imagery(){
		return i;
	}
	const double& imagery()const{
		return i;
	}
	double _length()const{
		return sqrt(r*r + i*i);
	}
	double _angle_rad()const{
		if(fabs(r)<1.0e-6){
			if(i>0.0)return acos(-1)/2;
			else return -acos(-1)/2;
		}
		if(r>= 0.0){
			return atan(i/r);
		}
		else{
			if(i>=0.0){
				return atan(i/r)+ _PI_;
			}
			else return atan(i/r)- _PI_;
		}
	}
	double _angle_deg()const{
		return _angle_rad()/_PI_*180;
	}
	complex_num& operator=(const complex_num& x){
		r = x.r;
		i = x.i;
		return *this;
	}
	complex_num& operator+=(const complex_num& x){
		r += x.r;
		i += x.i;
		return *this;
	}
	complex_num& operator*=(const complex_num& x){
		r = r*(x.r) - i*(x.i);
		i = r*(x.i) + i*(x.r);
		return *this;
	}
	complex_num& operator*=(double d){
		r*=d; i*=d;
		return *this;
	}
	
	std::ostream& report(std::ostream& os){
		os<<_length()<<" {"<<_angle_deg();
		return os;
	}
	
	operator double()const{ return this->r; }
};

//output a complex num.
std::ostream& operator<<(std::ostream& os,const complex_num& x);

#ifndef _GETSTR_
#define _GETSTR_
#define MAXLEN 60
std::istream& getstr(std::istream& is,char* pt);
#endif//_GETSTR_

//read in a complex number.
std::istream& operator>>(std::istream& is, complex_num& x);

//operator collections.
complex_num operator+(const complex_num& x, const complex_num& y);
complex_num operator-(const complex_num& x, const complex_num& y);
complex_num operator*(const complex_num& x, const complex_num& y);
complex_num operator*(const complex_num& x, double d);
complex_num operator/(const complex_num& x, const complex_num& y);
complex_num operator^(const complex_num& x, const complex_num& y);
complex_num operator~(const complex_num& x);

//convert a string into a complex number.
complex_num converter(const std::string& input);

#endif//COMPLEX_NUM_H
