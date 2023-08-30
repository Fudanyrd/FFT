#include "complex_num.h"

using std::string;

bool not_minus(char c){
	return !isspace(c) && c!='-';
}

std::istream& getstr(std::istream& is,char* pt){
	char ch; char* i=pt;
	while(is.get(ch)){
		if(ch!='\n'){
			*i++ = ch;
		}
		else break;
	}
	*i = '\0';
	return is;
}

//output a complex num.
std::ostream& operator<<(std::ostream& os,const complex_num& x){
	std::streamsize prec=os.precision();
	os<<std::setprecision(6)<<x.real()<<std::setprecision(prec);
	os<<" + ";
	os<<std::setprecision(6)<<x.imagery()<<std::setprecision(prec);
	os<<"j";
	return os;
}

//READ in a COMPLEX NUMBER from INPUT STREAM.
std::istream& operator>>(std::istream& is, complex_num& x){
	char str[MAXLEN];
	getstr(is,str);
	char *end = std::find(str,str+MAXLEN,'j');
	if(end== str+40){
		x = complex_num(atof(str),0);
	}
	else{
		char *avail = std::find(str,end, '+');
		if(avail == end){
			avail= std::find_if(str,end,not_minus);
			avail= std::find(avail,end,'-');
			 if(avail==end){//pure imagery number
			 	x= complex_num(0,atof(str));
			 }
			 else{
			 	char _r_[MAXLEN/2], _i_[MAXLEN/2];
			 	std::copy(str,avail,_r_);
			 	std::copy(avail+1,end,_i_);
			 	x = complex_num(atof(_r_), -atof(_i_));
			 }
		}
		else{
			char _r_[MAXLEN/2], _i_[MAXLEN/2];
			std::copy(str,avail,_r_);
			std::copy(avail+1,end,_i_);
			x = complex_num(atof(_r_),atof(_i_));
		}
	}
	return is;
}

//OPERATOR COLLECTIONS
complex_num operator+(const complex_num& x,const complex_num& y){
	return complex_num(x.real()+y.real(),x.imagery()+y.imagery());
}
complex_num operator-(const complex_num& x,const complex_num& y){
	return complex_num(x.real()-y.real(),x.imagery()-y.imagery());
}
complex_num operator*(const complex_num& x,const complex_num& y){
	return complex_num(x.real()* y.real() - x.imagery()* y.imagery(),
	x.imagery()*y.real() + x.real()* y.imagery());
}
complex_num operator*(const complex_num& x,double d){
	return complex_num(x.real()*d,x.imagery()*d);
}
complex_num operator/(const complex_num& x,const complex_num& y){
	double l= y._length()*y._length();
	return complex_num((x.real()*y.real()+x.imagery()*y.imagery())/l,
	(x.imagery()*y.real()-x.real()*y.imagery())/l);
}
complex_num operator^(const complex_num& x,const complex_num& y){
	return (x * y)/(x + y);
}
complex_num operator~(const complex_num& x){// get conjugate
	return complex_num(x.real(),-x.imagery());
}

//convert a string into a complex number.
//NOTE: it must be of form a+bj or a-bj rather than bj+a or bj-a.
//NOTE: use 1j or 1.0j instead of j.
complex_num converter(const string& input){
	string::const_iterator end = std::find(input.begin(),input.end(),'j');
	if(end == input.end()) //input is not complex number.
		return complex_num(atof(input.c_str()),0.0);
	else{
		string::const_iterator avail = std::find(input.begin(),end,'+');
		if(avail == end){
			avail = std::find_if(input.begin(),end,not_minus);
			avail = std::find(avail,end,'-');
			if(avail == end)
				return complex_num(0.0,atof(input.c_str()));
			else{
				string temp1(input.begin(),avail), temp2(avail+1,end);
				return complex_num(atof(temp1.c_str()),-atof(temp2.c_str()));
			}
		}
		else{
			string temp1(input.begin(),avail), temp2(avail+1,end);
			return complex_num(atof(temp1.c_str()), atof(temp2.c_str()));
		}
	}
}
