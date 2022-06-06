#pragma once
#include <algorithm>
#include <vector>
#include <string>
#include <regex>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include "constants.hpp"

// utils for arg parse
namespace argutil {	// got from stackoverflow
	char* getCmdOption(char ** begin, char ** end, const std::string & option);
	bool cmdOptionExists(char** begin, char** end, const std::string& option);
}


// utils for io processes
namespace ioutil {
	void split(std::vector<std::string> &VS, const std::string &S);	// split S by <space> and '\t', then store splitted strings in VS
	void split(std::vector<std::string> &VS, const std::string &S, const char sp);	// split S by char sp, then store splitted strings in VS
	void eq2space(std::string &S);	// replace '=' in S with <space>
	void qt2space(std::string &S);	// replace quatation marks to <space>
	void replacea2b(std::string &S, const std::string &a, const std::string &b);	// replace pattrens of a in S to b
	void trim_comment(std::string &S);	// trim comment contents of S (following '#' or '!')
	void trim_space(std::string &S); // trim <space> or '\t' at the head and end of str
	void cap2low(std::string &S);	// transform capital letter in S to lowercase
	bool stobool(const std::string &S);	// return bool value of S
	dcomplex stocx(const std::string &S, int n);	// return complex number from string S, S should be in the form "e^m" (meaning exp(im/n*2pi)), or just a real number
	bool read_index_from_str(const std::string &S, int &ibegin, int &iend);	// S should be in the form of single interger, or "ibegin-iend", return false for a bad input
	double simpexptod(const std::string &Exp);	// calculate the double value for specific expression
}


// utils for some math functions
namespace mathutil {
	void mod_1(double &val, const double prec = infinit_small);	// set val mod 1 (fall into [0,1))
	void mod_1(arma::dcolvec3 &V);
	int mod_1_re_int(double &val, const double prec = infinit_small);	// set val mod 1 (fall into [0,1)) and return integer part of val
	int get_rand_int(int min, int max);	// return random interger in [min, max]
	double get_rand_double(double min, double max);
	arma::cx_colvec proj(const arma::cx_colvec u, const arma::cx_colvec v);	// return projection of v on u
	int mod(const int val, const int interval_end, const int interval_begin = 0);
}

// utils for arma matrix
namespace armautil {
	bool is_eye(const arma::mat &M);
	bool is_eye(const arma::imat33 &M);
	void std33_to_dmat33(arma::dmat33 &M, const double M0[3][3]);
	void dmat33_to_std33(double M[3][3], const arma::dmat33 &M0);
	int find_dcolvec3_from_list(const arma::dcolvec3 &V, const std::vector<arma::dcolvec3> &Vlist, const double prec = infinit_small);
}


void print_version();

// templates for some frequently used math functions
template <typename T> int sgn(const T &val, const double threshold = infinit_small) {
	return (threshold < val) - (val < -threshold);
}

template <typename T> int notzero(const T &val, const double threshold = infinit_small) {
	return (threshold < val) + (val < -threshold);
}

int notzero(const dcomplex &val, const double threshold = infinit_small);
void roundzero(dcomplex &var, const double threshold = infinit_small);

#ifdef __DEBUG
// employ some debuggers
void debug_ccxu();
void debug_gxzhi();
#define DEBUG_CCXU debug_ccxu();
#define DEBUG_GXZHI debug_gxzhi();
#else
#define DEBUG_CCXU 
#define DEBUG_GXZHI
#endif
