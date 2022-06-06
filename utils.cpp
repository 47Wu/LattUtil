#include "utils.hpp"

#ifdef __DEBUG
// employ some debuggers
void debug_ccxu() {
	std::cerr << "ccxu: Hello there #" << ccxu_count << std::endl;
	ccxu_count++;
}

void debug_gxzhi() {
	std::cerr << "gxzhi: Hello there #" << gxzhi_count << std::endl;
	gxzhi_count++;
}
#endif


void print_version()
{
	using namespace std;

	cout << endl;
	cout << " | ------------------------------------------------------------------------------------------------ |" << endl;
	cout << " |              ======================== LattUtil v" 
		                                          << LUversion << " ========================                     |" << endl;
	cout << " |                                                                                                  |" << endl;
	cout << " | This is a utility program for handling VASP structure file (i.e. POSCAR).                        |" << endl;
	cout << " |                                                                                                  |" << endl;
	cout << " | List of capabilities:                                                                            |" << endl;
	cout << " | 1. Find and display the symmetry of the system (info of space group, etc.)                       |" << endl;
	cout << " | 2. Regulate the structure with a fine precision                                                  |" << endl;
	cout << " |    (this is sometimes useful since the output of some software is only in single float prec).    |" << endl;
	cout << " | 3. Re-define the lattice (you may use this function to build a supercell)                        |" << endl;
	cout << " |                                                                                                  |" << endl;
	cout << " | You may use the option \"-h\" to show usage.                                                       |" << endl;
	cout << " |                                                                                                  |" << endl;
	cout << " |                               Developed by WSQ in May 2022.                                      |" << endl;
	cout << " | ------------------------------------------------------------------------------------------------ |" << endl;
}

// some template specializations
int notzero(const dcomplex &val, const double threshold) {
	if ((threshold < val.real()) + (val.real() < -threshold)) {
		return 1;
	}
	else if ((threshold < val.imag()) + (val.imag() < -threshold)) {
		return 1;
	}
	else {
		return 0;
	}
}

void roundzero(dcomplex &var, const double threshold) {
	if (!notzero<double>(var.real(), threshold)) {
		var.real(0.0);
	}
	if (!notzero<double>(var.imag(), threshold)) {
		var.imag(0.0);
	}
	return;
}



/*-----------------------------------------------*/
/*---------- util func implementations ----------*/
/*-----------------------------------------------*/

void ioutil::split(std::vector<std::string> &VS, const std::string &S) {
	using namespace std;
	string tmp_str;
	auto it = S.begin();

	VS.clear();
	tmp_str.clear();

	do {
		while (isspace(*it) && it != S.end()) { ++it; }
		while (!isspace(*it) && it != S.end()) {
			tmp_str.push_back(*it);
			++it;
		}
		if (!tmp_str.empty()) { VS.push_back(tmp_str); }
		tmp_str.clear();
	} while (it != S.end());

	return;
}

void ioutil::split(std::vector<std::string>& VS, const std::string & S, const char sp)
{
	using namespace std;
	string tmp_str;
	auto it = S.begin();

	VS.clear();
	tmp_str.clear();

	do {
		while (*it == sp && it != S.end()) { ++it; }
		while (*it != sp && it != S.end()) {
			tmp_str.push_back(*it);
			++it;
		}
		if (!tmp_str.empty()) { VS.push_back(tmp_str); }
		tmp_str.clear();
	} while (it != S.end());

	return;
}

void ioutil::eq2space(std::string &S)
{
	std::string tmp_str;
	tmp_str = std::regex_replace(S, std::regex("="), " ");
	S = tmp_str;
	return;
}

void ioutil::qt2space(std::string &S)
{
	std::string tmp_str;
	tmp_str = std::regex_replace(S, std::regex("['\"]"), " ");
	S = tmp_str;
	return;
}

void ioutil::replacea2b(std::string &S, const std::string &a, const std::string &b)
{
	std::string tmp_str, reg_exp;
	reg_exp = "[" + a + "]";
	tmp_str = std::regex_replace(S, std::regex(reg_exp), b);
	S = tmp_str;
	return;
}

void ioutil::trim_comment(std::string & S)
{
	std::string tmp_str;
	std::regex comment("[#!].*$");
	tmp_str = std::regex_replace(S, comment, "");
	S = tmp_str;
	return;
}

void ioutil::trim_space(std::string & S)
{
	if (S.empty())return;

	using namespace std;
	string tmp_str;
	long unsigned int i1 = 0;
	long unsigned int i2 = S.length() - 1;

	tmp_str.clear();

	while (i1 < i2) {
		if (isspace(S[i1]))i1++;
		else break;
	}	// i1 -> 1st non-space character in S
	while (i2 > i1) {
		if (isspace(S[i2]))i2--;
		else break;
	}	// i2 -> last non-space character in S

	for (auto i = i1; i <= i2; i++) {
		tmp_str.push_back(S[i]);
	}

	S = tmp_str;

	return;
}

void ioutil::cap2low(std::string & S)
{
	const char A2a = 'A' - 'a';

	for (auto it = S.begin(); it < S.end(); it++) {
		if ((*it) >= 'A' && (*it) <= 'Z') {
			*it -= A2a;
		}
	}

	return;
}

bool ioutil::stobool(const std::string & S)
{
	std::string tmp_str = S;
	cap2low(tmp_str);
	replacea2b(tmp_str, ".", "");
	trim_space(tmp_str);
	switch (tmp_str[0])
	{
	case 't':
	case '1':
		return true;
	case 'f':
	case '0':
		return false;
	default:
		std::cout << "ioutil::stobool(): Error: cannot recognize boolean input: " << S << std::endl;
		exit(2);
	}

	return false;
}

dcomplex ioutil::stocx(const std::string &S, int n)
{
	// the valid string value is either a number or "e^m" or "-e^m"

	int m;
	double phase;
	dcomplex val(0.0, 0.0);

	if (S[0] == 'e' || S[1] == 'e') {
		if (S[0] == '-') { m = '0' - S[3]; }
		else { m = S[2] - '0'; }

		phase = 2.0 * pi * ((double)m) / ((double)n);
		val = std::polar(1.0, phase);	// it is said that std::polar can be 4.5x faster than std::exp in vectorized loops
	}

	else {
		val.real(std::stod(S));
	}

	roundzero(val);
	return val;
}

bool ioutil::read_index_from_str(const std::string & S, int & ibegin, int & iend)
{
	std::string tmp_S = S;
	std::vector<std::string> SPline;
	bool hyphen_flag;
	hyphen_flag = std::regex_match(tmp_S, std::regex("^.*[-].*$"));
	
	ioutil::replacea2b(tmp_S, "-", " ");
	ioutil::split(SPline, tmp_S);

	if (hyphen_flag && SPline.size() == 2) {
		ibegin = std::stoi(SPline[0]);
		iend = std::stoi(SPline[1]);
		return true;
	}

	else if (!hyphen_flag && SPline.size() == 1) {
		ibegin = std::stoi(SPline[0]);
		iend = ibegin;
		return true;
	}
	
	ibegin = -1; iend = -1;
	return false;
}

double ioutil::simpexptod(const std::string & Exp)
{
	/*
	 * return value for simple input expression
	 * supported constant: "pi", "sqrt$a" ($a is a number)
	 * supported operations: + - * /
	 * e.g. Exp = "1/5 * sqrt3"
	 * braket is not supported
	 */

	using namespace std;
	using namespace ioutil;

	string tmp_Exp = Exp;
	vector<string> SPexp;
	vector<char> ope_stack;
	vector<double> val_stack;
	double res, val_new;
	bool isval;

	replacea2b(tmp_Exp, "+", " + ");
	replacea2b(tmp_Exp, "-", " - ");
	replacea2b(tmp_Exp, "*", " * ");
	replacea2b(tmp_Exp, "/", " / ");

	split(SPexp, tmp_Exp);

	ope_stack.clear();
	val_stack.clear();

	for (const string term : SPexp) {
		switch (term[0]){
		case '+':
		case '-':
		case '*':
		case '/':
			ope_stack.push_back(term[0]);
			isval = false;
			break;
		case 'p':	// pi
			val_new = pi;
			isval = true;
			break;
		case 's':	// sqrt$a
			val_new = stod(term.substr(4));
			val_new = sqrt(val_new);
			isval = true;
			break;
		default:
			val_new = stod(term);
			isval = true;
			break;
		}
		if (isval) {
			if (ope_stack.empty()) {
				val_stack.push_back(val_new);
			}
			else if (ope_stack.back() == '+' || ope_stack.back() == '-') {
				val_stack.push_back(val_new);
			}
			else if (ope_stack.back() == '*') {
				val_stack.back() *= val_new;
				ope_stack.pop_back();
			}
			else if (ope_stack.back() == '/') {
				val_stack.back() /= val_new;
				ope_stack.pop_back();
			}
		}
	}

	vector<char>::const_iterator itr_ope;
	vector<double>::const_iterator itr_val;

	if (ope_stack.size() == val_stack.size()) {
		res = 0.0;
		itr_ope = ope_stack.begin();
		itr_val = val_stack.begin();
	}
	else if (ope_stack.size() == val_stack.size() - 1) {
		res = val_stack.front();
		itr_ope = ope_stack.begin();
		itr_val = next(val_stack.begin());
	}
	else {
		cout << "ioutil::expressiontod(): Error: cannot resolve expression: " << Exp << endl;
		exit(2);
	}

	while (itr_ope != ope_stack.end() && itr_val != val_stack.end()) {
		if (*itr_ope == '+') {
			res += *itr_val;
		}
		else if (*itr_ope == '-') {
			res -= *itr_val;
		}
		itr_ope++;
		itr_val++;
	}

	return res;
}


void mathutil::mod_1(double &val, const double prec)
{
	// in the case that val = 0.999999...
	val = val - std::floor(val + prec);
	return;
}

void mathutil::mod_1(arma::dcolvec3 & V)
{
	mod_1(V(0));
	mod_1(V(1));
	mod_1(V(2));
}

int mathutil::mod_1_re_int(double & val, const double prec)
{
	int ival;
	ival = std::floor(val + prec);
	val = val - ival;
	return ival;
}

int mathutil::get_rand_int(int min, int max)
{
	// put the following code into main() to initialize seed for generator
	// std::srand(std::time(nullptr));
	if (max == min) { return min; }
	return (int)(std::floor(((double)rand()) / ((double)RAND_MAX + 1) * (max + 1 - min))) + min;
}

double mathutil::get_rand_double(double min, double max)
{
	if (max == min) { return min; }
	return ((double)rand()) / ((double)RAND_MAX + 1) * (max - min) + min;
}

arma::cx_colvec mathutil::proj(const arma::cx_colvec u, const arma::cx_colvec v)	// project v on u
{	
	using arma::cdot;
	using arma::norm;

	if (u.n_elem != v.n_elem) {
		std::cout << "mathutil::proj(): Error: projecting v on u with different dimension." << std::endl;
		exit(5);
	}

	return cdot(u, v) / (norm(u)*norm(u)) * u;
}

int mathutil::mod(const int val, const int interval_end, const int interval_begin)
{
	// return the corresponding value of moving val periodically to the interval [begin, end)
	// by default, the function acts as a modulo operation (interval_begin=0)

	if (interval_end <= interval_begin) {
		std::cerr << "mathutil::mod(): Error: bad interval." << std::endl;
		exit(ec_math);
	}

	int res;

	if (interval_begin == 0) {
		res = val % interval_end;
		if (res < 0) {
			res += interval_end;
		}
		return res;
	}
	else {
		return mod(val - interval_begin, interval_end - interval_begin) + interval_begin;
	}
}

char * argutil::getCmdOption(char ** begin, char ** end, const std::string & option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return nullptr;
}

bool argutil::cmdOptionExists(char ** begin, char ** end, const std::string & option)
{
	return std::find(begin, end, option) != end;
}

bool armautil::is_eye(const arma::mat & M)
{
	if (M.n_rows != M.n_cols) {
		return false;
	}

	for (int i = 0; i < M.n_rows; i++) {
		for (int j = 0; j < M.n_cols; i++) {
			if (i == j) {
				if (notzero(M(i, j) - 1)) {
					return false;
				}
			}
			else {
				if (notzero(M(i, j))) {
					return false;
				}
			}
		}
	}

	return true;
}

bool armautil::is_eye(const arma::imat33 & M)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) {
				if (M(i, j) != 1) { return false; }
			}
			else {
				if (M(i, j) != 0) { return false; }
			}
		}
	}

	return true;
}

void armautil::std33_to_dmat33(arma::dmat33 & M, const double M0[3][3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			M(i, j) = M0[i][j];
		}
	}
}

void armautil::dmat33_to_std33(double M[3][3], const arma::dmat33 & M0)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			M[i][j] = M0(i, j);
		}
	}
}

int armautil::find_dcolvec3_from_list(const arma::dcolvec3 & V, const std::vector<arma::dcolvec3>& Vlist, const double prec)
{
	for (int i = 0; i < Vlist.size(); i++) {
		if (!notzero(Vlist[i](0) - V(0), prec) &&
			!notzero(Vlist[i](1) - V(1), prec) &&
			!notzero(Vlist[i](2) - V(2), prec)) 
		{
			return i;
		}
	}

	return -1;	// not found
}
