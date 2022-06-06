#pragma once
#include <iomanip>
#include "utils.hpp"

class inputdataset
{
public:
	int argc;
	char **argv;
	std::string str_filename;

	// tags for workflow control
	bool reg_str;                           // regulate structure?
	bool reg_to_primitive;                  // convert to primitive?
	bool reg_to_std;                        // convert to spg standardized cell?
	bool reg_to_reconstruct;                // reconstruct lattice?
	latttype reg_to_latttype;               // output lattice type (frac or cart)
	double symprec;                         // precision for finding symmetry

	// struct info
	int natom;
	int natomtype;
	arma::dmat33 L;                                  // lattice matrix (in spglib convention)
	std::vector<std::string> atom_species;           // atom species
	std::vector<int> atom_types;                     // start from 1
	std::vector<arma::dcolvec3> atom_positions;      // fraction coordinates of atom_positions

	// other info
	arma::imat33 T_reconstruct;             // transform matrix for lattice reconstruct

	// constructor & destructor
	inputdataset(int argc, char ** argv) :
		argc(argc),
		argv(argv),
		str_filename("POSCAR"),
		reg_str(false),
		reg_to_primitive(false),
		reg_to_std(false),
		reg_to_reconstruct(false),
		reg_to_latttype(frac),
		symprec(1E-2),
		natom(0),
		natomtype(0),
		L(arma::fill::eye),
		atom_species(),
		atom_types(),
		atom_positions(),
		T_reconstruct(arma::fill::eye),
		tag_lines()
	{}
	inputdataset() = delete;
	~inputdataset(){}

	// functions
	void parse_opts();
	void read_pos();
	void print_settings();
	void print_struct();

private:
	std::vector<std::string> tag_lines;

	void print_usage();
};
