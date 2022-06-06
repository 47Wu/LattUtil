#pragma once
#include "input.hpp"

extern "C" {
#include "spglib.h"
}

void find_equiv_sites(const arma::dmat33 &P, const arma::dcolvec3 &os, const int natom_out, const int natom_in,
	                  std::vector<int> &types_out, std::vector<arma::dcolvec3> &positions_out, 
	                  const std::vector<int> &types_in, const std::vector<arma::dcolvec3> &positions_in,
	                  const double prec);

class str_regulator 
{
public:
	SpglibDataset *spgdataset;

	// settings
	const bool reg_str;                           // regulate structure?
	const bool reg_to_primitive;                  // convert to primitive?
	const bool reg_to_std;                        // convert to spg standardized cell?
	const latttype reg_to_latttype;               // output lattice type (frac or cart)
	const double symprec;

	// regulated structure
	const int natom;
	arma::dmat33 L;                                  // regulated lattice matrix
	const std::vector<std::string> atom_species;     // atom species
	std::vector<int> atom_types;                     // start from 1
	std::vector<arma::dcolvec3> atom_positions;      // fraction coordinates of atom_positions
	std::vector<int> natom_each_type;

	// constructor & destructor
	str_regulator(const inputdataset * LUinput_ptr) :
		LUinput_ptr(LUinput_ptr),
		spgdataset(nullptr),
		reg_str(LUinput_ptr->reg_str),
		reg_to_primitive(LUinput_ptr->reg_to_primitive),
		reg_to_std(LUinput_ptr->reg_to_std),
		reg_to_latttype(LUinput_ptr->reg_to_latttype),
		symprec(LUinput_ptr->symprec),
		natom(LUinput_ptr->natom),
		L(LUinput_ptr->L),
		atom_species(LUinput_ptr->atom_species),
		atom_types(LUinput_ptr->atom_types),
		atom_positions(LUinput_ptr->atom_positions),
		natom_each_type(atom_species.size(), 0)
	{
		for (int i = 0; i < natom; i++) {
			natom_each_type[atom_types[i] - 1]++;
		}
	}
	str_regulator() = delete;
	~str_regulator() { if (spgdataset != nullptr) spg_free_dataset(spgdataset); }

	void get_reg_struct();
	void print_symmetry() const;
	void write_reg_struct() const;

private:
	const inputdataset *LUinput_ptr;
};
