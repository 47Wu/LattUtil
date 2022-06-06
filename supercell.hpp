#pragma once
#include "input.hpp"
#include "regstruct.hpp"

class supercell
{
public:
	const arma::imat33 T_reconstruct;

	// structure data for redefined lattice
	int natom_sc;
	arma::dmat33 L_sc;
	const std::vector<std::string> atom_species;
	std::vector<int> atom_types_sc;
	std::vector<arma::dcolvec3> atom_positions_sc;
	std::vector<int> natom_each_type_sc;

	// constructor & destructor
	supercell(const inputdataset * LUinput_ptr, const str_regulator * ST_reg_ptr) :
		LUinput_ptr(LUinput_ptr),
		ST_reg_ptr(ST_reg_ptr),
		T_reconstruct(LUinput_ptr->T_reconstruct),
		natom_sc(0),
		L_sc(arma::fill::eye),
		atom_species(LUinput_ptr->atom_species),
		atom_types_sc(),
		atom_positions_sc(),
		natom_each_type_sc(atom_species.size(), 0)
	{}
	supercell() = delete;
	~supercell(){}

	void get_supercell();
	void write_supercell() const;

private:
	const inputdataset *LUinput_ptr;
	const str_regulator *ST_reg_ptr;
};
