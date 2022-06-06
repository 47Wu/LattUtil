#include "supercell.hpp"

void supercell::get_supercell()
{
	const int natom_reg = ST_reg_ptr->natom;
	const arma::dmat33 &L_reg = ST_reg_ptr->L;
	const double V_reg = arma::det(L_reg);
	const std::vector<int> &types_reg = ST_reg_ptr->atom_types;
	const std::vector<arma::dcolvec3> &positions_reg = ST_reg_ptr->atom_positions;

	arma::dmat33 P_inv; for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) { P_inv(i, j) = T_reconstruct(i, j); } }
	const arma::dmat33 P = P_inv.i();
	const arma::dcolvec3 os = { 0.0, 0.0, 0.0 };

	const double symprec = LUinput_ptr->symprec;

	// get L_sc and natom_sc
	L_sc = L_reg * P_inv;
	const double V_sc = arma::det(L_sc);
	const int sc_multp = (int)std::round(V_sc / V_reg);
	natom_sc = natom_reg * sc_multp;
	
	// get atom types and atom positions
	find_equiv_sites(P, os, natom_sc, natom_reg, atom_types_sc, atom_positions_sc, types_reg, positions_reg, symprec);

	// get natom_each_type_sc
	for (const atype : atom_types_sc) {
		natom_each_type_sc[atype - 1]++;
	}

	return;
}

void supercell::write_supercell() const
{
	using namespace std;

	string filename;
	fstream fout;
	arma::dcolvec3 x;
	const int natomtype = atom_species.size();

	filename = "POSCAR";
	filename += LUinput_ptr->reg_str ? "_reg" : "";
	filename += "_reconstruct";
	fout.open(filename, fstream::out);
	//
	fout << "# Reconstructed Struct from LattUtil" << endl;
	fout << "    1.0" << endl;
	//
	for (int i = 0; i < 3; i++) {
		fout << right << fixed << setprecision(16) << setw(23) << this->L_sc(0, i);
		fout << right << fixed << setprecision(16) << setw(23) << this->L_sc(1, i);
		fout << right << fixed << setprecision(16) << setw(23) << this->L_sc(2, i);
		fout << endl;
	}
	//
	for (const string &species : this->atom_species) {
		fout << right << setw(5) << species;
	}
	fout << endl;
	//
	for (const int n : this->natom_each_type_sc) {
		fout << right << setw(5) << n;
	}
	fout << endl;
	//
	if (LUinput_ptr->reg_to_latttype == cart) {
		fout << "Cartesian" << endl;
	}
	else {
		fout << "Direct" << endl;
	}
	//
	for (int atype = 1; atype <= natomtype; atype++) {
		for (int i = 0; i < natom_sc; i++) {
			if (atom_types_sc[i] == atype) {
				if (LUinput_ptr->reg_to_latttype == cart) {  // cart
					x = L_sc * atom_positions_sc[i];
				}
				else {                          // frac
					x = atom_positions_sc[i];
				}
				fout << right << fixed << setprecision(16) << setw(20) << x(0);
				fout << right << fixed << setprecision(16) << setw(20) << x(1);
				fout << right << fixed << setprecision(16) << setw(20) << x(2);
				fout << endl;
			}
		}
	}
	//
	fout.close();
	cout << "Structure reconstructed to " << filename << "." << endl;

	return;
}
