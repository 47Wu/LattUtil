#include "regstruct.hpp"

void str_regulator::get_reg_struct()
{
	using namespace std;
	using namespace armautil;

	// the regulation of the structure is done by back-transforming the idealized structure
	double lattice[3][3];
	double position[MAXNATOM][3];
	int types[MAXNATOM];
	//
	dmat33_to_std33(lattice, L);
	//
	for (size_t i = 0; i < this->atom_positions.size(); i++) {
		for (int j = 0; j < 3; j++) {
			position[i][j] = this->atom_positions[i](j);
		}
		types[i] = this->atom_types[i];
	}
	//
	this->spgdataset = spg_get_dataset(lattice, position, types, natom, symprec);

	// by construct, the structure info is not regulated
	if (!reg_str) {
		return;
	}

	// get transformation matrix and origin shift (see spglib doc)
	arma::dmat33 P;
	arma::dmat33 R;
	arma::dcolvec3 os;
	//
	std33_to_dmat33(P, spgdataset->transformation_matrix);
	std33_to_dmat33(R, spgdataset->std_rotation_matrix);
	for (int i = 0; i < 3; i++) { os(i) = spgdataset->origin_shift[i]; }

#ifdef __DEBUG
	cout << " transformation matrix:" << endl;
	cout << P << endl;
	cout << " origin shift:" << endl;
	cout << "(" << os(0) << "  " << os(1) << "  " << os(2) << ")" << endl << endl;
	cout << " rotation induced by idealize:" << endl;
	cout << R;
	cout << " det(R) = " << arma::det(R) << endl;
	//cout << " mapping table to primitive:" << endl;
	//for (int i = 0; i < natom; i++) {
	//	cout << i << "->" << spgdataset->mapping_to_primitive[i] << endl;
	//}
	cout << endl;
#endif
	
	// get idealized structure
	int natom_std;
	natom_std = spg_standardize_cell(lattice, position, types, natom, /* to_pri */ 0, /* no_ideal */ 0, symprec);
	
	// update L atom_types atom_positions
	std33_to_dmat33(L, lattice);
	L = R.i() * L * P;
	//
	// we want to get x_origin from x_std, therefore we need transformation:
	//    x_origin = P_inv * x_std + P_inv * os
	// in: std struct, out: origin struct
	std::vector<int> types_std(natom_std);
	std::vector<arma::dcolvec3> positions_std(natom_std);
	for (int i = 0; i < natom_std; i++) {
		types_std[i] = types[i];
		positions_std[i](0) = position[i][0];
		positions_std[i](1) = position[i][1];
		positions_std[i](2) = position[i][2];
	}
	find_equiv_sites(P.i(), P.i()*os, natom, natom_std, this->atom_types, this->atom_positions, types_std, positions_std, symprec);

	// Now the structure members of (*this) has been regulated.

	return;
}

void str_regulator::print_symmetry() const
{
	using namespace std;

	cout << "\n\n-------------------  SYMMETRY  -----------------------\n\n";
	cout << "Hall Symbol: " << this->spgdataset->hall_symbol << endl;
	cout << "International Symbol: " << this->spgdataset->international_symbol << endl;
	cout << "International Number: " << this->spgdataset->spacegroup_number << endl;
	cout << "Point Group: " << this->spgdataset->pointgroup_symbol << endl;
	cout << endl;

	return;
}

void str_regulator::write_reg_struct() const
{
	using namespace std;

	if (!reg_to_primitive && !reg_to_std && !reg_str) { return; }	// nothing will be done if no reg_flags are defined

	string filename;
	fstream fout;
	arma::dcolvec3 x;
	const int natomtype = atom_species.size();

	if (reg_str) {
		filename = "POSCAR_reg_origin";
		fout.open(filename, fstream::out);
		//
		fout << "# Regulated Struct from LattUtil (origin)" << endl;
		fout << "    1.0" << endl;
		//
		for (int i = 0; i < 3; i++) {
			fout << right << fixed << setprecision(16) << setw(23) << this->L(0, i);
			fout << right << fixed << setprecision(16) << setw(23) << this->L(1, i);
			fout << right << fixed << setprecision(16) << setw(23) << this->L(2, i);
			fout << endl;
		}
		//
		for (const string &species : this->atom_species) {
			fout << right << setw(5) << species;
		}
		fout << endl;
		//
		for (const int n : this->natom_each_type) {
			fout << right << setw(5) << n;
		}
		fout << endl;
		//
		if (reg_to_latttype == cart) {
			fout << "Cartesian" << endl;
		}
		else {
			fout << "Direct" << endl;
		}
		//
		for (int atype = 1; atype <= natomtype; atype++) {
			for (int i = 0; i < natom; i++) {
				if (atom_types[i] == atype) {
					if (reg_to_latttype == cart) {  // cart
						x = L * atom_positions[i];
					}
					else {                          // frac
						x = atom_positions[i];
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
		cout << "Structure regulated to " << filename << "." << endl;
	}

	if (reg_to_primitive) {
		double latt_pri[3][3];
		double pos_pri[MAXNATOM][3];
		int types_pri[MAXNATOM];
		int natom_pri;
		arma::dmat33 L_pri;

		armautil::dmat33_to_std33(latt_pri, L);
		for (int i = 0; i < natom; i++) {
			types_pri[i] = atom_types[i];
			pos_pri[i][0] = atom_positions[i](0);
			pos_pri[i][1] = atom_positions[i](1);
			pos_pri[i][2] = atom_positions[i](2);
		}

		natom_pri = spg_standardize_cell(latt_pri, pos_pri, types_pri, natom, /* to_primitive = */ 1, /* no_idealize = */ 0, symprec);
		armautil::std33_to_dmat33(L_pri, latt_pri);

		// output to file
		filename = "POSCAR_reg_pri";
		fout.open(filename, fstream::out);
		//
		fout << "# Regulated Struct from LattUtil (primitive)" << endl;
		fout << "    1.0" << endl;
		//
		for (int i = 0; i < 3; i++) {
			fout << right << fixed << setprecision(16) << setw(23) << latt_pri[0][i];
			fout << right << fixed << setprecision(16) << setw(23) << latt_pri[1][i];
			fout << right << fixed << setprecision(16) << setw(23) << latt_pri[2][i];
			fout << endl;
		}
		//
		for (const string &species : this->atom_species) {
			fout << right << setw(5) << species;
		}
		fout << endl;
		//
		for (const int n : this->natom_each_type) {
			fout << right << setw(5) << (n * natom_pri) / natom;
		}
		fout << endl;
		//
		if (reg_to_latttype == cart) {
			fout << "Cartesian" << endl;
		}
		else {
			fout << "Direct" << endl;
		}
		//
		for (int atype = 1; atype <= natomtype; atype++) {
			for (int i = 0; i < natom; i++) {
				if (types_pri[i] == atype) {
					if (reg_to_latttype == cart) {  // cart
						x(0) = pos_pri[i][0]; x(1) = pos_pri[i][1]; x(0) = pos_pri[i][2];
						x = L_pri * x;
					}
					else {                          // frac
						x(0) = pos_pri[i][0]; x(1) = pos_pri[i][1]; x(0) = pos_pri[i][2];
					}
					fout << right << fixed << setprecision(16) << setw(20) << x(0);
					fout << right << fixed << setprecision(16) << setw(20) << x(1);
					fout << right << fixed << setprecision(16) << setw(20) << x(2);
					fout << endl;
				}
			}
		}
		fout.close();
		cout << "Structure regulated to " << filename << "." << endl;
	}

	if (reg_to_std) {
		double latt_std[3][3];
		double pos_std[MAXNATOM][3];
		int types_std[MAXNATOM];
		int natom_std;
		arma::dmat33 L_std;

		armautil::dmat33_to_std33(latt_std, L);
		for (int i = 0; i < natom; i++) {
			types_std[i] = atom_types[i];
			pos_std[i][0] = atom_positions[i](0);
			pos_std[i][1] = atom_positions[i](1);
			pos_std[i][2] = atom_positions[i](2);
		}

		natom_std = spg_standardize_cell(latt_std, pos_std, types_std, natom, /* to_primitive = */ 0, /* no_idealize = */ 0, symprec);
		armautil::std33_to_dmat33(L_std, latt_std);

		// output to file
		filename = "POSCAR_reg_std";
		fout.open(filename, fstream::out);
		//
		fout << "# Regulated Struct from LattUtil (spg stdard)" << endl;
		fout << "    1.0" << endl;
		//
		for (int i = 0; i < 3; i++) {
			fout << right << fixed << setprecision(16) << setw(23) << latt_std[0][i];
			fout << right << fixed << setprecision(16) << setw(23) << latt_std[1][i];
			fout << right << fixed << setprecision(16) << setw(23) << latt_std[2][i];
			fout << endl;
		}
		//
		for (const string &species : this->atom_species) {
			fout << right << setw(5) << species;
		}
		fout << endl;
		//
		for (const int n : this->natom_each_type) {
			fout << right << setw(5) << (n * natom_std) / natom;
		}
		fout << endl;
		//
		if (reg_to_latttype == cart) {
			fout << "Cartesian" << endl;
		}
		else {
			fout << "Direct" << endl;
		}
		//
		for (int atype = 1; atype <= natomtype; atype++) {
			for (int i = 0; i < natom; i++) {
				if (types_std[i] == atype) {
					if (reg_to_latttype == cart) {  // cart
						x(0) = pos_std[i][0]; x(1) = pos_std[i][1]; x(0) = pos_std[i][2];
						x = L_std * x;
					}
					else {                          // frac
						x(0) = pos_std[i][0]; x(1) = pos_std[i][1]; x(0) = pos_std[i][2];
					}
					fout << right << fixed << setprecision(16) << setw(20) << x(0);
					fout << right << fixed << setprecision(16) << setw(20) << x(1);
					fout << right << fixed << setprecision(16) << setw(20) << x(2);
					fout << endl;
				}
			}
		}
		fout.close();
		cout << "Structure regulated to " << filename << "." << endl;
	}

	return;
}

void find_equiv_sites(const arma::dmat33 & P, const arma::dcolvec3 & os, const int natom_out, const int natom_in, 
	                  std::vector<int>& types_out, std::vector<arma::dcolvec3>& positions_out, 
	                  const std::vector<int>& types_in, const std::vector<arma::dcolvec3>& positions_in,
	                  const double prec)
{
	/* 
	 * Given the transformation_matrix P and origin_shift os that defines:
	 *        x_out = P * x_in + os
	 * For all x_in given by position_in, find all equivalent x_out and fill them into x_out
	 * Note that x_out may NOT be always one-to-one mapped by x_in, we need to set up an supercell for x_out finding.
	 * The supercell indices are determined by:
	 *    note that the final x_out should be in [0, 1), we can thus have:
	 *        x_in = P^(-1) * (x_out - os)
	 *    the range for each component of x_in could be determined as:
	 *    e.g.:   x_in[0]_min >= sum { min ( P^(-1)_[0][i] * (x_out[i] - os[i]) ) }
	 */

	const arma::dmat33 P_inv = P.i();

	arma::dcolvec3 x_in, x_out, R_in;
	arma::icolvec3 sc_index_min, sc_index_max;

	types_out.clear();
	positions_out.clear();

	// find sc_index_min/max
	for (int i = 0; i < 3; i++) {
		double x_in_min = 0.0, x_in_max = 0.0;
		for (int j = 0; j < 3; j++) {
			x_in_min += std::min(P_inv(i, j)*(-os(j)), P_inv(i, j)*(1 - os(j)));
			x_in_max += std::max(P_inv(i, j)*(-os(j)), P_inv(i, j)*(1 - os(j)));
		}
		sc_index_min(i) = (int)std::floor(x_in_min);
		sc_index_max(i) = (int)std::floor(x_in_max);
	}

	// the site coordinate under redefined lattice and origin lattice is not simply 1-to-n corresponded
	for (int i = 0; i < natom_in; i++) {
		//
		for (int R0 = sc_index_min(0); R0 <= sc_index_max(0); R0++) {
			for (int R1 = sc_index_min(1); R1 <= sc_index_max(1); R1++) {
				for (int R2 = sc_index_min(2); R2 <= sc_index_max(2); R2++) {
					R_in(0) = R0; R_in(1) = R1; R_in(2) = R2;
					//
					x_in = positions_in[i] + R_in;
					x_out = P * x_in + os;
					mathutil::mod_1(x_out);
					//
					if (armautil::find_dcolvec3_from_list(x_out, positions_out, prec) < 0) {	// new site
						positions_out.push_back(x_out);
						types_out.push_back(types_in[i]);
					}
					//
				}
			}
		}
		//
	}

	// check natom_out
	if (natom_out != positions_out.size()) {
		std::cerr << "find_equiv_sites(): Error: find " << positions_out.size() << "sites, it ought to be " << natom_out << std::endl;
		exit(ec_math);
	}
	
	return;
}
