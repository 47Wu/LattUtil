#include "input.hpp"

void inputdataset::parse_opts()
{
	// parse options of the program
	using namespace std;
	using namespace argutil;
	using namespace ioutil;

	string line;
	vector<string> SPline;

	// -h: print usage
	if (cmdOptionExists(argv, argv + argc, "-h")) {
		print_usage();
		exit(ec_io);
	}

	// -r: regulate structure
	if (cmdOptionExists(argv, argv + argc, "-r")) {
		this->reg_str = true;
	}

	// -p: regulate to primitive
	if (cmdOptionExists(argv, argv + argc, "-p")) {
		this->reg_str = true;
		this->reg_to_primitive = true;
	}

	// -s: regulate to std
	if (cmdOptionExists(argv, argv + argc, "-s") || cmdOptionExists(argv, argv + argc, "--std")) {
		this->reg_str = true;
		this->reg_to_std = true;
	}

	// -f: structure filename
	if (cmdOptionExists(argv, argv + argc, "-f")) {
		this->str_filename = getCmdOption(argv, argv + argc, "-f");		
	}

	// -t: output lattice type
	if (cmdOptionExists(argv, argv + argc, "-t")) {
		string type = getCmdOption(argv, argv + argc, "-t");
		if (type == "f" || type == "F") {
			this->reg_to_latttype = frac;
		}
		else if (type == "c" || type == "C") {
			this->reg_to_latttype = cart;
		}
		else {
			cerr << "inputdataset::parse_opts(): Error: unrecognized lattice type: " << type << endl;
			print_usage();
			exit(ec_io);
		}
	}

	// --sc: create supercell
	if (cmdOptionExists(argv, argv + argc, "--sc")) {
		//
		if (cmdOptionExists(argv, argv + argc, "--rd")) {
			cerr << "inputdataset::parse_opts(): Error: cannot running with both \"--sc\" and \"--rd\" options." << endl;
			print_usage();
			exit(ec_io);
		}
		//
		line = getCmdOption(argv, argv + argc, "--sc");
		split(SPline, line);
		this->T_reconstruct.eye();
		this->T_reconstruct(0, 0) = stoi(SPline[0]);
		this->T_reconstruct(1, 1) = stoi(SPline[1]);
		this->T_reconstruct(2, 2) = stoi(SPline[2]);
		//
	}

	// --rd: redefine lattice
	if (cmdOptionExists(argv, argv + argc, "--rd")) {
		string filename = getCmdOption(argv, argv + argc, "--rd");
		fstream fin;
		fin.open(filename, fstream::in);
		fin >> T_reconstruct(0, 0) >> T_reconstruct(1, 0) >> T_reconstruct(2, 0);
		fin >> T_reconstruct(0, 1) >> T_reconstruct(1, 1) >> T_reconstruct(2, 1);
		fin >> T_reconstruct(0, 2) >> T_reconstruct(1, 2) >> T_reconstruct(2, 2);
		fin.close();
	}
	
	if (!armautil::is_eye(T_reconstruct)) {
		this->reg_to_reconstruct = true;
	}
	else {
		this->reg_to_reconstruct = false;
	}
	
	// --sp: symprec
	if (cmdOptionExists(argv, argv + argc, "--sp")) {
		line = getCmdOption(argv, argv + argc, "--sp");
		this->symprec = stod(line);
	}

	return;
}

void inputdataset::read_pos()
{
	using namespace std;
	using namespace ioutil;

	fstream fin;
	string line;
	vector<string> SPline;
	double factor = 1.0;
	int ispecie = 0;
	int icount = 0;
	char coord_type;

	this->atom_types.clear();
	this->atom_positions.clear();

	fin.open(str_filename, fstream::in);
	if (fin.fail()) {
		cerr << "inputdataset::read_pos(): Error: cannot open \"" << str_filename << "\"." << endl;
		cerr << endl;
		exit(ec_io);
	}
	getline(fin, line);	// first line: dummy info

	getline(fin, line); split(SPline, line);	// second line: factor;
	factor = stod(SPline[0]);

	getline(fin, line);	split(SPline, line); //a-axis
	L(0, 0) = stod(SPline[0]) * factor;
	L(1, 0) = stod(SPline[1]) * factor;
	L(2, 0) = stod(SPline[2]) * factor;

	getline(fin, line);	split(SPline, line); //b-axis
	L(0, 1) = stod(SPline[0]) * factor;
	L(1, 1) = stod(SPline[1]) * factor;
	L(2, 1) = stod(SPline[2]) * factor;

	getline(fin, line);	split(SPline, line); //c-axis
	L(0, 2) = stod(SPline[0]) * factor;
	L(1, 2) = stod(SPline[1]) * factor;
	L(2, 2) = stod(SPline[2]) * factor;

	getline(fin, line);	split(this->atom_species, line); //atom species
	this->natomtype = this->atom_species.size();

	// check repeated name of species
	for (auto it1 = this->atom_species.begin(); it1 < this->atom_species.end(); it1++) {
		for (auto it2 = this->atom_species.begin(); it2 < it1; it2++) {
			if (*it2 == *it1) {
				cout << "inputdataset::read_pos(): Error: same species name \"" << *it1 << "\" ";
				cout << "for species [" << std::distance(this->atom_species.begin(), it2) + 1 << "] and [" << std::distance(this->atom_species.begin(), it1) + 1 << "]." << endl;
				exit(ec_io);
			}
		}
	}

	getline(fin, line); split(SPline, line);	//num of each species
	this->natom = 0; icount = 0; ispecie = 0;
	//
	for (auto it = SPline.begin(); it < SPline.end(); it++) {
		ispecie++;
		int nspecie = stoi(*it);
		for (int ii = 0; ii < nspecie; ii++) {
			this->atom_types.push_back(ispecie);
			icount++;
		}
		this->natom += nspecie;
	}

	if (this->natom > MAXNATOM) {
		cout << "inputdataset::read_pos(): Error: MAXNATOM exceeded." << endl;
		cout << "   Current MAXNATOM: " << MAXNATOM << ", Needed value: " << natom << endl;
		cout << "   You may define a larger MAXNATOM and recompile LattUtil to workaround this error." << endl;
		exit(ec_io);
	}

	getline(fin, line); split(SPline, line);	//Direct or Cart
	coord_type = SPline[0][0];

	if (coord_type == 'C' || coord_type == 'c') {
		//
		arma::dcolvec3 r_cart, r_frac;
		//
		for (int isite = 0; isite < natom; isite++) {
			getline(fin, line); split(SPline, line);	//Atom coordinates
			r_cart(0) = stod(SPline[0]);
			r_cart(1) = stod(SPline[1]);
			r_cart(2) = stod(SPline[2]);
			r_frac = L.i() * r_cart;
			this->atom_positions.push_back(r_frac);
		}
	}
	//
	else if (coord_type == 'D' || coord_type == 'd') {
		//
		arma::dcolvec3  r_frac;
		//
		for (int ii = 0; ii < this->natom; ii++) {	//Atom coordinates
			getline(fin, line); split(SPline, line);
			r_frac(0) = stod(SPline[0]);
			r_frac(1) = stod(SPline[1]);
			r_frac(2) = stod(SPline[2]);
			this->atom_positions.push_back(r_frac);
		}
	}
	//
	else {
		cout << "inputdataset::read_pos(): Error: Bad POSCAR - unrecognized coord_type: " << coord_type << "." << endl;
		exit(ec_io);
	}

	fin.close();
}

void inputdataset::print_settings()
{
	using namespace std;
	
	cout << "\n-------------------  SETTINGS  -----------------------\n\n";
	cout << "              regulate structure: " << (reg_str ? "true" : "false") << endl;
	cout << " regulate structure to primitive: " << (reg_to_primitive ? "true" : "false") << endl;
	cout << "             reconstruct lattice: " << (reg_to_reconstruct ? "true" : "false") << endl;
	
	if (reg_to_reconstruct) {
		cout << endl << " Redefined Lattice:" << endl;
		cout << "A1 = " << setw(3) << T_reconstruct(0, 0) << " " << setw(3) << T_reconstruct(1, 0) << " " << setw(3) << T_reconstruct(2, 0) << endl;
		cout << "A2 = " << setw(3) << T_reconstruct(0, 1) << " " << setw(3) << T_reconstruct(1, 1) << " " << setw(3) << T_reconstruct(2, 1) << endl;
		cout << "A3 = " << setw(3) << T_reconstruct(0, 2) << " " << setw(3) << T_reconstruct(1, 2) << " " << setw(3) << T_reconstruct(2, 2) << endl;
	}

	cout << endl;

	return;
}

void inputdataset::print_struct()
{
	using namespace std;

	cout << "\n-------------------  STRUCTURE  -----------------------\n\n";

	cout << "Lattice parameter:" << endl;
	for (int i = 0; i < 3; i++) {
		cout << "a_" << i + 1 << ":  ";
		cout << left << fixed << setprecision(6) << setw(10) << L(0,i);
		cout << left << fixed << setprecision(6) << setw(10) << L(1,i);
		cout << left << fixed << setprecision(6) << setw(10) << L(2,i) << endl;
	}

	cout << "\nAtomic positions:" << endl;
	for (int i = 0; i < natom; i++) {
		cout << right << setw(5) << i + 1 << "-" << left << setw(5) << (atom_species[atom_types[i] - 1] + ":");
		cout << right << fixed << setprecision(6) << setw(10) << atom_positions[i](0) << " ";
		cout << right << fixed << setprecision(6) << setw(10) << atom_positions[i](1) << " ";
		cout << right << fixed << setprecision(6) << setw(10) << atom_positions[i](2) << endl;
	}

	cout << defaultfloat;
	cout << endl;

	return;
}

void inputdataset::print_usage()
{
	using namespace std;

	cout << endl;
	cout << "Usage:" << endl;
	cout << "  LattUtil.x -r -p -s -f <filename> -t <type> --sc <supercell index> --rd <lattice redef file> --sp <symprec>" << endl;
	cout << "Options:" << endl;
	cout << "   -h: help, print this message" << endl;
	cout << "   -r: regulate structure" << endl;
	cout << "   -p: regulate to primitive" << endl;
	cout << "   -s: --std, regulate to spg standardized cell" << endl;
	cout << "   -f: filename, the next argument is the filename of structure file (i.e., POSCAR)" << endl;
	cout << "   -t: type of output lattice, the next argument is \"c\" for cart, \"f\" for frac, the default is f" << endl;
	cout << " --sc: create supercell, e.g., --sc \"2 2 1\"" << endl;
	cout << " --rd: redefine lattice, the next argument is the filename of redefined lattice" << endl
		<< "       which contains an matrix of <A11 A12 A13> <A21 A22 A23> <A31 A32 A33>" << endl
		<< "       here new lattice vector a1_new = A11*a1_old + A12*a2_old + A13*a3_old" << endl
		<< "       all A_ij's must be integer" << endl;
	cout << " --sp: symprec, the next argument is a small number specifying precision for symmetry finding," << endl
		 << "       the default is 1E-2" << endl;
	cout << endl;

	return;
}
