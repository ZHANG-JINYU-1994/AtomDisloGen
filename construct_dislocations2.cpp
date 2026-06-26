// construct_dislocations2.cpp
//
// Generate or import a cubic crystal lattice and insert user-defined dislocation
// lines or loops into the atomic configuration. The displacement field is
// evaluated from the solid angle subtended by the dislocation surface, and
// additional atoms are inserted near the cut surface when required.
//
// Supported lattice types in the current implementation: sc, bcc, and fcc.
// Supported dislocation geometries include circles, ellipses, regular polygons,
// rectangles, helices, and user-defined polylines.
//
// If you use this code, please cite:
// Jin-Yu Zhang and Wen-Zheng Zhang, Modelling Simul. Mater. Sci. Eng. 27 (2019) 035008.
# include <iostream>
# include <algorithm>
# include <cmath>
# include <fstream>
# include <iomanip>
# include <sstream>
# include <stdlib.h>
# include <stdio.h>
# include <cstring>
//# include <lapacke.h>
//# include <lapacke_config.h>

using namespace std;

//****************************************************************
// Global constants and data structures
const double eps = 1.0e-8;
const double pi = 3.141592653589793;
// Lattice information
struct lattice
// Lattice information
{
	double lattice_para;		//cubic lattice parameter
	char* lattice_type;		//lattice type
	double sizebox[3][2];	//lower and upper bounds of the simulation box
	double Mb2o[3][3];		//orientation matrix from crystal coordinates to box coordinates
	int bbcc_num;				//number of nearest-neighbor vectors
	double(*bre)[3];			//reciprocal vectors defining the Wigner-Seitz cell
	int atomnum;				//number of atoms
	int* points_type;			//atom type
	double(*pointsb)[3];	//atom coordinates
};
lattice lat, lat_disl; //lat stores the full simulation structure; lat_disl stores atoms identified as belonging to dislocation cores

//****************************************************************
// Basic 3D vector operations
double dot(double a[3], double b[3])
// Dot product
{
	double c;
	c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return c;
}

void cross(double a[3], double b[3], double c_cross[3])
// Cross product
{
	c_cross[0] = a[1] * b[2] - a[2] * b[1];
	c_cross[1] = a[2] * b[0] - a[0] * b[2];
	c_cross[2] = a[0] * b[1] - a[1] * b[0];
}

double norm(double a[3])
// Vector norm
{
	double c;
	c = sqrt(dot(a, a));
	return c;
}

void orthogonalization(double xdir[3], double zdir[3])
// Check whether two input directions are nearly orthogonal and orthogonalize the second one
{
	double length_x = norm(xdir), length_z = norm(zdir), dot_xz = dot(xdir, zdir);
	if (length_x < eps || length_z < eps)
	{
		cout << "the input vector should not be zero vector" << endl;
		exit(1);
	}
	else
	{
		if (fabs(dot_xz / (length_x * length_z)) < 0.0175) //if the angle is larger than about 89 degrees, orthogonalize zdir
		{
			for (int i = 0; i < 3; i++)
			{
				zdir[i] -= dot_xz * xdir[i] / (length_x * length_x);
			}
		}
		else //otherwise, the two directions are not sufficiently orthogonal
		{
			cout << "Error: two directions should be orthogonal" << endl;
			exit(1);
		}
	}
}

double solidangle_tri(double a[3], double b[3], double c[3])
// Solid-angle contribution from a triangular surface element
{
	double al, bl, cl, denominator, numerator, s_tri;
	al = norm(a);
	bl = norm(b);
	cl = norm(c);
	denominator = al * bl * cl + al * dot(b, c) + bl * dot(c, a) + cl * dot(a, b);
	double bc_cross[3];
	cross(b, c, bc_cross);
	numerator = dot(bc_cross, a);
	s_tri = atan(numerator / denominator);
	if (denominator < 0)
	{
		if (numerator < 0)
			s_tri -= pi;
		else
			s_tri += pi;
	}
	if (abs(denominator) < eps && abs(numerator) < eps)
	{
		s_tri = 0.5 * pi;
	}
	if (al < eps)
	{
		double temp0 = dot(b, c) / (bl * cl);
		if (temp0 > 1.0)
		{
			s_tri = 0.0;
		}
		else if (temp0 < -1.0)
		{
			s_tri = 0.5 * pi;
		}
		else
		{
			s_tri = 0.5 * acos(temp0);
		}
	}
	if (bl < eps)
	{
		double temp0 = dot(c, a) / (cl * al);
		if (temp0 > 1.0)
		{
			s_tri = 0.0;
		}
		else if (temp0 < -1.0)
		{
			s_tri = 0.5 * pi;
		}
		else
		{
			s_tri = 0.5 * acos(temp0);
		}
	}
	if (cl < eps)
	{
		double temp0 = dot(a, b) / (al * bl);
		if (temp0 > 1.0)
		{
			s_tri = 0.0;
		}
		else if (temp0 < -1.0)
		{
			s_tri = 0.5 * pi;
		}
		else
		{
			s_tri = 0.5 * acos(temp0);
		}
	}
	return s_tri;
}
//****************************************************************
// Debug helper 1: print a vector
void print_vector(double vector[], int m)
// Print a one-dimensional array
{
	printf("the vector is \n");
	for (int i = 0; i < m; i++)
	{
		printf("%f ", vector[i]);
	}
	printf("\n");
	fflush(stdout);
}
// Debug helper 2: print an n-by-3 matrix
void print_matrix(double matrix[][3], int m)
// Print an n-by-3 array
{
	printf("the matrix is \n");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
	fflush(stdout);
}
// Debug helper 3: write atom coordinates and an auxiliary per-atom property
void output_for_test(char filename[], int num_of_point, int point_type[], double point_position[][3], double point_property[])
{
	//write atom types, coordinates, and an auxiliary per-atom property
	int* outbox, atomcount = num_of_point, atomnum_out, * points_type_out, type = 0;
	double(*pointsb_out)[3], * points_property_out;
	outbox = new int[num_of_point];
	memset(outbox, 0, num_of_point * sizeof(int));
	// discard atoms outside the region defined by lat.sizebox
	for (int i = 0; i < num_of_point; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if ((point_position[i][j] < lat.sizebox[j][0] - eps) || (point_position[i][j] > lat.sizebox[j][1] + eps))
			{
				outbox[i] = 1;
				atomcount--;
				break;
			}
		}
	}
	pointsb_out = new double[atomcount][3];
	points_type_out = new int[atomcount];
	points_property_out = new double[atomcount];
	memset(pointsb_out, 0, atomcount * 3 * sizeof(double));
	int count = 0;
	for (int i = 0; i < num_of_point; i++)
	{
		if (outbox[i] != 1)
		{
			if (type < point_type[i] + 1)
			{
				type = point_type[i] + 1;
			}
			points_type_out[count] = point_type[i];
			points_property_out[count] = point_property[i];
			for (int j = 0; j < 3; j++)
			{
				pointsb_out[count][j] = point_position[i][j];
			}
			count++;
		}
	}
	atomnum_out = count;
	delete[] outbox;
	//write the output file
	std::stringstream stream;
	stream << "ITEM: TIMESTEP" << "\n";
	stream << "0" << "\n";
	stream << "ITEM: NUMBER OF ATOMS" << "\n";
	stream << atomnum_out << "\n";
	stream << "ITEM: BOX BOUNDS pp pp pp" << "\n";
	stream << "\t\t" << lat.sizebox[0][0] << "\t" << lat.sizebox[0][1] << "\txlo xhi" << "\n";
	stream << "\t\t" << lat.sizebox[1][0] << "\t" << lat.sizebox[1][1] << "\tylo yhi" << "\n";
	stream << "\t\t" << lat.sizebox[2][0] << "\t" << lat.sizebox[2][1] << "\tzlo zhi" << "\n";
	stream << "ITEM: ATOMS id type x y z property" << "\n";
	// points_type_out stores atom types for visualization; pointsb_out stores atom coordinates; sizebox defines the orthogonal box bounds
	for (int i = 0; i < atomnum_out; i++)
	{
		stream << std::setw(9) << (i + 1) << std::setw(3) << points_type_out[i] + 1 << std::setw(15) << pointsb_out[i][0] << std::setw(15) << pointsb_out[i][1] << std::setw(15) << pointsb_out[i][2] << std::setw(15) << points_property_out[i] << "\n";
	}
	std::string str;
	str = stream.str();
	std::ofstream fout;
	fout.open(filename);
	fout << str << std::endl;
	fout.close();
	delete[] pointsb_out;
	delete[] points_type_out;
	delete[] points_property_out;
}
//****************************************************************
// Lattice construction and import routines
void init_lat(double lattice_para, char lattice_type[], double sizebox[3][2], double xdir[3], double zdir[3])
// Read user-specified lattice parameters
// Initialize the global lattice type
{
	//lattice parameter
	if (lattice_para > 0)
	{
		lat.lattice_para = lattice_para;
	}
	else
	{
		cout << "Error: lattice parameter should be larger than 0" << endl;
		exit(1);
	}
	//lattice type
	lat.lattice_type = lattice_type;
	//simulation box limits
	for (int i = 0; i < 3; i++)
	{
		if (sizebox[i][0] < sizebox[i][1])
		{
			for (int j = 0; j < 2; j++)
			{
				lat.sizebox[i][j] = sizebox[i][j];
			}
		}
		else
		{
			cout << "Error: upper limits of the box should be larger than the lower limits" << endl;
			exit(1);
		}
	}
	//crystal orientation
	orthogonalization(xdir, zdir);
	double ydir[3];
	cross(zdir, xdir, ydir);
	for (int i = 0; i < 3; i++)
	{
		lat.Mb2o[i][0] = xdir[i] / norm(xdir);
		lat.Mb2o[i][1] = ydir[i] / norm(ydir);
		lat.Mb2o[i][2] = zdir[i] / norm(zdir);
	}
}

void import_lattice(char filename[], double lattice_para, char lattice_type[], double sizebox[3][2], double xdir[3], double zdir[3]) 
// Import an existing structure
// Import a lattice from a LAMMPS-style data file
{
	init_lat(lattice_para, lattice_type, sizebox, xdir, zdir);
	const char type1[] = "sc", type2[] = "bcc", type3[] = "fcc";
	double (*bbcc)[3], (*bre0)[3], length_b;
	const double bbcc_sc[3][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
	const double bbcc_bcc[7][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.5, -0.5, 0.5 }, { 0.5, 0.5, -0.5 } };
	const double bbcc_fcc[6][3] = { {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5} , {0.0, 0.5, 0.5}, {-0.5, 0.5, 0.0}, {-0.5, 0.0, 0.5} , {0.0, -0.5, 0.5} };
	// choose the nearest-neighbor vectors and basis atoms
	if (!strcmp(lat.lattice_type, type1))
	{
		lat.bbcc_num = 3;
		bbcc = new double[lat.bbcc_num][3];
		memcpy(bbcc, bbcc_sc, 3 * lat.bbcc_num * sizeof(double));
	}
	else if (!strcmp(lat.lattice_type, type2))
	{
		lat.bbcc_num = 7;
		bbcc = new double[lat.bbcc_num][3];
		memcpy(bbcc, bbcc_bcc, 3 * lat.bbcc_num * sizeof(double));
	}
	else if (!strcmp(lat.lattice_type, type3))
	{
		lat.bbcc_num = 6;
		bbcc = new double[lat.bbcc_num][3];
		memcpy(bbcc, bbcc_fcc, 3 * lat.bbcc_num * sizeof(double));
	}
	else
	{
		cout << "Input lattice type is not included" << endl;
		exit(1);
	}
	// compute the facet normals of the Wigner-Seitz cell
	bre0 = new double[lat.bbcc_num][3];
	for (int i = 0; i < lat.bbcc_num; i++)
	{
		length_b = dot(bbcc[i], bbcc[i]);
		for (int j = 0; j < 3; j++)
		{
			bre0[i][j] = bbcc[i][j] / (length_b * lat.lattice_para);
		}
	}
	delete[] bbcc;
	//transform vectors into the box coordinate system
	lat.bre = new double[lat.bbcc_num][3];
	memset(lat.bre, 0, lat.bbcc_num * 3 * sizeof(double));
	for (int i = 0; i < lat.bbcc_num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				lat.bre[i][j] += bre0[i][k] * lat.Mb2o[k][j];
			}
		}
	}
	delete[] bre0;
	//import atomnum and pointsb from the input file
	fstream in;
	int flag_Atoms = 0, Atoms_row = -1, flag_atoms = 0, flag_x = 0, flag_y = 0, flag_z = 0, count_row = 0;
	string buf_s, temp_s;
	in.open(filename);
	string strBuff1;
	getline(in, strBuff1);//skip the first line
	while (!in.eof())
	{
		count_row++;
		string strBuff;
		getline(in, strBuff); //read line by line
		if (flag_Atoms == 0) //check whether atom-coordinate import has started
		{
			int position_space = 0;
			while (position_space < int(strBuff.length())) //scan all characters in the current line
			{
				if (isspace(strBuff[position_space]) == 0)
				{
					if (strBuff.at(position_space) != '#') //if the line is not a comment, search for keywords
					{
						if (strBuff.find("atoms") != strBuff.npos)
						{
							stringstream ss(strBuff);
							ss >> lat.atomnum >> temp_s;
							ss.str("");
							if (temp_s == "atoms" && lat.atomnum > 0)
							{
								lat.pointsb = new double[lat.atomnum][3];
								flag_atoms = 1;
								break;
							}
							else
							{
								break;
							}
						}
						if (strBuff.find("xlo xhi") != strBuff.npos)
						{
							stringstream ss(strBuff);
							if (lat.sizebox[0][1] - lat.sizebox[0][0] < eps) //if the box size was not fully specified, read it from the file
							{
								ss >> lat.sizebox[0][0] >> lat.sizebox[0][1] >> temp_s;
								ss.str("");
								cout << temp_s << endl;
								if (lat.sizebox[0][1] - lat.sizebox[0][0] > eps && temp_s == "xlo xhi")
								{
									flag_x = 1;
								}
							}
							else
							{
								flag_x = 1;
							}
							break;
						}
						if (strBuff.find("ylo yhi") != strBuff.npos)
						{
							stringstream ss(strBuff);
							if (lat.sizebox[1][1] - lat.sizebox[1][0] < eps)
							{
								ss >> lat.sizebox[1][0] >> lat.sizebox[1][1] >> temp_s;
								ss.str("");
								if (lat.sizebox[1][1] - lat.sizebox[1][0] > eps && temp_s == "ylo yhi")
								{
									flag_y = 1;
								}
							}
							else
							{
								flag_y = 1;
							}
							break;
						}
						if (strBuff.find("zlo zhi") != strBuff.npos)
						{
							stringstream ss(strBuff);
							if (lat.sizebox[2][1] - lat.sizebox[2][0] < eps)
							{
								ss >> lat.sizebox[2][0] >> lat.sizebox[2][1] >> temp_s;
								ss.str("");
								if (lat.sizebox[2][1] - lat.sizebox[2][0] > eps && temp_s == "zlo zhi")
								{
									flag_z = 1;
								}
							}
							else
							{
								flag_z = 1;
							}
							break;
						}
						if (strBuff.find("Atoms") != strBuff.npos) //start reading atom coordinates
						{
							flag_Atoms = 1;
							break;
						}
						break;
					}
					else
					{
						break; // lines starting with # are treated as comments
					}
				}
				else
				{
					position_space++;
					continue; //skip spaces and tabs
				}
			}
		}
		else
		{
			if (Atoms_row == -1) //skip the first row after the Atoms header
			{
				Atoms_row = 0;
			}
			else
			{
				stringstream ss(strBuff); //parse atom-coordinate entries
				int temp_id, temp_type;
				ss >> temp_id >> temp_type >> lat.pointsb[Atoms_row][0] >> lat.pointsb[Atoms_row][1] >> lat.pointsb[Atoms_row][2];
				ss.str("");
				Atoms_row++;
				if (Atoms_row == lat.atomnum)
				{
					break;
				}
			}
		}
	}
	in.close();
	if (flag_Atoms == 0)
	{
		cout << "Error: import file in section Atoms" << endl;
		exit(1);
	}
	if (flag_atoms == 0)
	{
		cout << "Error: import file in keyword atoms" << endl;
		exit(1);
	}
	if (flag_x == 0)
	{
		cout << "Error: import file in keyword xlo xhi" << endl;
		exit(1);
	}
	if (flag_y == 0)
	{
		cout << "Error: import file in keyword ylo yhi" << endl;
		exit(1);
	}
	if (flag_z == 0)
	{
		cout << "Error: import file in keyword zlo zhi" << endl;
		exit(1);
	}
	if (lat.atomnum != Atoms_row)
	{
		cout << "Error: number of atoms imported is not consistent with expected" << endl;
		exit(1);
	}
}

void create_lattice(double lattice_para, char lattice_type[], double sizebox[3][2], double xdir[3], double zdir[3]) 
// Generate a lattice from ideal crystallographic data
// Create the initial crystal lattice
{
	init_lat(lattice_para, lattice_type, sizebox, xdir, zdir);
	const char type1[] = "sc", type2[] = "bcc", type3[] = "fcc";
	double(*basis)[3], (*bbcc)[3], (*bre0)[3], length_b, (*pointsb0)[3], (*pointsb1)[3], big_ratio = 5.0;
	int basis_num, count, latticenum, atomcount = 0, (*outbox), minlattice[3], maxlattice[3];
	const double basis_sc[1][3] = { {0.0, 0.0, 0.0} };
	const double bbcc_sc[3][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
	const double basis_bcc[2][3] = { {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5} };
	const double bbcc_bcc[7][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.5, -0.5, 0.5 }, { 0.5, 0.5, -0.5 } };
	const double basis_fcc[4][3] = { {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5} , {0.0, 0.5, 0.5} };
	const double bbcc_fcc[6][3] = { {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5} , {0.0, 0.5, 0.5}, {-0.5, 0.5, 0.0}, {-0.5, 0.0, 0.5} , {0.0, -0.5, 0.5} };
	double boxedge[8][3], boxedgelattice[3][8];
	// nearest-neighbor vectors and basis atoms of the lattice
	if (!strcmp(lat.lattice_type, type1))
	{
		basis_num = 1;
		lat.bbcc_num = 3;
		basis = new double[basis_num][3];
		bbcc = new double[lat.bbcc_num][3];
		memcpy(basis, basis_sc, 3 * basis_num * sizeof(double));
		memcpy(bbcc, bbcc_sc, 3 * lat.bbcc_num * sizeof(double));
	}
	else if (!strcmp(lat.lattice_type, type2))
	{
		basis_num = 2;
		lat.bbcc_num = 7;
		basis = new double[basis_num][3];
		bbcc = new double[lat.bbcc_num][3];
		memcpy(basis, basis_bcc, 3 * basis_num * sizeof(double));
		memcpy(bbcc, bbcc_bcc, 3 * lat.bbcc_num * sizeof(double));
	}
	else if (!strcmp(lat.lattice_type, type3))
	{
		basis_num = 4;
		lat.bbcc_num = 6;
		basis = new double[basis_num][3];
		bbcc = new double[lat.bbcc_num][3];
		memcpy(basis, basis_fcc, 3 * basis_num * sizeof(double));
		memcpy(bbcc, bbcc_fcc, 3 * lat.bbcc_num * sizeof(double));
	}
	else
	{
		cout << "Input lattice type is not included" << endl;
		exit(1);
	}
	// facet normals of the Wigner-Seitz cell
	bre0 = new double[lat.bbcc_num][3];
	for (int i = 0; i < lat.bbcc_num; i++)
	{
		length_b = dot(bbcc[i], bbcc[i]);
		for (int j = 0; j < 3; j++)
		{
			bre0[i][j] = bbcc[i][j] / (length_b * lat.lattice_para);
		}
	}
	delete[] bbcc;
	// transform into the box coordinate system
	lat.bre = new double[lat.bbcc_num][3];
	memset(lat.bre, 0, lat.bbcc_num * 3 * sizeof(double));
	for (int i = 0; i < lat.bbcc_num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				lat.bre[i][j] += bre0[i][k] * lat.Mb2o[k][j];
			}
		}
	}
	delete[] bre0;
	// eight corners of the box
	count = 0;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				boxedge[count][0] = lat.sizebox[0][i];
				boxedge[count][1] = lat.sizebox[1][j];
				boxedge[count][2] = lat.sizebox[2][k];
				if (i == 0)
				{
					boxedge[count][0] -= big_ratio * lat.lattice_para;
				}
				else
				{
					boxedge[count][0] += big_ratio * lat.lattice_para;
				}
				if (j == 0)
				{
					boxedge[count][1] -= big_ratio * lat.lattice_para;
				}
				else
				{
					boxedge[count][1] += big_ratio * lat.lattice_para;
				}
				if (k == 0)
				{
					boxedge[count][2] -= big_ratio * lat.lattice_para;
				}
				else
				{
					boxedge[count][2] += big_ratio * lat.lattice_para;
				}
				count++;
			}
		}
	}
	// box corners in crystal coordinates
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			boxedgelattice[j][i] = dot(boxedge[i], lat.Mb2o[j]) / lat.lattice_para;
		}
	}
	// crystal-coordinate range needed to cover the box
	for (int i = 0; i < 3; i++)
	{
		minlattice[i] = (int)floor(*min_element(boxedgelattice[i], boxedgelattice[i] + 8));
		maxlattice[i] = (int)ceil(*max_element(boxedgelattice[i], boxedgelattice[i] + 8));
	}
	// generate the crystal
	latticenum = (maxlattice[0] - minlattice[0] + 1) * (maxlattice[1] - minlattice[1] + 1) * (maxlattice[2] - minlattice[2] + 1);
	lat.atomnum = latticenum * basis_num;
	pointsb0 = new double[lat.atomnum][3];
	pointsb1 = new double[lat.atomnum][3];
	memset(pointsb1, 0, lat.atomnum * 3 * sizeof(double));
	count = 0;
	for (int i = 0; i < basis_num; i++)
	{
		for (int j = minlattice[0]; j <= maxlattice[0]; j++)
		{
			for (int k = minlattice[1]; k <= maxlattice[1]; k++)
			{
				for (int l = minlattice[2]; l <= maxlattice[2]; l++)
				{
					pointsb0[count][0] = (j + basis[i][0]) * lat.lattice_para;
					pointsb0[count][1] = (k + basis[i][1]) * lat.lattice_para;
					pointsb0[count][2] = (l + basis[i][2]) * lat.lattice_para;
					for (int m = 0; m < 3; m++)
					{
						for (int n = 0; n < 3; n++)
						{
							pointsb1[count][m] += pointsb0[count][n] * lat.Mb2o[n][m];
						}
					}
					count++;
				}
			}
		}
	}
	delete[] basis;
	delete[] pointsb0;
	// delete atoms outside the expanded box
	outbox = new int[lat.atomnum];
	memset(outbox, 0, lat.atomnum * sizeof(int));
	atomcount = lat.atomnum;
	for (int i = 0; i < lat.atomnum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if ((pointsb1[i][j] < lat.sizebox[j][0] - big_ratio * lat.lattice_para - eps) || (pointsb1[i][j] > lat.sizebox[j][1] + big_ratio * lat.lattice_para + eps))
			{
				outbox[i] = 1;
				atomcount--;
				break;
			}
		}
	}
	lat.pointsb = new double[atomcount][3];
	memset(lat.pointsb, 0, atomcount * 3 * sizeof(double));
	count = 0;
	for (int i = 0; i < lat.atomnum; i++)
	{
		if (outbox[i] != 1)
		{
			for (int j = 0; j < 3; j++)
			{
				lat.pointsb[count][j] = pointsb1[i][j];
			}
			count++;
		}
	}
	lat.atomnum = count;
	delete[] pointsb1;
	delete[] outbox;
}
//****************************************************************
// Dislocation data structures and operations
struct dislocation
// Dislocation information
{
	char* type;
	int numP;						// number of dislocation-line segments
	double b[3];					// Burgers vector of the dislocation
	double(*Points)[3];		// vertices of the dislocation line
	double trace_of_cut[3];
};
	//****************************
	// Define commonly used dislocation geometries
dislocation polygon(double center[3], double axisa[3], double a, double axisc[3], double c, double b[3], int n)
{
	dislocation disl;
	double axisa_real[3], axisc_real[3], theta;
	orthogonalization(axisa, axisc);
	theta = 2 * pi / n;
	disl.numP = n;
	disl.Points = new double[disl.numP][3];
	for (int i = 0; i < 3; i++)
	{
		disl.b[i] = b[i];
		axisa_real[i] = a * axisa[i] / norm(axisa);
		axisc_real[i] = c * axisc[i] / norm(axisc);
	}
	for (int i = 0; i < disl.numP; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			disl.Points[i][j] = center[j] + axisa_real[j] * cos(i * theta) + axisc_real[j] * sin(i * theta);
		}
	}
	return disl;
}

dislocation create_ellipse(double center[3], double axisa[3], double a, double axisc[3], double c, double b[3])
// Elliptical dislocation loop
{
	double theta0;
	int n;
	dislocation disl;
	if (a < eps || c < eps)
	{
		cout << "the input value should be larger than zero" << endl;
		exit(1);
	}
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	orthogonalization(axisa, axisc);
	theta0 = 2.0 * lat.lattice_para / max(a, c);
	n = (int)ceil(2 * pi / theta0);
	disl = polygon(center, axisa, a, axisc, c, b, n);
	disl.type = "ellipse";
	memset(disl.trace_of_cut, 0, 3 * sizeof(double));
	return disl;
}

dislocation create_circle( double center[3], double normal[3], double radius, double b[3])
// Circular dislocation loop
{
	double axisa[3], axisc[3];
	dislocation disl;
	if (radius < eps)
	{
		cout << "the input value should be larger than zero" << endl;
		exit(1);
	}
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	if (norm(normal) < eps)
	{
		cout << "the input vector should not be zero vector" << endl;
		exit(1);
	}
	if (normal[0] > eps)
	{
		axisa[0] = -normal[1];
		axisa[1] = normal[0];
		axisa[2] = 0.0;
	}
	else if (normal[1] > eps)
	{
		axisa[1] = -normal[2];
		axisa[2] = normal[1];
		axisa[0] = 0.0;
	}
	else if (normal[2] > eps)
	{
		axisa[2] = -normal[0];
		axisa[0] = normal[2];
		axisa[1] = 0.0;
	}
	else
	{
		cout << "the normal of the circle should not be zero vector" << endl;
		exit(1);
	}
	cross(normal, axisa, axisc);
	disl = create_ellipse(center, axisa, radius, axisc, radius, b);
	disl.type = "circle";
	memset(disl.trace_of_cut, 0, 3 * sizeof(double));
	return disl;
}

dislocation create_polygon(double center[3], double normal[3], double edge_dir[3], double edge_len, int n, double b[3])
// Regular polygonal dislocation loop
{
	dislocation disl;
	double radius, theta, rotation_matrix[3][3], axisa[3], length_normal;
	if (edge_len < eps)
	{
		cout << "the input value should be larger than zero" << endl;
		exit(1);
	}
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	if (n < 3)
	{
		cout << "the number of edges should be at least three" << endl;
		exit(1);
	}
	orthogonalization(normal, edge_dir);
	theta = pi / n;
	radius = edge_len / (2 * sin(theta));
	length_normal = norm(normal);
	for (int i = 0; i < 3; i++)
	{
		normal[i] /= length_normal;
	}
	rotation_matrix[0][0] = sin(theta) + normal[0] * normal[0] * (1 - sin(theta));
	rotation_matrix[0][1] = normal[0] * normal[1] * (1 - sin(theta)) - normal[2] * cos(theta);
	rotation_matrix[0][2] = normal[0] * normal[2] * (1 - sin(theta)) + normal[1] * cos(theta);
	rotation_matrix[1][0] = normal[1] * normal[0] * (1 - sin(theta)) + normal[2] * cos(theta);
	rotation_matrix[1][1] = sin(theta) + normal[1] * normal[1] * (1 - sin(theta));
	rotation_matrix[1][2] = normal[1] * normal[2] * (1 - sin(theta)) - normal[0] * cos(theta);
	rotation_matrix[2][0] = normal[2] * normal[0] * (1 - sin(theta)) - normal[1] * cos(theta);
	rotation_matrix[2][1] = normal[2] * normal[1] * (1 - sin(theta)) + normal[0] * cos(theta);
	rotation_matrix[2][2] = sin(theta) + normal[2] * normal[2] * (1 - sin(theta));
	for (int i = 0; i < 3; i++)
	{
		axisa[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			axisa[i] += rotation_matrix[i][j] * edge_dir[j];
		}
	}
	double axisc[3];
	cross(normal, axisa, axisc);
	disl = polygon(center, axisa, radius, axisc, radius, b, n);
	disl.type = "regular polygon";
	memset(disl.trace_of_cut, 0, 3 * sizeof(double));
	return disl;
}

dislocation create_rectangle(double center[3], double axisa[3], double a, double axisc[3], double c, double b[3])
// Rectangular dislocation loop
{
	dislocation disl;
	double axisa_real[3], axisc_real[3];
	disl.type = "rectangle";
	memset(disl.trace_of_cut, 0, 3 * sizeof(double));
	if (a < eps || c < eps)
	{
		cout << "the input value should be larger than zero" << endl;
		exit(1);
	}
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	orthogonalization(axisa, axisc);
	disl.numP = 4;
	disl.Points = new double[disl.numP][3];
	for (int i = 0; i < 3; i++)
	{
		disl.b[i] = b[i];
		axisa_real[i] = a * axisa[i] / norm(axisa);
		axisc_real[i] = c * axisc[i] / norm(axisc);
		disl.Points[0][i] = center[i] + axisa_real[i] + axisc_real[i];
		disl.Points[1][i] = center[i] - axisa_real[i] + axisc_real[i];
		disl.Points[2][i] = center[i] - axisa_real[i] - axisc_real[i];
		disl.Points[3][i] = center[i] + axisa_real[i] - axisc_real[i];
	}
	return disl;
}

dislocation create_helix(double center[3], double normal[3], double axisa[3], double radius, double pitch, double turns, double b[3])
// Helical dislocation line
{
	double theta0, theta, axisc[3], xdir[3], ydir[3], zdir[3];
	int n;
	dislocation disl;
	disl.type = "helix";
	if (radius < eps || pitch < eps || turns < eps)
	{
		cout << "the input value should be larger than zero" << endl;
		exit(1);
	}
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	orthogonalization(normal, axisa);
	theta0 = 2.0 * lat.lattice_para / sqrt(radius * radius + pitch * pitch / (4 * pi * pi));
	n = (int) ceil(turns * ceil(0.5 * 2 * pi / theta0));
	disl.numP = (2 * n + 1) + 6;
	theta = turns * 2 * pi / (2.0 * n + 0.0);
	cross(normal, axisa, axisc);
	for (int i = 0; i < 3; i++)
	{
		disl.b[i] = b[i];
		xdir[i] = axisa[i] / norm(axisa);
		ydir[i] = axisc[i] / norm(axisc);
		zdir[i] = normal[i] / norm(normal);
	}
	memcpy(disl.trace_of_cut, ydir, 3 * sizeof(double));
	disl.Points = new double[disl.numP][3];
	for (int i = -n; i <= n ; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			disl.Points[i + n][j] = center[j] + radius * cos(i * theta) * xdir[j] + radius * sin(i * theta) * ydir[j] + pitch / (2 * pi) * (i * theta) * zdir[j];
		}
	}
	for (int i = 0; i < 3; i++)
	{
		disl.Points[2 * n + 1][i] = center[i] + pitch / (2 * pi) * (n * theta) * zdir[i];
		disl.Points[2 * n + 2][i] = center[i] + pitch * zdir[i] / (2 * pi) * (n * theta) + 1000000.0 * zdir[i];
		disl.Points[2 * n + 3][i] = center[i] + pitch * zdir[i] / (2 * pi) * (n * theta) + 1000000.0 * zdir[i] + 1000000.0 * ydir[i];
		disl.Points[2 * n + 4][i] = center[i] - pitch * zdir[i] / (2 * pi) * (n * theta) - 1000000.0 * zdir[i] + 1000000.0 * ydir[i];
		disl.Points[2 * n + 5][i] = center[i] - pitch * zdir[i] / (2 * pi) * (n * theta) - 1000000.0 * zdir[i];
		disl.Points[2 * n + 6][i] = center[i] - pitch / (2 * pi) * (n * theta) * zdir[i];
	}
	return disl;
}

dislocation create_pointlist(double points[][3], int n, double b[3], double traceofcut[3])
// User-defined polyline dislocation
{
	dislocation disl;
	disl.type = "pointlist";
	if (norm(b) < eps)
	{
		cout << "the input Burgers vector should not be zero vector" << endl;
		exit(1);
	}
	memcpy(disl.b, b, 3 * sizeof(double));
	if (norm(traceofcut) < eps)
	{
		memset(disl.trace_of_cut, 0, 3 * sizeof(double));
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				disl.trace_of_cut[i] += b[j] * lat.Mb2o[j][i];
			}
			disl.trace_of_cut[i] /= norm(b);
		}
	}
	else
	{
		memcpy(disl.trace_of_cut, traceofcut, 3 * sizeof(double));
		double length_tr = norm(traceofcut);
		for (int i = 0; i < 3; i++)
		{
			disl.trace_of_cut[i] /= length_tr;
		}
	}
	if (n > 2)
	{
		disl.numP = n;
	}
	else
	{
		cout << "Error: the number of points defining the dislocation line should be at least three\n" << endl;
		exit(1);
	}
	disl.Points = new double[disl.numP][3];
	memcpy(disl.Points, points, disl.numP * 3 * sizeof(double));
	return disl;
}
	//****************************
void construct_dislocations(dislocation disl)
// Insert a dislocation into the current lattice
{
	int count = 0, Maxneigh, idx;
	double  pstart[3], pend[3], pcenter[3], b_unit[3] = { 0 }, additional_vector[3] = { 0 }, add_v[3] = { 0 }, pointsb_a[3], pointsb_b[3], pointsb_c[3], * solidangle_total, (*pointsb_add)[3], box_l[3];
	double sizebox_big[3][2] = { {1.0e7, -1.0e7}, {1.0e7, -1.0e7}, {1.0e7, -1.0e7} };
	int(*bin_of_poi_xyz)[3], * bin_of_poi_index, ** poi_in_bin_list, * poi_in_bin_list_count;
	int Nbin[3], Nbin_total, binxyz[3][3], bin_index, idn = -1, * point_type0, * point_type, atomnum_add, (*delete_flag);
	double  vector_temp[3], maxdp, lneigh = 1.1 * lat.lattice_para, lneigh2 = lneigh * lneigh, lneigh_repr = 1.0 / lneigh;
	//direction of the Burgers vector in box coordinates
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b_unit[i] += disl.b[j] * lat.Mb2o[j][i];
		}
		b_unit[i] /= norm(disl.b);
	}
	//special treatment for FCC partial dislocations
	if (!strcmp(lat.lattice_type, "fcc"))
	{
		double shockley_b[24][3] = { {1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}, {1.0 / 6.0, 1.0 / 6.0, -2.0 / 6.0}, \
														{-1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}, { 1.0 / 6.0, -1.0 / 6.0, 2.0 / 6.0 }, \
														{1.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0}, { 1.0 / 6.0, -2.0 / 6.0, 1.0 / 6.0 }, \
														{-1.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0}, { 1.0 / 6.0, 2.0 / 6.0, -1.0 / 6.0 }, \
														{2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}, { -2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 }, \
														{2.0 / 6.0, -1.0 / 6.0, 1.0 / 6.0}, { 2.0 / 6.0, 1.0 / 6.0, -1.0 / 6.0 } };
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				shockley_b[i + 12][j] = -shockley_b[i][j];
			}
		}
		double frank_b[8][3] = { {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, {-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, \
												{1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0}, { 1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0 }, \
												{-1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0}, { 1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0 }, \
												{-1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0}, { -1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0 } };
		double addtional_v_frank_b[8][3] = { {-1.0 / 6.0, -1.0 / 6.0, 2.0 / 6.0}, {1.0 / 6.0, -1.0 / 6.0, 2.0 / 6.0}, \
																	{-1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}, {-1.0 / 6.0, -1.0 / 6.0, -2.0 / 6.0 }, \
																	{1.0 / 6.0, 1.0 / 6.0, -2.0 / 6.0}, { -1.0 / 6.0, 1.0 / 6.0, -2.0 / 6.0 }, \
																	{1.0 / 6.0, -1.0 / 6.0, -2.0 / 6.0}, { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0 } };
		//Shockley partial dislocation
		double deviation_vector[3];
		for (int i = 0; i < 24; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				deviation_vector[j] = disl.b[j] - shockley_b[i][j];
			}
			if (norm(deviation_vector) < 0.1)
			{
				cout << "Warning: The partial dislocation line should be parallel to the stacking fault plane and should not contain exactly an atomic plane" << endl;
				fflush(stdout);
				memcpy(disl.b, shockley_b[i], 3 * sizeof(double));
				break;
			}
		}
		//Frank partial dislocation
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				deviation_vector[j] = disl.b[j] - frank_b[i][j];
			}
			if (norm(deviation_vector) < 0.1)
			{
				cout << "Warning: The partial dislocation line should be parallel to the stacking fault plane and should not contain exactly an atomic plane\n" << endl;
				memcpy(disl.b, frank_b[i], 3 * sizeof(double));
				memcpy(additional_vector, addtional_v_frank_b[i], 3 * sizeof(double));
				break;
			}
		}
	}
	//compute the solid angle
	solidangle_total = new double[lat.atomnum];
	memset(solidangle_total, 0, lat.atomnum * sizeof(double));
	for (int i = 0; i < 3; i++)
	{
		pcenter[i] = 0;
		for (int j = 0; j < disl.numP; j++)
		{
			pcenter[i] += disl.Points[j][i];
		}
		pcenter[i] /= disl.numP;
	}
	for (int i = 0; i < disl.numP; i++)
	{
		// compute the solid angle
		memcpy(pstart, disl.Points[i], 3 * sizeof(double));
		if (i == disl.numP - 1)
		{
			memcpy(pend, disl.Points[0], 3 * sizeof(double));
		}
		else
		{
			memcpy(pend, disl.Points[i + 1], 3 * sizeof(double));
		}
		for (int j = 0; j < lat.atomnum; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				pointsb_a[k] = pstart[k] - lat.pointsb[j][k];
				pointsb_b[k] = pend[k] - lat.pointsb[j][k];
				pointsb_c[k] = pcenter[k] - lat.pointsb[j][k];
			}
			if ((!strcmp(disl.type, "helix")) || (!strcmp(disl.type, "pointlist")))
			{
				solidangle_total[j] += solidangle_tri(pointsb_a, pointsb_b, disl.trace_of_cut);
			}
			else
			{
				solidangle_total[j] += solidangle_tri(pointsb_a, pointsb_b, pointsb_c);
			}
		}
	}
	for (int i = 0; i < lat.atomnum; i++)
	{
		solidangle_total[i] *= 2.0;
	}

	//bounding box and dimensions of the current atomic configuration
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < lat.atomnum; j++)
		{
			if (sizebox_big[i][0] > lat.pointsb[j][i])
				sizebox_big[i][0] = lat.pointsb[j][i] - 0.1;
			if (sizebox_big[i][1] < lat.pointsb[j][i])
				sizebox_big[i][1] = lat.pointsb[j][i] + 0.1;
		}
		box_l[i] = sizebox_big[i][1] - sizebox_big[i][0];
	}
	// classify atoms using a binning algorithm
	// build each atom bin index and each bin-to-atom list
	for (int i = 0; i < 3; i++)
	{
		Nbin[i] = (int)ceil(box_l[i] * lneigh_repr);
	}
	Nbin_total = Nbin[0] * Nbin[1] * Nbin[2];
	bin_of_poi_xyz = new int[lat.atomnum][3];
	bin_of_poi_index = new int[lat.atomnum];

	for (int i = 0; i < lat.atomnum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			bin_of_poi_xyz[i][j] = (int)floor((lat.pointsb[i][j] - sizebox_big[j][0]) * lneigh_repr);
		}
		bin_of_poi_index[i] = bin_of_poi_xyz[i][0] + Nbin[0] * bin_of_poi_xyz[i][1] + Nbin[0] * Nbin[1] * bin_of_poi_xyz[i][2];
	}
	poi_in_bin_list_count = new int[Nbin_total];
	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
	for (int i = 0; i < lat.atomnum; i++)
	{
		idx = bin_of_poi_index[i];
		poi_in_bin_list_count[idx]++;
	}
	Maxneigh = *max_element(poi_in_bin_list_count, poi_in_bin_list_count + Nbin_total);
	poi_in_bin_list = new int* [Nbin_total];
	for (int i = 0; i < Nbin_total; i++)
	{
		poi_in_bin_list[i] = new int[Maxneigh];
	}
	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
	for (int i = 0; i < lat.atomnum; i++)
	{
		idx = bin_of_poi_index[i];
		poi_in_bin_list[idx][poi_in_bin_list_count[idx]] = i;
		poi_in_bin_list_count[idx]++;
	}
	delete[] bin_of_poi_index;
	// build neighbor lists and identify atoms on the upper side of the cut surface
	point_type0 = new int[lat.atomnum];
	memset(point_type0, 0, lat.atomnum * sizeof(int));
	count = 0;
	for (int i = 0; i < lat.atomnum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			binxyz[j][0] = bin_of_poi_xyz[i][j] - 1;
			binxyz[j][1] = bin_of_poi_xyz[i][j];
			binxyz[j][2] = bin_of_poi_xyz[i][j] + 1;
		}
		for (int j = 0; j < 3; j++)
		{
			if (point_type0[i]) break;
			for (int k = 0; k < 3; k++)
			{
				if (point_type0[i]) break;
				for (int l = 0; l < 3; l++)
				{
					if (point_type0[i])	break;
					if ((binxyz[0][j] >= 0) && (binxyz[0][j] < Nbin[0]) && (binxyz[1][k] >= 0) && (binxyz[1][k] < Nbin[1]) && (binxyz[2][l] >= 0) && (binxyz[2][l] < Nbin[2]))
					{
						bin_index = binxyz[0][j] + Nbin[0] * binxyz[1][k] + Nbin[0] * Nbin[1] * binxyz[2][l];
						for (int m = 0; m < poi_in_bin_list_count[bin_index]; m++)
						{
							if (Maxneigh >= 2)
							{
								idn = poi_in_bin_list[bin_index][m];
							}
							for (int n = 0; n < 3; n++)
							{
								vector_temp[n] = lat.pointsb[i][n] - lat.pointsb[idn][n];
							}
							if ((vector_temp[0] * vector_temp[0] + vector_temp[1] * vector_temp[1] + vector_temp[2] * vector_temp[2] < lneigh2) && (idn != i))
							{
								if (solidangle_total[i] - solidangle_total[idn] >  2.0 * pi * 0.95)
								{
									point_type0[i] = 1;
									count++;
									break;
								}
							}
						}
					}
				}
			}
		}
	}
	delete[] poi_in_bin_list_count;
	delete[] poi_in_bin_list;
	delete[] bin_of_poi_xyz;
	//output_for_test("solidangle", lat.atomnum, point_type0, lat.pointsb, solidangle_total);
	// displace atoms, add atoms where required, and assign temporary atom categories
	atomnum_add = lat.atomnum + count;
	count = 0;
	pointsb_add = new double[atomnum_add][3];
	point_type = new int[atomnum_add];
	memset(point_type, 0, atomnum_add * sizeof(int));
	double b[3];
	for (int i = 0; i < 3; i++)
	{
		b[i] = 0;
		add_v[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			b[i] += disl.b[j] * lat.Mb2o[j][i] * lat.lattice_para;
			add_v[i] += additional_vector[j] * lat.Mb2o[j][i] * lat.lattice_para;
		}
	}
	for (int i = 0; i < lat.atomnum; i++)
	{
		//apply the dislocation displacement field
		for (int j = 0; j < 3; j++)
		{
			pointsb_add[count][j] = lat.pointsb[i][j] + solidangle_total[i] * b[j] / (-4 * pi);
		}
		point_type[count] = point_type0[i];
		count++;
		//add atoms by translating atoms on the upper side of the cut surface
		if ((point_type0[i] == 1) && (atomnum_add > 2))
		{
			for (int j = 0; j < 3; j++)
			{
				pointsb_add[count][j] = pointsb_add[count - 1][j] + b[j] + add_v[j];
			}
			point_type[count] = 2;
			count++;
		}
	}
	delete[] lat.pointsb;
	delete[] point_type0;
	delete[] solidangle_total;

	// update the bounding box
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < atomnum_add; j++)
		{
			if (sizebox_big[i][0] > pointsb_add[j][i])
				sizebox_big[i][0] = pointsb_add[j][i] - 0.1;
			if (sizebox_big[i][1] < pointsb_add[j][i])
				sizebox_big[i][1] = pointsb_add[j][i] + 0.1;
		}
		box_l[i] = sizebox_big[i][1] - sizebox_big[i][0];
	}
	//assign atoms to bins
	bin_of_poi_xyz = new int[atomnum_add][3];
	bin_of_poi_index = new int[atomnum_add];
	for (int i = 0; i < 3; i++)
	{
		Nbin[i] = (int)ceil(box_l[i] * lneigh_repr);
	}
	Nbin_total = Nbin[0] * Nbin[1] * Nbin[2];
	for (int i = 0; i < atomnum_add; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			bin_of_poi_xyz[i][j] = (int)floor((pointsb_add[i][j] - sizebox_big[j][0]) * lneigh_repr);
		}
		bin_of_poi_index[i] = bin_of_poi_xyz[i][0] + Nbin[0] * bin_of_poi_xyz[i][1] + Nbin[0] * Nbin[1] * bin_of_poi_xyz[i][2];
	}
	poi_in_bin_list_count = new int[Nbin_total];
	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
	for (int i = 0; i < atomnum_add; i++)
	{
		idx = bin_of_poi_index[i];
		poi_in_bin_list_count[idx]++;
	}
	Maxneigh = *max_element(poi_in_bin_list_count, poi_in_bin_list_count + Nbin_total);
	poi_in_bin_list = new int* [Nbin_total];
	for (int i = 0; i < Nbin_total; i++)
	{
		poi_in_bin_list[i] = new int[Maxneigh];
	}
	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
	for (int i = 0; i < atomnum_add; i++)
	{
		idx = bin_of_poi_index[i];
		poi_in_bin_list[idx][poi_in_bin_list_count[idx]] = i;
		poi_in_bin_list_count[idx]++;
	}
	delete[] bin_of_poi_index;
	// build neighbor lists
	delete_flag = new int[atomnum_add];
	memset(delete_flag, 0, atomnum_add * sizeof(int));
	count = 0;
	for (int i = 0; i < atomnum_add; i++)
	{
		if (point_type[i] != 0)
		{
			for (int j = 0; j < 3; j++)
			{
				binxyz[j][0] = bin_of_poi_xyz[i][j] - 1;
				binxyz[j][1] = bin_of_poi_xyz[i][j];
				binxyz[j][2] = bin_of_poi_xyz[i][j] + 1;
			}
			for (int j = 0; j < 3; j++)
			{
				if (delete_flag[i]) break;
				for (int k = 0; k < 3; k++)
				{
					if (delete_flag[i]) break;
					for (int l = 0; l < 3; l++)
					{
						if (delete_flag[i]) break;
						if ((binxyz[0][j] >= 0) && (binxyz[0][j] < Nbin[0]) && (binxyz[1][k] >= 0) && (binxyz[1][k] < Nbin[1]) && (binxyz[2][l] >= 0) && (binxyz[2][l] < Nbin[2]))
						{
							bin_index = binxyz[0][j] + Nbin[0] * binxyz[1][k] + Nbin[0] * Nbin[1] * binxyz[2][l];
							for (int m = 0; m < poi_in_bin_list_count[bin_index]; m++)
							{
								if (delete_flag[i]) break;
								idn = poi_in_bin_list[bin_index][m];
								for (int n = 0; n < 3; n++)
								{
									vector_temp[n] = pointsb_add[i][n] - pointsb_add[idn][n];
								}
								if ((vector_temp[0] * vector_temp[0] + vector_temp[1] * vector_temp[1] + vector_temp[2] * vector_temp[2] < lneigh2) && (point_type[idn] < point_type[i]))
								{
									maxdp = 0.0;
									for (int p = 0; p < lat.bbcc_num; p++)
									{
										if (maxdp < fabs(dot(vector_temp, lat.bre[p])))
											maxdp = fabs(dot(vector_temp, lat.bre[p]));
									}
									if (maxdp < 0.5 + eps)
									{
										delete_flag[i] = 1;
										count++;
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	delete[] bin_of_poi_xyz;
	delete[] poi_in_bin_list_count;
	delete[] poi_in_bin_list;

	//remove overlapping or redundant atoms
	lat.atomnum = atomnum_add - count;
	lat.pointsb = new double[lat.atomnum][3];
	lat.points_type = new int[lat.atomnum];
	count = 0;
	for (int i = 0; i < atomnum_add; i++)
	{
		if (delete_flag[i] == 0)
		{
			for (int j = 0; j < 3; j++)
			{
				lat.pointsb[count][j] = pointsb_add[i][j];
				lat.points_type[count] = 0;
			}
			count++;
		}
	}
	delete[] point_type;
	delete[] delete_flag;
	delete[] pointsb_add;
	delete[] disl.Points;
}
//
//void identify_disl()
//// identify dislocations using the Nye tensor
//{
//	double lneigh, box_l[3], (*vector_temp)[3], vector_temp_unit[3], vector_temp_l, (*ref_v_lattice)[3], (*ref_v_box)[3], * ref_v_box_l, (*ref_v_box_unit)[3], (*v_ref), (*v_real);
//	double sizebox_big[3][2] = { {1.0e7, -1.0e7}, {1.0e7, -1.0e7}, {1.0e7, -1.0e7} };
//	double dA_devided_by_b, vector_tmp[3], (*deformation)[9], (*deformation_real);
//	double etensor[3][3][3], deltaG[3][3][3], nye_tensor[9], (*bre_unit)[3];
//	double strain_tol = 0.02, coplanar_max, sva_tol = 0.1, length_lower_tol = 0.85, length_upper_tol = 1.15, cos25_tol = 0.906307787, cos15_tol = 0.965925826;
//	int ref_num, (*bin_of_poi_xyz)[3], * bin_of_poi_index, ** poi_in_bin_list, * poi_in_bin_list_count, *vector_temp_flag;
//	int Nbin[3], Nbin_total, binxyz[3][3], bin_index, idn, idx, Maxneigh, *dislocation_type;
//	// kronecker tensor
//	double eye[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
//	// 
//	memset(etensor, 0, 3 * 3 * 3 * sizeof(double));
//	etensor[0][1][2] = 1.0;
//	etensor[1][2][0] = 1.0;
//	etensor[2][0][1] = 1.0;
//	etensor[1][0][2] = -1.0;
//	etensor[2][1][0] = -1.0;
//	etensor[0][2][1] = -1.0;
//	double shockley_b[12][3] = { {1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}, {1.0 / 6.0, 1.0 / 6.0, -2.0 / 6.0}, \
//												{-1.0 / 6.0, 1.0 / 6.0, 2.0 / 6.0}, { 1.0 / 6.0, -1.0 / 6.0, 2.0 / 6.0 }, \
//												{1.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0}, { 1.0 / 6.0, -2.0 / 6.0, 1.0 / 6.0 }, \
//												{-1.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0}, { 1.0 / 6.0, 2.0 / 6.0, -1.0 / 6.0 }, \
//												{2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}, { -2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 }, \
//												{2.0 / 6.0, -1.0 / 6.0, 1.0 / 6.0}, { 2.0 / 6.0, 1.0 / 6.0, -1.0 / 6.0 } }, shockley_b_unit[12][3];
//	double frank_b[4][3] = { {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, {-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, \
//											{1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0}, { 1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0 } }, frank_b_unit[4][3];
//	memset(shockley_b_unit, 0, 12 * 3 * sizeof(double));
//	memset(frank_b_unit, 0, 4 * 3 * sizeof(double));
//	for (int i = 0; i < 12; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			for (int k = 0; k < 3; k++)
//			{
//				shockley_b_unit[i][j] += shockley_b[i][k] * lat.Mb2o[k][j] * sqrt(6.0);
//			}
//		}
//	}
//	for (int i = 0; i < 4; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			for (int k = 0; k < 3; k++)
//			{
//				frank_b_unit[i][j] += frank_b[i][k] * lat.Mb2o[k][j] * sqrt(3.0);
//			}
//		}
//	}
//	// position of all neighboring atoms of lattices
//	double ref_sc[18][3] = {	{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},\
//											{-1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -1.0},\
//											{1.0, 1.0, 0.0},  {1.0, -1.0, 0.0},  {-1.0, 1.0, 0.0},  {-1.0, -1.0, 0.0},\
//											{1.0, 0.0, 1.0}, { 1.0, 0.0, -1.0 }, { -1.0, 0.0, 1.0 }, { -1.0, 0.0, -1.0 },\
//											{0.0, 1.0, 1.0}, { 0.0, 1.0, -1.0 }, { 0.0, -1.0, 1.0 }, { 0.0, -1.0, -1.0 } };
//	double ref_bcc[14][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.5, -0.5, 0.5 } , { 0.5, 0.5, -0.5 },\
//											{-1.0, 0.0, 0.0}, { 0.0, -1.0, 0.0 }, { 0.0, 0.0, -1.0 }, { -0.5, -0.5, -0.5 }, { 0.5, -0.5, -0.5 }, { -0.5, 0.5, -0.5 }, { -0.5, -0.5, 0.5 } };
//	double ref_fcc[12][3] = { {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5} , {0.0, 0.5, 0.5}, {-0.5, 0.5, 0.0}, {-0.5, 0.0, 0.5} , {0.0, -0.5, 0.5},\
//											{-0.5, -0.5, 0.0}, { -0.5, 0.0, -0.5 }, { 0.0, -0.5, -0.5 }, { 0.5, -0.5, 0.0 }, { 0.5, 0.0, -0.5 }, { 0.0, 0.5, -0.5 } };
//	// copy the info of the crystal
//	lat_disl.lattice_para = lat.lattice_para;
//	lat_disl.lattice_type = lat.lattice_type;
//	memcpy(lat_disl.sizebox, lat.sizebox, 3 * 2 * sizeof(double));
//	memcpy(lat_disl.Mb2o, lat.Mb2o, 3 * 3 * sizeof(double));
//	lat_disl.bbcc_num = lat.bbcc_num;
//	// unitilize the vectors
//	bre_unit = new double[lat.bbcc_num][3];
//	for (int i = 0; i < lat.bbcc_num; i++)
//	{
//		double bre_l = norm(lat.bre[i]);
//		for (int j = 0; j < 3; j++)
//		{
//			bre_unit[i][j] = lat.bre[i][j] / bre_l;
//		}
//	}
//
//	// which lattice type
//	if (!strcmp(lat.lattice_type, "sc"))
//	{
//		ref_num = 18;
//		lneigh = 1.8 * lat.lattice_para;
//		dA_devided_by_b = lat.lattice_para;
//		ref_v_lattice = new double[ref_num][3];
//		memcpy(ref_v_lattice, ref_sc, ref_num * 3 * sizeof(double));
//		coplanar_max = 8;
//	}
//	else if ((!strcmp(lat.lattice_type, "bcc")))
//	{
//		ref_num = 14;
//		dA_devided_by_b = sqrt(6.0) / 2.0 * lat.lattice_para;
//		lneigh = 1.1 * lat.lattice_para;
//		ref_v_lattice = new double[ref_num][3];
//		memcpy(ref_v_lattice, ref_bcc, ref_num * 3 * sizeof(double));
//		coplanar_max = 6;
//	}
//	else if ((!strcmp(lat.lattice_type, "fcc")))
//	{
//		ref_num = 12;
//		dA_devided_by_b = sqrt(6.0) / 3.0 * lat.lattice_para;
//		lneigh = 1.1 * lat.lattice_para;
//		ref_v_lattice = new double[ref_num][3];
//		memcpy(ref_v_lattice, ref_fcc, ref_num * 3 * sizeof(double));
//		coplanar_max = 6;
//	}
//	else
//	{
//		cout << "Input lattice type is not included" << endl;
//		exit(1);
//	}
//	double lneigh_repr = 1.0 / lneigh;
//
//	// transformed into box axis, find length and direction
//	ref_v_box = new double[ref_num][3];
//	ref_v_box_l = new double[ref_num];
//	ref_v_box_unit = new double[ref_num][3];
//	memset(ref_v_box, 0, ref_num * 3 * sizeof(double));
//	for (int i = 0; i < ref_num; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			for (int k = 0; k < 3; k++)
//			{
//				ref_v_box[i][j] += ref_v_lattice[i][k] * lat.Mb2o[k][j] * lat.lattice_para;
//			}
//		}
//		//length 
//		ref_v_box_l[i] = norm(ref_v_box[i]);
//		//direction
//		for (int j = 0; j < 3; j++)
//		{
//			ref_v_box_unit[i][j] = ref_v_box[i][j] / ref_v_box_l[i];
//		}
//	}
//	delete[] ref_v_lattice;
//	// system size
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < lat.atomnum; j++)
//		{
//			if (sizebox_big[i][0] > lat.pointsb[j][i])
//			{
//				sizebox_big[i][0] = lat.pointsb[j][i] - 0.1;
//			}
//			if (sizebox_big[i][1] < lat.pointsb[j][i])
//			{
//				sizebox_big[i][1] = lat.pointsb[j][i] + 0.1;
//			}
//		}
//		box_l[i] = sizebox_big[i][1] - sizebox_big[i][0];
//	}
//
//	////////////////////////////////////////////////////////////////////////////////////////
//	// num of boxes
//	for (int i = 0; i < 3; i++)
//	{
//		Nbin[i] = (int)ceil(box_l[i] * lneigh_repr);
//	}
//	Nbin_total = Nbin[0] * Nbin[1] * Nbin[2];
//	bin_of_poi_xyz = new int[lat.atomnum][3];
//	bin_of_poi_index = new int[lat.atomnum];
//	// id and coordinate of box
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			bin_of_poi_xyz[i][j] = (int)floor((lat.pointsb[i][j] - sizebox_big[j][0]) * lneigh_repr);
//		}
//		bin_of_poi_index[i] = bin_of_poi_xyz[i][0] + Nbin[0] * bin_of_poi_xyz[i][1] + Nbin[0] * Nbin[1] * bin_of_poi_xyz[i][2];
//	}
//	// max number of atoms in a box
//	poi_in_bin_list_count = new int[Nbin_total];
//	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		idx = bin_of_poi_index[i];
//		poi_in_bin_list_count[idx]++;
//	}
//	Maxneigh = *max_element(poi_in_bin_list_count, poi_in_bin_list_count + Nbin_total);
//	poi_in_bin_list = new int* [Nbin_total];
//	for (int i = 0; i < Nbin_total; i++)
//	{
//		poi_in_bin_list[i] = new int[Maxneigh];
//	}
//	// determine atom num in each box, build neighboring list
//	memset(poi_in_bin_list_count, 0, Nbin_total * sizeof(int));
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		idx = bin_of_poi_index[i];
//		poi_in_bin_list[idx][poi_in_bin_list_count[idx]] = i;
//		poi_in_bin_list_count[idx]++;
//	}
//	delete[] bin_of_poi_index;
//	////////////////////////////////////////////////////////////////////////////////////////
//	// atomic deformation matrix
//	deformation = new double[lat.atomnum][9];
//
//	double* strain;
//	strain = new double[lat.atomnum];
//	memset(strain, 0, lat.atomnum * sizeof(double));
//
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			binxyz[j][0] = bin_of_poi_xyz[i][j] - 1;
//			binxyz[j][1] = bin_of_poi_xyz[i][j];
//			binxyz[j][2] = bin_of_poi_xyz[i][j] + 1;
//		}
//		// 
//		int count_temp = 0;
//		int *neighbor_bin_atoms_index;
//		// 
//		// 
//		neighbor_bin_atoms_index = new int[27 * Maxneigh];
//		memset(neighbor_bin_atoms_index, -1, 27 * Maxneigh * sizeof(int));
//		for (int j = 0; j < 3; j++)
//		{
//			for (int k = 0; k < 3; k++)
//			{
//				for (int l = 0; l < 3; l++)
//				{
//					if ((binxyz[0][j] >= 0) && (binxyz[0][j] < Nbin[0]) && (binxyz[1][k] >= 0) && (binxyz[1][k] < Nbin[1]) && (binxyz[2][l] >= 0) && (binxyz[2][l] < Nbin[2]))
//					{
//						bin_index = binxyz[0][j] + Nbin[0] * binxyz[1][k] + Nbin[0] * Nbin[1] * binxyz[2][l];
//						for (int m = 0; m < poi_in_bin_list_count[bin_index]; m++)
//						{
//							if (Maxneigh >= 2)
//							{
//								neighbor_bin_atoms_index[count_temp] = poi_in_bin_list[bin_index][m];
//							}
//							count_temp++;
//						}
//					}
//				}
//			}
//		}
//		// find correspondence
//		vector_temp = new double[count_temp][3];
//		vector_temp_flag = new int[count_temp];
//		memset(vector_temp_flag, -1, count_temp * sizeof(int));
//		int count = 0;
//		int count_flag = 0;
//		for (int j = 0; j < count_temp; j++)
//		{
//			idn = neighbor_bin_atoms_index[j];
//			for (int n = 0; n < 3; n++)
//			{
//				vector_temp[j][n] = lat.pointsb[idn][n] - lat.pointsb[i][n];
//			}
//			vector_temp_l = norm(vector_temp[j]);
//			for (int n = 0; n < 3; n++)
//			{
//				vector_temp_unit[n] = vector_temp[j][n] / vector_temp_l;
//			}
//			if ((vector_temp_l < lneigh) && (idn != i))
//			{
//				// 
//				for (int p = 0; p < ref_num; p++)
//				{
//					if (vector_temp_l < length_upper_tol * ref_v_box_l[p] && vector_temp_l > length_lower_tol * ref_v_box_l[p] && dot(vector_temp_unit, ref_v_box_unit[p]) > cos25_tol)
//					{
//						// 
//						vector_temp_flag[j] = p;
//						// 
//						count_flag++;
//						break;
//					}	
//				}
//			}
//		}
//		if (count_flag > coplanar_max)
//		{
//			// 
//			v_ref = new double[count_flag * 3];
//			v_real = new double[count_flag * 3];
//			int count = 0;
//			for (int j = 0; j < count_temp; j++)
//			{
//				if (vector_temp_flag[j] != -1)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						v_ref[count * 3 + k] = ref_v_box[vector_temp_flag[j]][k];
//						v_real[count * 3 + k] = vector_temp[j][k];
//					}
//					count++;
//				}
//			}
//			int info;
//			//info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR, count_flag, 3, 3, v_real, 3, v_ref, 3, s, rcond, rank);
//			//info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR, count_flag, 3, 3, v_real, 3, v_ref, 3, s, rcond, rank);
//			info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', count_flag, 3, 3, v_real, 3, v_ref, 3);
//			if (info > 0)
//			{
//				cout << "Calculation of deformation encounters an error" << endl;
//				exit(1);
//			}
//			memcpy(deformation[i], v_ref, 9 * sizeof(double));
//			delete[] v_ref;
//			delete[] v_real;
//		}
//		else
//		//
//		{
//			// 
//			v_ref = new double[count_flag * 3];
//			v_real = new double[count_flag * 3];
//			int count = 0;
//			for (int j = 0; j < count_temp; j++)
//			{
//				if (vector_temp_flag[j] != -1)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						v_ref[count * 3 + k] = ref_v_box[vector_temp_flag[j]][k];
//						v_real[count * 3 + k] = vector_temp[j][k];
//					}
//					count++;
//				}
//			}
//			double rcond = 0.0;
//			int s[3] = { 0 }, rank[1];
//			int info;
//			//info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR, count_flag, 3, 3, v_real, 3, v_ref, 3, s, rcond, rank);
//			info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR, count_flag, 3, 3, v_real, 3, v_ref, 3, s, rcond, rank);
//			//info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', count_flag, 3, 3, v_real, 3, v_ref, 3);
//			if (info > 0)
//			{
//				cout << "Calculation of deformation encounters an error" << endl;
//				exit(1);
//			}
//			if (rank[0] >= 3)
//			{
//				memcpy(deformation[i], v_ref, 9 * sizeof(double));
//			}
//			else
//			{
//				memcpy(deformation[i], eye, 9 * sizeof(double));
//			}
//			delete[] v_ref;
//			delete[] v_real;
//			
//		}
//
//		for (int p = 0; p < 3; p++)
//		{
//			for (int q = 0; q < 3; q++)
//			{
//				if (p == q)
//				{
//					strain[i] += (deformation[i][3 * p + q] - 1.0) * (deformation[i][3 * p + q] - 1.0);
//				}
//				else
//				{
//					strain[i] += deformation[i][3 * p + q] * deformation[i][3 * p + q];
//				}
//			}
//		}
//		strain[i] = sqrt(strain[i] / 9);
//		//
//		delete[] vector_temp;
//		delete[] vector_temp_flag;
//		delete[] neighbor_bin_atoms_index;
//	}
//	delete[] ref_v_box;
//	delete[] ref_v_box_l;
//	delete[] ref_v_box_unit;
//
//	////////////////////////////////////////////////////////////////////////////////////////
//	// compute the gradient and curl of the deformation-gradient matrix (Nye tensor)
//	dislocation_type = new int[lat.atomnum];
//	int dislocation_count = 0;
//	memset(dislocation_type, -1, lat.atomnum * sizeof(int));
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		if (strain[i] < strain_tol)
//		{
//			continue;
//		}
//		else
//		{
//			for (int j = 0; j < 3; j++)
//			{
//				binxyz[j][0] = bin_of_poi_xyz[i][j] - 1;
//				binxyz[j][1] = bin_of_poi_xyz[i][j];
//				binxyz[j][2] = bin_of_poi_xyz[i][j] + 1;
//			}
//			//
//			int count_temp = 0;
//			int* neighbor_bin_atoms_index;
//			neighbor_bin_atoms_index = new int[27 * Maxneigh];
//			memset(neighbor_bin_atoms_index, -1, 27 * Maxneigh * sizeof(int));
//			for (int j = 0; j < 3; j++)
//			{
//				for (int k = 0; k < 3; k++)
//				{
//					for (int l = 0; l < 3; l++)
//					{
//						if ((binxyz[0][j] >= 0) && (binxyz[0][j] < Nbin[0]) && (binxyz[1][k] >= 0) && (binxyz[1][k] < Nbin[1]) && (binxyz[2][l] >= 0) && (binxyz[2][l] < Nbin[2]))
//						{
//							bin_index = binxyz[0][j] + Nbin[0] * binxyz[1][k] + Nbin[0] * Nbin[1] * binxyz[2][l];
//							for (int m = 0; m < poi_in_bin_list_count[bin_index]; m++)
//							{
//								idn = poi_in_bin_list[bin_index][m];
//								for (int n = 0; n < 3; n++)
//								{
//									vector_tmp[n] = lat.pointsb[i][n] - lat.pointsb[idn][n];
//								}
//								if ((norm(vector_tmp) < lneigh) && (idn != i))
//								{
//									neighbor_bin_atoms_index[count_temp] = idn;
//									count_temp++;
//								}
//							}
//						}
//					}
//				}
//			}
//			// Nye tensor
//			v_real = new double[count_temp * 3];
//			deformation_real = new double[count_temp * 9];
//			// 
//			for (int j = 0; j < count_temp; j++)
//			{
//				idn = neighbor_bin_atoms_index[j];
//				for (int n = 0; n < 3; n++)
//				{
//					v_real[j * 3 + n] = lat.pointsb[i][n] - lat.pointsb[idn][n];
//				}
//				for (int n = 0; n < 9; n++)
//				{
//					deformation_real[j * 9 + n] = deformation[i][n] - deformation[idn][n];
//				}
//			}
//			if (count_temp > coplanar_max)
//			{
//				// 
//				//double rcond = 0.0;
//				//int s[3] = { 0 }, rank[1];
//				int info;
//				//info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR, count_temp, 3, 9, v_real, 3, deformation_real, 9, s, rcond, rank);
//				//info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR, count_temp, 3, 9, v_real, 3, deformation_real, 9, s, rcond, rank);
//				info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', count_temp, 3, 9, v_real, 3, deformation_real, 9);
//				if (info > 0)
//				{
//					cout << "Calculation of the gradient of deformation matrix encounters an error" << endl;
//					exit(1);
//				}
//				for (int j = 0; j < 3; j++)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						for (int l = 0; l < 3; l++)
//						{
//							deltaG[j][k][l] = deformation_real[9 * j + 3 * k + l];
//						}
//					}
//				}
//				memset(nye_tensor, 0, 9 * sizeof(double));
//				for (int j = 0; j < 3; j++)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						for (int m = 0; m < 3; m++)
//						{
//							for (int n = 0; n < 3; n++)
//							{
//								nye_tensor[3 * j + k] -= etensor[j][m][n] * deltaG[m][n][k];
//							}
//						}
//					}
//				}
//			}
//			else
//			{
//				// 
//				double rcond = 0.0;
//				int s[3] = { 0 }, rank[1];
//				int info;
//				//info = LAPACKE_dgelsd(LAPACK_ROW_MAJOR, count_temp, 3, 9, v_real, 3, deformation_real, 9, s, rcond, rank);
//				info = LAPACKE_dgelsy(LAPACK_ROW_MAJOR, count_temp, 3, 9, v_real, 3, deformation_real, 9, s, rcond, rank);
//				//info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', count_temp, 3, 9, v_real, 3, deformation_real, 9);
//				if (info > 0)
//				{
//					cout << "Calculation of the gradient of deformation matrix encounters an error" << endl;
//					exit(1);
//				}
//				for (int j = 0; j < 3; j++)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						for (int l = 0; l < 3; l++)
//						{
//							deltaG[j][k][l] = deformation_real[9 * j + 3 * k + l];
//						}
//					}
//				}
//				memset(nye_tensor, 0, 9 * sizeof(double));
//				if (rank[0] >= 3)
//				{
//					for (int j = 0; j < 3; j++)
//					{
//						for (int k = 0; k < 3; k++)
//						{
//							for (int m = 0; m < 3; m++)
//							{
//								for (int n = 0; n < 3; n++)
//								{
//									nye_tensor[3 * j + k] -= etensor[j][m][n] * deltaG[m][n][k];
//								}
//							}
//						}
//					}
//				}
//
//			}
//			delete[] neighbor_bin_atoms_index;
//			// 
//			int info;
//			double sva[3];
//			double u[9], vt[9];
//			double maxdp = 0.0;
//			double vt_unit[3];
//			int maxtype = 8;
//			info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', 3, 3, nye_tensor, 3, sva, u, 3, vt, 3);
//			if (info > 0)
//			{
//				cout << "Calculation of singular values encounters an error" << endl;
//				exit(1);
//			}
//			if (sva[0] * dA_devided_by_b > sva_tol)
//			{
//				// 
//				double vt_lr = 1.0 / sqrt(vt[0] * vt[0] + vt[1] * vt[1] + vt[2] * vt[2]);
//				for (int j = 0; j < 3; j++)
//				{
//					vt_unit[j] = vt[j] * vt_lr;
//				}
//				for (int j = 0; j < lat.bbcc_num; j++)
//				{
//					if (maxdp < fabs(dot(vt_unit, bre_unit[j])))
//					{
//						maxdp = fabs(dot(vt_unit, bre_unit[j]));
//						maxtype = j;
//					}
//				}
//				if (!strcmp(lat.lattice_type, "fcc"))
//				{
//					for (int j = 0; j < 12; j++)
//					{
//						if (maxdp < fabs(dot(vt_unit, shockley_b_unit[j])))
//						{
//							maxdp = fabs(dot(vt_unit, shockley_b_unit[j]));
//							maxtype = 6;
//						}
//					}
//					for (int j = 0; j < 4; j++)
//					{
//						if (maxdp < fabs(dot(vt_unit, frank_b_unit[j])))
//						{
//							maxdp = fabs(dot(vt_unit, frank_b_unit[j]));
//							maxtype = 7;
//						}
//					}
//				}
//				if (maxdp > cos15_tol)
//				{
//					dislocation_type[i] = maxtype;
//				}
//				else
//				{
//					dislocation_type[i] = 8;
//				}
//				dislocation_count++;
//			}
//			delete[] v_real;
//			delete[] deformation_real;
//		}
//	}
//	delete[] deformation;
//	delete[] strain;
//	delete[] bin_of_poi_xyz;
//	delete[] poi_in_bin_list_count;
//	delete[] poi_in_bin_list;
//	delete[] bre_unit;
//
//	lat_disl.points_type = new int[dislocation_count];
//	lat_disl.atomnum = dislocation_count;
//	lat_disl.pointsb = new double[dislocation_count][3];
//	int count = 0;
//	for (int i = 0; i < lat.atomnum; i++)
//	{
//		if (dislocation_type[i] != -1)
//		{
//			lat_disl.points_type[count] = dislocation_type[i];
//			memcpy(lat_disl.pointsb[count], lat.pointsb[i], 3 * sizeof(double));
//			count++;
//		}
//	}
//	delete[] dislocation_type;
//}

//****************************************************************
// Output routines
void output(lattice lat, char filename[])
// Write the generated structure
{
	int *outbox, atomcount = lat.atomnum, atomnum_out, *points_type_out, type = 0;
	double (*pointsb_out)[3];
	outbox = new int[lat.atomnum];
	memset(outbox, 0, lat.atomnum * sizeof(int));
	// discard atoms outside the region defined by lat.sizebox
	for (int i = 0; i < lat.atomnum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if ((lat.pointsb[i][j] < lat.sizebox[j][0] - eps) || (lat.pointsb[i][j] > lat.sizebox[j][1] + eps))
			{
				outbox[i] = 1;
				atomcount--;
				break;
			}
		}
	}
	pointsb_out = new double[atomcount][3];
	points_type_out = new int[atomcount];
	memset(pointsb_out, 0, atomcount * 3 * sizeof(double));
	int count = 0;
	for (int i = 0; i < lat.atomnum; i++)
	{
		if (outbox[i] != 1)
		{
			if (type < lat.points_type[i] + 1)
			{
				type = lat.points_type[i] + 1;
			}
			points_type_out[count] = lat.points_type[i];
			for (int j = 0; j < 3; j++)
			{
				pointsb_out[count][j] = lat.pointsb[i][j];
			}
			count++;
		}
	}
	atomnum_out = count;
	delete[] outbox;
	//write the output file
	std::stringstream stream;
	stream << "#Initial atom positions. Genernated by Jinyu Zhang." << "\n";
	stream << "\t\t" <<  atomnum_out << " \tatoms" << "\n";
	stream << "\t\t " << type << "\tatom types" << "\n";
	stream << "\t\t0.0 0.0 0.0\txy xz yz" << "\n";
	stream << "\t\t" << lat.sizebox[0][0] << "\t" << lat.sizebox[0][1] << "\txlo xhi" << "\n";
	stream << "\t\t" << lat.sizebox[1][0] << "\t" << lat.sizebox[1][1] << "\tylo yhi" << "\n";
	stream << "\t\t" << lat.sizebox[2][0] << "\t" << lat.sizebox[2][1] << "\tzlo zhi" << "\n";
	stream << " " << "\n";
	stream << " Atoms" << "\n";
	stream << " " << "\n";
	// points_type_out stores atom types for visualization; pointsb_out stores atom coordinates; sizebox defines the orthogonal box bounds
	for (int i = 0; i < atomnum_out; i++)
	{
		stream << std::setw(9) << (i + 1) << std::setw(3) << points_type_out[i] + 1 << std::setw(15) << pointsb_out[i][0] << std::setw(15) << pointsb_out[i][1] << std::setw(15) << pointsb_out[i][2] << "\n";
	}
	std::string str;
	str = stream.str();

	std::ofstream fout;
	fout.open(filename);
	fout << str << std::endl;
	fout.close();
	delete[] pointsb_out;
	delete[] points_type_out;
}

void visualization(lattice lat)
{
	int* outbox, atomcount = lat.atomnum, atomnum_out, * points_type_out, type = 0;
	double(*pointsb_out)[3];
	outbox = new int[lat.atomnum];
	memset(outbox, 0, lat.atomnum * sizeof(int));
	// discard atoms outside the region defined by lat.sizebox
	for (int i = 0; i < lat.atomnum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if ((lat.pointsb[i][j] < lat.sizebox[j][0] - eps) || (lat.pointsb[i][j] > lat.sizebox[j][1] + eps))
			{
				outbox[i] = 1;
				atomcount--;
				break;
			}
		}
	}
	pointsb_out = new double[atomcount][3];
	points_type_out = new int[atomcount];
	memset(pointsb_out, 0, atomcount * 3 * sizeof(double));
	int count = 0;
	for (int i = 0; i < lat.atomnum; i++)
	{
		if (outbox[i] != 1)
		{
			if (type < lat.points_type[i] + 1)
			{
				type = lat.points_type[i] + 1;
			}
			points_type_out[count] = lat.points_type[i];
			for (int j = 0; j < 3; j++)
			{
				pointsb_out[count][j] = lat.pointsb[i][j];
			}
			count++;
		}
	}
	atomnum_out = count;
	delete[] outbox;

	//additional code for visualization output
	//.........................................
	// atomnum_out stores the number of atoms; points_type_out stores atom types for visualization; pointsb_out stores atom coordinates; sizebox defines the orthogonal box bounds

	delete[] pointsb_out;
	delete[] points_type_out;
}

//****************************************************************
int main()
//
{
	double lattice_para = 2.86;																										// lattice parameter
	char lattice_type[] = "bcc";																										// lattice type
	double sizebox[3][2] = { {-29.67, 29.47}, {-7.02, 7.00}, {-24.28, 24.26} };										// simulation box range
	//double xdir[3] = { 1, -1, 0 }, zdir[3] = { 1, 1, 1 };																	// orientation of crystal 1
	double xdir[3] = {-1, 1, 1 }, zdir[3] = { 1, 0, 1 };																		// orientation of crystal 2

	// import user defined lattice
	//char filename_in[] = "dislocations_atom_initial.txt";
	// import_lattice(filename_in, lattice_para, lattice_type, sizebox, xdir, zdir);

	// create lattice
	create_lattice(lattice_para, lattice_type, sizebox, xdir, zdir);
	dislocation disl;
	//***************************************************
	// point list
	//double b[3] = { 0.5, 0.5, 0.5 };
	//double p[4][3] = { { -70.1, -70.0, 30.0 } , { -25.0, 40.0, -20.0 }, {29.9, 30.0, 0.0 }, {60.0, -44.0, -30.0} };
	////double p[4][3] = { { 100.0, -100.0, 0.1 } , { 100.0, 100.0, 0.1 }, { -100.0, 100.0, 0.1 }, { -100.0, -100.0, 0.1 } };
	//double tr[3] = { 0, 0, 0 };
	//disl = create_pointlist(p, 4, b, tr);
	//construct_dislocations(disl);

	// ellipse
	//double b[3] = { 0.5, 0.5, 0 };
	//double center[3] = { 0.0, 0.0, 0.0 }, axisa[3] = { 0.0, 0.5, -0.5 }, axisc[3] = { 1.0, 0.0, 0.0 }, a = 20.0, c = 40.0;
	//disl = create_ellipse(center, axisa, a, axisc, c, b);
	//construct_dislocations(disl);

	// circle
	//double b[3] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
	////double b[3] = { -1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0 };
	////double b[3] = { 1.0 / 6.0, 1.0 / 6.0, -2.0 / 6.0 };
	//double center[3] = { 0.1, 0.3, -0.5 }, normal[3]= { 0.0, 0.0, 1.0 }, radius = 50.0;
	//disl = create_circle(center, normal, radius, b);
	//construct_dislocations(disl);

	// regular polygon
	// double b[3] = { 0.5, -0.5, 0.0 };
	// double center[3] = { 0.0, 0.0, 0.0 }, normal[3] = { 0.0, 1.0, 1.0 }, axisa[3] = {0.0, 1.0, -1.0}, radius = 40.0;
	// int n = 7;
	// disl = create_polygon(center, normal, axisa, radius, n, b);
	// construct_dislocations(disl);

	//rectangle
	double b[3] = {-0.5, 0.5, 0.5 };
	double center[3] = { 0.4, 0.07, 0.0 }, axisa[3] = { 0.0, 0.0, 1.0 }, axisc[3] = { 0.0, 1.0, 0.0 }, a = 9.0, c = 600.0;
	disl = create_rectangle(center, axisa, a, axisc, c, b);
	construct_dislocations(disl);
	
    // optionally generate a second dislocation here
	//double b1[3] = { -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0 };
	//double center1[3] = { -0.0, -0.0, -50.1 }, normal1[3] = { 1.0, -1.0, 1.0 }, axisa1[3] = { 0.0, 1.0, 1.0 }, radius1 = 40.0;
	//int n1 = 7;
	////disl = create_polygon(center1, normal1, axisa1, radius1, n1, b1);
	//disl = create_circle(center1, normal1, radius1, b1);
	//construct_dislocations(disl);

	//// helical dislocation
	//double b[3] = { 0.5, 0.5, 0.0 };																							// 
	//double center[3] = { 0.0, 0.0, 40.0 }, normal[3] = { 1.0, 1.0, 0.0 }, axisa[3] = { 1.0, -1.0, 0.0 }, radius = 33.1, pitch = 27.0, turns = 6.0;
	//////																																		//
	//disl = create_helix(center, normal, axisa, radius, pitch, turns, b);
	//construct_dislocations(disl);

	//*********************************************
	//
	char filename[] = "dislocations_atom_initial.txt";
	output(lat, filename);
	//*********************************************
	// identify dislocations using the Nye tensor
	char filename_disl[] = "dislocations_atom_initial_disl.txt";
	//identify_disl();
	//
	output(lat_disl, filename_disl);
	//
	delete[] lat.pointsb;
	delete[] lat.points_type;
	delete[] lat.bre;
	delete[] lat_disl.pointsb;
	delete[] lat_disl.points_type;
	return 0;
}