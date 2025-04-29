#include <iostream>
#ifdef __unix
#include <string.h>
#include <unistd.h>
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),  (mode)))==NULL
#else
#include <io.h>
#include <direct.h>
#endif
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <cstdio>
#include <tuple>
#include <random>
#include <limits.h>
#include <iomanip>
#include <unordered_map>
#include <numeric>
#include <stack>
#include <unordered_set>
#include <omp.h>
#include "eigen-3.4.0/Eigen/Eigen"

using namespace std;
using Position = std::vector<double>;       // 表示单个粒子的三维坐标 (x, y, z)
using Configuration = std::vector<Position>; // 表示某一时间点的所有粒子坐标
using Trajectory = std::vector<Configuration>; // 表示所有构型
#include "head.h"
#ifdef gromacs_traj
extern "C" {
#include "xdrfile.h"
#include "xdrfile_xtc.h"
}
#endif
#include "force.h"
#include "rdf.h"
#include "structure_factor.h"
#include "diffusion.h"
#include "ISF.h"
#include "angle.h"
#include "z.h"
#include "eigen.h"
#include "non_affine_D.h"
#include "tetrahedron.h"
#include "bond.h"
#include "matrix.h"
#include "position.h"
#include "read_traj.h"
#include "chi_4.h"
#include "cluster.h"

int main(int argc, char* argv[])
{
	int Tbegin = static_cast<int> (time(NULL));
	project = argv[1];
	if (project == "SiO2") {
		dt = 0;
		first_layer = "FL";
	}else if (project == "AlSm") {
		N_micro = 6 * 21 * 21 * 21;
		micro_form = "6";
		decompose_way = "SVD";
	}
	else if (project == "super_critical") { //"1000" "3375" "8000" "15625" "64000"	
		N_micro = 1000; //矩阵的微观态数量
		micro_form = 1; //矩阵形式
		matrix_name = "lattice_nrho"; //矩阵名称
		decompose_way = "SVD"; //分解方式
		std::cout << "N_micro: " << N_micro << ". micro_form: " << micro_form << ". matrix_name: " << matrix_name
			<< ". decompose_way: " << decompose_way << endl;
	}
	string physics = argv[2];
	if (physics == "MSD" || physics == "sq" || physics == "rdf" || physics == "SISF" || physics == "momentum" || physics == "chi4" || physics == "L" || physics == "ISF" ||
		physics == "local_NA" || physics == "p_u" || physics == "NA" || physics == "posi" || physics == "z" || physics == "cluster" ||
		physics == "matrix" || physics == "eigen") {
		string dir_data = argv[3];
		filename_input = argv[4]; //"traj_comp.xtc"; "traj_angular.atom";
		atom_kind = argv[5];
		total_config = stoi(argv[6]);
		S_k_max = stof(argv[7]);
		distance_cut = stof(argv[8]);
		string i_case = argv[9];
		if (physics == "chi4") {
			distance_cut = distance_cut / 3.0;
		}
		else {
			distance_cut = distance_cut * distance_cut; //一般用距离的平方
		}
		dir_write = "../../data/" + physics + "/";
		if (i_case == "none") {
			dir_read = "../../data/" + dir_data + "/";
			dir_out = dir_write + dir_data + "_";
		}
		else {
			dir_read = "../../data/" + dir_data + "/" + i_case + "/";
			dir_out = dir_write + dir_data + "_" + i_case + "_";
		}
		std::cout << dir_read << ". dir_out: " << dir_out << ". S_k_max: " << S_k_max << ". distance_cut: " << distance_cut << endl;
		if (total_config <= 0) {
			string file_traj = dir_read + filename_input;
			//total_config = CountLines((char*)file_traj.c_str()) / (N + N_infor);
			total_config = count_lines(file_traj) / (10009) - 1;
		}
		matrix_name = "../../data/matrix/" + dir_data + "_" + i_case + "_matrix_" + atom_kind + ".txt";
		eigen_name = "../../data/eigen/" + dir_data + "_" + i_case + "_eigen_" + atom_kind + ".txt";
	}
	else if (physics == "matrix_s" || physics == "eigen_s") {
		N = stoi(argv[2]);
		cout << N << endl;
		string N_str = "N_" + to_string(N);
		string Pressure = argv[3];
		string Temperature = argv[4];
		string Simu_case = argv[5];
		total_config = stoi(argv[6]);
		filename_input = argv[7];
		string out_label = argv[8];
		dir_read = "../../data/Widom_line/" + N_str + "/" + Pressure + "/" + Temperature + "/" + Simu_case + "/";
		string dir_matrix = "../../data/Widom_line/matrix/";
		string dir_eigen = "../../data/Widom_line/eigen/";
		matrix_name = dir_matrix + to_string(N) + "_" + Pressure + "_" + Temperature + "_" + out_label + ".txt";
		eigen_name = dir_eigen + to_string(N) + "_" + Pressure + "_" + Temperature + "_" + out_label + "_eigen.txt";
		std::cout << dir_read << ". filename_input: " << filename_input << ". matrix_name: " << matrix_name << ". eigen_name: " << eigen_name << endl;
		N_micro = N;
		GRID_SIZE = std::round(pow(N, 1.0 / 3.0)); //四舍五入
	}
	int isCreate = createDir(dir_write);
	if (isCreate == 0) {
		cout << "Created path: " << dir_write << endl;
	}
	else {
		printf("error: create path failed!\n");
		exit(-1);
	}
	const string file_traj = dir_read + filename_input;
	std::cout << "ini_config:" << ini_config << ". total_config: " << total_config << ". d_S: " << d_S << ". file_traj: " << file_traj << endl;

#ifdef gen_rdf
	vector<double> rdf_11(N_Bins_r, 0.0); //11
	vector<double> rdf_12(N_Bins_r, 0.0); //12 21
	vector<double> rdf_22(N_Bins_r, 0.0); //22
	vector<double> rdf_All(N_Bins_r, 0.0);
#endif

#if defined(gen_L) || defined(gen_sq) || defined(gen_ISF) || defined(gen_SISF)
	double L_mean[Dim] = { 0.0 };
	get_L_mean(file_traj, L_mean); //得到平均的盒子尺寸
#endif
#if defined(gen_L)
	string outname = dir_out + "L.txt";
	std::ofstream outputFile(outname, ios::out);
	for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
		outputFile << L_mean[i_Dim] << endl;
	}
	outputFile.close();
	return 0;
#endif
#ifdef gen_sq
#ifdef read_wavevectors
	std::vector<Vec3_wave_vector> wavevectors_all = read_wavevector();
#else
	std::vector<Vec3_wave_vector> wavevectors_all = get_wavevectors(L_mean);
#endif
	Structure_Factor Sq;
	init_structure_factor(wavevectors_all, Sq);
#endif

#if defined(gen_ISF) || defined(gen_SISF)
	std::vector<Vec3_wave_vector> wavevectors = max_wavevectors(L_mean);
	string outname_ISF = dir_out + "ISF_" + atom_kind + ".txt";
	string outname_SISF = dir_out + "SISF_" + atom_kind + ".txt";
#endif

#if defined(gen_force_chain) || defined(gen_bond)
	initialization_force();
#endif

#ifdef gen_eigen
	get_eigen(matrix_name, eigen_name);
#else


#ifdef gromacs_traj //得到N, Num_1, Num_2
	char filename_c[512];  //
	std::strcpy(filename_c, file_traj.c_str());
	XDRFILE* xtc = xdrfile_open(filename_c, "r");
	if (xtc == NULL) {
		printf("Error! Failed to open xtc file: %s\n", filename_c);
		return -1;
	}
	int read_return = read_xtc_natoms(filename_c, &N);
	if (read_return != 0) {
		printf("Error! reading number of atoms from xtc file.\n");
		return -1;
	}
	N_atom_mol = 15; //甲苯
	N_mol = N / N_atom_mol;
	rvec* posi_gromacs = (rvec*)calloc(N, sizeof(rvec)); //构型
	Num_1 = N_mol * 7;
	Num_2 = N_mol * 8;
#else
	ifstream read_traj(file_traj);
	std::cout << "File input: " << file_traj << endl;
	if (read_traj.fail()) {
		std::cout << ". File opening failed! " << endl;
		exit(-1);
	}
	Stat(file_traj, Num_1, Num_2, N);
#endif
	std::cout << "N: " << N << ".Num_1 : " << Num_1 << ".Num_2 : " << Num_2 << endl;

#if defined(gen_MSD) || defined(gen_SISF) || defined(gen_matrix) || defined(gen_trajectory) || defined(gen_momentum) || defined(gen_p_hop) || defined(gen_p_u) || defined(gen_chi4) || defined(gen_non_affine_D) //分析连续的轨迹
	int N_particle;
	if (atom_kind == "1") {
		N_particle = Num_1;
	}
	else if (atom_kind == "2") {
		N_particle = Num_2;
	}
	else if (atom_kind == "C2") { //目前只适用于gromacs
		N_particle = N_mol;
	}
	else {
		N_particle = N;
	}
	cout << "N_particle: " << N_particle << endl;
	Trajectory trajectory(total_config, Configuration(N_particle, Position(3, 0.0)));
	std::vector<vector<double>> L_traj(total_config, vector<double>(Dim, 0.0));
#ifdef gromacs_traj 
	read_gromacs_trajectory(file_traj, trajectory, L_traj);
#else
	read_lammps_trajectory(read_traj, L_traj, trajectory);
#endif
#ifdef gen_non_affine_D
	get_non_affine_D(trajectory, L_traj);
#endif
#ifdef gen_chi4
	get_chi4(trajectory, L_traj, atom_kind);
#endif
#ifdef gen_p_u
	get_p_u(trajectory, L_traj, atom_kind);
#endif
#ifdef gen_p_hop
	get_p_hop(trajectory, L_traj, atom_kind);
#endif
#ifdef gen_momentum
	get_momentum(trajectory, L_traj, atom_kind);
#endif
#ifdef gen_MSD
	get_MSD(trajectory, L_traj, atom_kind);
#endif
#ifdef gen_SISF 		
	vector<long double> F_t;
	F_t = get_SISF(trajectory, L_traj, wavevectors);
	std::ofstream out_SISF(outname_SISF, ios::out);
	for (int i = 0; i < F_t.size(); i++) {
		out_SISF << F_t[i] << endl;
	}
	out_SISF.close();
#endif
#ifdef gen_matrix
	write_matrix(trajectory, L_traj, atom_kind);
#endif

#ifdef gen_trajectory
	get_trajectory(trajectory, L_traj);
#endif

#else //分析每个时刻的构型

	vector<Vec3> Positions_ini, Positions_last(N), positions_NA_last(N), positions_NA_new(N);
	int number_sample = 0;
	double L_ini[Dim] = { 0.0 }, L[Dim] = { 0.0 }, L_last[Dim] = { 0.0 };
#ifdef gromacs_traj //跳过前ini_config个构型
	int step;
	float lambda, t;
	matrix box;
	for (int i_config = 0; i_config < ini_config; i_config++) {
		read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
		if (read_return != 0)
		{
			printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
			return -1;
		}
	}
#else
	skip_lammps_trajectory(read_traj);
#endif

	for (int i_config = 0; i_config < total_config; i_config++) {
		if ((i_config % d_S == i_d_S)) {
			vector<Vec3> Positions;
#ifdef gromacs_traj //得到当前帧的构型坐标和盒子尺寸
			read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
			if (read_return != 0) {
				printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
				return -1;
			}
			for (int j = 0; j < Dim; j++) {
				L[j] = box[j][j] * 10.0;
			}
			for (int i = 0; i < N; i++) {			// 处理坐标
				for (int j = 0; j < Dim; j++) {
					float x_shifted = posi_gromacs[i][j] * 10.0;
					while (x_shifted < 0.0) x_shifted = x_shifted + L[j];
					while (x_shifted >= L[j]) x_shifted = x_shifted - L[j];
					posi_gromacs[i][j] = x_shifted;
				}
			}

			for (int i = 0; i < N; i++) {
				int serial = i % atom_kind_Toluene.size();
				int kind = atom_kind_Toluene[serial];
				Positions.push_back({ posi_gromacs[i][0], posi_gromacs[i][1], posi_gromacs[i][2], kind, i });
			}
#else
			read_lammps_trajectory(read_traj, Positions, L);
#endif
			std::cout << dir_read << i_config << ". L_x: " << L[0] << endl;

			if (i_config == 0) {
				for (int j = 0; j < N; j++) {
					Positions_ini.push_back({ Positions[j].x, Positions[j].y , Positions[j].z, Positions[j].kind, Positions[j].serial });
				}
				for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
					L_ini[i_Dim] = L[i_Dim];
				}
			}

			number_sample++;
#ifdef gen_cluster
			get_cluster(Positions, L, i_config);
#endif
#ifdef gen_position
			get_position(Positions, L);
#endif
#ifdef gen_bond
			get_bond(Positions, L, whe_initial);
#endif
#ifdef gen_tetrahedron
			get_tetrahedron(Positions, L);
#endif

#ifdef gen_D_NA
			if (i_config > 0) {
				for (int i = 0; i < Positions.size(); ++i) { //减去仿射变换的偏移量
					positions_NA_new[i].x = Positions[i].x - ((L[0]) / L_last[0]) * Positions_last[i].x;
					positions_NA_new[i].y = Positions[i].y - ((L[1]) / L_last[1]) * Positions_last[i].y;
					positions_NA_new[i].z = Positions[i].z - ((L[2]) / L_last[2]) * Positions_last[i].z;
				}
				cout << L_last[0] << " " << Positions_last[0].x << " " << L[0] << " " << Positions[0].x << endl;
			}
			if (i_config > 1) {
				get_global_non_affine_D(Positions, L, positions_NA_last, positions_NA_new, Positions_ini, L_ini, i_config);
			}
			if (i_config > 0) {
				positions_NA_last = positions_NA_new;
			}
#endif
#ifdef gen_matrix_water
			get_matrix_water(Positions, L, GRID_SIZE);
#endif
#ifdef gen_rdf
			get_rdf(Positions, L, rdf_11, rdf_12, rdf_22, rdf_All);
#endif
#ifdef gen_sq
			get_structure_factor(Positions, L, wavevectors_all, Sq);
#endif
#ifdef gen_ISF 				
			vector<double> F_t_x(Dim, 0.0);
			vector<double> F_t_y(Dim, 0.0);
			get_ISF(Positions, F_t_x, F_t_y, wavevectors);
			std::ofstream out_ISF(outname_ISF, ios::app);
			out_ISF << F_t_x[0] << " " << F_t_x[1] << " " << F_t_x[2] << " " << F_t_y[0] << " " << F_t_y[1] << " " << F_t_y[2] << endl;
			out_ISF.close();
#endif
#ifdef gen_angle
			get_angle(Positions, L);
#endif
#ifdef gen_O_Si
			get_O_Si(Positions, L);
#endif
#ifdef gen_z
			get_z(Positions, L, Positions_ini, L_ini);
#endif
#ifdef gen_force_chain
			compute_force_chain(Positions, L);
#endif

			Positions_last = Positions;
			for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
				L_last[i_Dim] = L[i_Dim];
			}

		}
		else {
#ifdef gromacs_traj //跳过当前帧
			read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
			if (read_return != 0) {
				printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
				return -1;
			}
#else
			for (int j = 0; j < N + N_infor; j++) {
				getline(read_traj, tmp);
			}
#endif
		}
	}

#ifdef gen_rdf
	write_rdf(rdf_11, rdf_12, rdf_22, rdf_All, number_sample);
#endif
#ifdef gen_sq
	write_structure_factor(Sq, number_sample);
#endif

#endif

#ifdef gromacs_traj
	free(posi_gromacs);
	xdrfile_close(xtc);
#endif

#endif

	double minute = (static_cast<int>(time(NULL)) - Tbegin) / 60.0;
	std::cout << "Simulation time: " << minute << " min. " << endl;
	return 0;
}