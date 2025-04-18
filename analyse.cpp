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
#include "matrix_water.h"
#include "chi_4.h"
#include "cluster_NA.h"

int main(int argc, char* argv[])
{
	int Tbegin = static_cast<int> (time(NULL));
	functions = argv[1];
	if (functions == "toluene_MSD" || functions == "toluene_sq" || functions == "toluene_rdf" || functions == "toluene_SISF") {
		for (int i = 0; i < 24; i++) {
			Temp_list.emplace_back(to_string(150 + i*10));
		}
		atom_kind = argv[2];
		int i_temperature = stoi(argv[3]);
		total_config = stoi(argv[4]);
		string dir_save = argv[5];
		S_k_max = stof(argv[6]);
		traj_filename = "traj_comp.xtc";
		dir_read = "../../data/" + Temp_list[i_temperature] + "K/";
		dir_out = "../../data/" + dir_save + "/" + Temp_list[i_temperature] + "K_";
		std::cout << dir_read << ". dir_out: " << dir_out << endl;
	}
	else if (functions == "supercritical_matrix" || functions == "supercritical_eigen") {
		N = stoi(argv[2]);
		cout << N << endl;
		N_str = "N_" + to_string(N);
		Pressure = argv[3];
		Temperature = argv[4];
		Simu_case = argv[5];
		total_config = stoi(argv[6]);
		traj_filename = argv[7];
		out_label = argv[8];
		dir_read = "../../data/Widom_line/" + N_str + "/" + Pressure + "/" + Temperature + "/" + Simu_case + "/";
		dir_matrix = "../../data/Widom_line/matrix/";
		dir_eigen = "../../data/Widom_line/eigen/";
		matrix_name = dir_matrix + to_string(N) + "_" + Pressure + "_" + Temperature + "_" + out_label + ".txt";
		eigen_name = dir_eigen + to_string(N) + "_" + Pressure + "_" + Temperature + "_" + out_label + "_eigen.txt";
		std::cout << dir_read << ". traj_filename: " << traj_filename << ". matrix_name: " << matrix_name << ". eigen_name: " << eigen_name << endl;
		N_micro = N;
		const int GRID_SIZE = std::round(pow(N, 1.0 / 3.0)); //四舍五入
	}
	else if (functions == "Cu50Zr50") {
		//atom_kind = argv[2];
		//dir_data = argv[3];
		//S_k_max = stof(argv[4]);
		total_config = 1001;
		ini_config = 0;
		dir_read = "../../data/Cu50Zr50/";
		dir_out = dir_read;
		traj_filename = "new_long_eq_traj_T900.atom";
		std::cout << dir_out << endl;
	}
	else if (functions == "SiO2_local_NA" || functions == "SiO2_p_u" || functions == "SiO2_NA" || functions == "SiO2_posi" || functions == "SiO2_z" || functions == "SiO2_cluster") { //SiO2拉伸
		ini_config = 0;
		atom_kind = argv[2];
		dir_system = argv[3];
		total_config = stoi(argv[4]);
		dt = stoi(argv[5]);
		O_O_d2 = stof(argv[6]);
		first_layer = argv[7];
		min_serial = stoi(argv[8]);
		string dir_save = argv[9];
		O_O_d2 = O_O_d2 * O_O_d2;
		dir_out = "../../data/" + dir_save + "/" + dir_system + "_";
		dir_read = "../../data/" + dir_system + "/";
		std::cout << "dir_out: " << dir_out << ". dir_read: " << dir_read << ". dt: " << dt << ". O_O_d2: " << O_O_d2 << ". min_serial: " << min_serial << endl;
	}
	else if (functions == "AlSm_momentum" || functions == "AlSm_chi4") {
		ini_config = 0;
		atom_kind = argv[2];
		dir_system = argv[3];
		total_config = stoi(argv[4]);
		traj_filename = argv[5];	
		Al_Al_d = stof(argv[6]) / 3.0;
		string dir_save = argv[7];
		dir_read = "../../data/" + dir_system + "/";
		dir_out = "../../data/" + dir_save + "/" + dir_system + "_";
		std::cout << "dir_out: " << dir_out << ". dir_read: " << dir_read << ". Al_Al_d: " << Al_Al_d << endl;
	}
	else if (functions == "AlSm_matrix" || functions == "AlSm_eigen" || functions == "AlSm_posi") { 
		atom_kind = argv[2];
		Temperature = argv[3];
		S_k_max = stof(argv[4]);
		smoothing = (argv[5]);
		N_micro = 6 * 21 * 21 * 21;
		micro_form = "6";
		ini_config = 0;
		total_config = 10000 - ini_config;
		dir_read = "../../data/" + Temperature + "/";
		std::cout << "atom_kind: " << atom_kind << ", time_interval: " << time_interval << endl;		
		decompose_way = "none";
		dir_matrix = "../../data/matrix/";
		dir_eigen = "../../data/eigen/";
		matrix_name = dir_matrix + Temperature + "_" + smoothing + ".txt";
		eigen_name = dir_eigen + Temperature + "_" + smoothing + "_eigen.txt";
		std::cout << dir_read << ". traj_filename: " << traj_filename << ". matrix_name: " << matrix_name << ". eigen_name: " << eigen_name << endl;
	}
	else if (functions == "supercooled") { //supercooled
		dir_system = argv[2];
		dir_data = argv[3];
		dir_read = "../../data/raw/supercooled/" + dir_data + "/";
		dir_out = dir_read;
		M = stoi(argv[4]); //构型数
		N_micro = stoi(argv[5]); //矩阵的微观态数量
		micro_form = argv[6]; //矩阵形式
		part = argv[7]; //前半部分，后半部分，全部
		matrix_kind = argv[8]; //矩阵的种类
		finite_size = "non_finite_size"; //是否研究有限尺寸
		if_relative = "relative"; //是否减去平均值
		matrix_name = "lattice_Gauss_vel_rot_Theta"; //矩阵名称  microstatepar_vel_rot_Theta
		decompose_way = "SVD"; //分解方式
		for (int i = 0; i < 10; i++) {Temp_list.emplace_back(to_string(i * 10 + 200));}
		Temp_list.emplace_back(to_string(235));
		std::cout << dir_read << ". M: " << M << ". N_micro: " << N_micro << ". micro_form: " << micro_form << ". part: " << part << ". matrix_kind: " << matrix_kind << ". matrix_name: " << matrix_name << endl;
	}
	else if (functions == "super_critical") {
		dir_system = argv[2];
		dir_data = argv[3];
		dir_read = "../../data/raw/" + dir_data + "/";
		dir_out = dir_read;
		M = stoi(argv[4]); //构型数
		N_micro = stoi(argv[5]); //矩阵的微观态数量
		micro_form = argv[6]; //矩阵形式
		finite_size = argv[7]; //是否研究有限尺寸
		if_relative = argv[8]; //是否减去平均值
		matrix_name = argv[9]; //矩阵名称
		decompose_way = "Eigen"; //分解方式
		for (int i = 0; i < 41; i++) {
			Temp_list.emplace_back(to_string(i * 5 + 600));
		}
		vector<string> Temp_list_175 = { "648.7", "648.9", "649.1", "649.3", "649.5", "649.7", "649.9", "650.1", "650.3", "651", "652", "653", "654" };
		if (dir_data == "175bar") {
			for (int i = 0; i < Temp_list_175.size(); i++) {
				Temp_list.emplace_back(Temp_list_175[i]);
			}
		}
		std::cout << dir_read << ". M: " << M << ". N_micro: " << N_micro << ". micro_form: " << micro_form << ". finite_size: " << finite_size << ". if_relative: " << if_relative << ". matrix_name: " << matrix_name << endl;
	}
	else if (functions == "super_critical_finite_size") { //sbatch critical_finite_size.sh "1000" "3375" "8000" "15625" "64000"		
		dir_system = argv[2];
		dir_data = argv[3]; //压强
		dir_read = "../../data/raw/WL/";
		dir_out = dir_read;
		M = stoi(argv[4]); //构型数
		N_micro = stoi(argv[5]); //矩阵的微观态数量
		micro_form = argv[6]; //矩阵形式
		finite_size = argv[7]; //是否研究有限尺寸
		if_relative = argv[8]; //是否减去平均值
		matrix_name = argv[9]; //矩阵名称
		decompose_way = argv[10]; //分解方式
		std::cout << dir_read << ". M: " << M << ". N_micro: " << N_micro << ". micro_form: " << micro_form << ". finite_size: " << finite_size << ". if_relative: " << if_relative << ". matrix_name: " << matrix_name
			<< ". decompose_way: " << decompose_way << endl;
		int isCreate = createDir(dir_out);
		if (isCreate == 0) {
			cout << "created path: " << dir_out << endl;
		}
		else {
			printf("error: create path failed!\n");
			exit(-1);
		}
	}
	//total_config = CountLines((char*)file.c_str()) / (N + N_infor);
	//total_config = count_lines(file_traj) / (N + N_infor);
	const string file_traj = dir_read + traj_filename;

	std::cout << "ini_config:" << ini_config << ". total_config: " << total_config << ". d_S: " << d_S << ". file_traj: " << file_traj << endl;

#ifdef gen_rdf
	vector<double> rdf_11(N_Bins_r, 0.0); //11
	vector<double> rdf_12(N_Bins_r, 0.0); //12 21
	vector<double> rdf_22(N_Bins_r, 0.0); //22
	vector<double> rdf_All(N_Bins_r, 0.0);
#endif

#if defined(gen_sq) || defined(gen_ISF) || defined(gen_SISF)
	double L_mean[Dim] = { 0.0 };
	get_L_mean(file_traj, L_mean); //得到平均的盒子尺寸
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
	string outname_ISF = dir_out + "ISF.txt";
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
	vector<double> F_t;
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