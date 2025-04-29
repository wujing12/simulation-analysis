#pragma once
inline void Stat(string file_traj, int& Num_1, int& Num_2, int& N) {
	int serial;
	ifstream read_traj(file_traj);
	for (int j = 0; j < 3; j++) {
		getline(read_traj, tmp);
	}
	read_traj >> N;
	std::cout << "N in traj:" << N << endl;
	read_traj.seekg(0, ios::beg);
	vector<int> element(N);
	for (int j = 0; j < N_infor; j++) {
		getline(read_traj, tmp);
	}
	for (int j = 0; j < N; j++) {
		read_traj >> serial >> element[j];
		getline(read_traj, tmp);
	}

	Num_1 = 0;
	Num_2 = 0;
	for (int j = 0; j < N; j++) {
		if (element[j] == 1) {
			Num_1 += 1;
		}
		else {
			Num_2 += 1;
		}
	}
	read_traj.close();
}
inline void get_L_mean(const string& file_traj, double(&L_mean)[Dim]) {
#ifdef gromacs_traj
	char filename_c[512];
	std::strcpy(filename_c, file_traj.c_str());

	// 打开 .xtc 文件
	XDRFILE* xtc = xdrfile_open(filename_c, "r");
	if (xtc == NULL) {
		printf("Failed to open xtc file: %s\n", filename_c);
		return;
	}
	int step;
	float lambda, t;
	matrix box;
	int read_return;
	rvec* posi_gromacs = (rvec*)calloc(N, sizeof(rvec));
	for (int i_config = 0; i_config < ini_config; i_config++) {
		read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
		if (read_return != 0)
		{
			printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
			exit(-1);
		}
	}
	for (int i_config = 0; i_config < total_config; i_config++) {
		// 读取每个配置
		read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
		if (read_return != 0) {
			printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
			exit(-1);
		}
		for (int j = 0; j < 3; j++) {
			L_mean[j] += box[j][j] * 10.0;
		}
	}

	// 释放内存
	free(posi_gromacs);
	// 关闭 xtc 文件
	xdrfile_close(xtc);
#else

	ifstream read_traj(file_traj);
	for (int i = 0; i < ini_config; i++) { //读取前置的不需要的构型
		for (int j = 0; j < N + N_infor; j++) {
			getline(read_traj, tmp);
		}
	}
	for (int i = 0; i < total_config; i++) {
		for (int j = 0; j < N_head; j++) {
			getline(read_traj, tmp);
		}
		double Lleft, Lright;
		for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
			read_traj >> Lleft >> Lright;
			L_mean[i_Dim] += Lright - Lleft;
		}
		getline(read_traj, tmp);
		getline(read_traj, tmp);
		for (int j = 0; j < N; j++) {
			getline(read_traj, tmp);
		}
	}
	read_traj.close();
#endif

	for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
		L_mean[i_Dim] = L_mean[i_Dim] / total_config;
		std::cout << L_mean[i_Dim] << " ";
	}
	std::cout << endl;

}
inline void skip_lammps_trajectory(std::ifstream& read_traj) {
	for (int i = 0; i < ini_config; i++) { //读取前置的不需要的构型
		for (int j = 0; j < N + N_infor; j++) {
			getline(read_traj, tmp);
		}
	}
}
inline void read_gromacs_trajectory(const string& file_traj, Trajectory& trajectory, std::vector<vector<double>>& L_traj) {
#ifdef gromacs_traj

	char filename_c[512];
	std::strcpy(filename_c, file_traj.c_str());

	// 打开 .xtc 文件
	XDRFILE* xtc = xdrfile_open(filename_c, "r");
	if (xtc == NULL) {
		printf("Failed to open xtc file: %s\n", filename_c);
		return;
	}
	int step;
	float lambda, t;
	matrix box;
	int read_return;
	rvec* posi_gromacs = (rvec*)calloc(N, sizeof(rvec));
	for (int i_config = 0; i_config < ini_config; i_config++) {
		read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
		if (read_return != 0)
		{
			printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
			exit(-1);
		}
	}
	int count = 0;
	for (int i_config = 0; i_config < total_config; i_config++) {
		// 读取每个配置
		read_return = read_xtc(xtc, N, &step, &t, box, posi_gromacs, &lambda);
		if (read_return != 0) {
			printf("Error! xtc file only contains %d configurations \n ", i_config + 1);
			exit(-1);
		}
		if ((i_config % d_S == i_d_S)) {
			for (int j = 0; j < Dim; j++) {
				L_traj[count][j] = box[j][j] * 10.0;
			}
			std::cout << dir_read << count << ". L_x: " << L_traj[count][0] << endl;
			// 处理坐标
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < Dim; j++) {
					float x_shifted = posi_gromacs[i][j] * 10.0;
					while (x_shifted < 0.0) x_shifted = x_shifted + L_traj[count][j];
					while (x_shifted >= L_traj[count][j]) x_shifted = x_shifted - L_traj[count][j];
					posi_gromacs[i][j] = x_shifted;
				}
			}

			// 存储处理后的位置
			vector<Vec3> Positions_1, Positions_2, Positions_All;

			for (int i = 0; i < N; i++) {
				int serial = i % atom_kind_Toluene.size();
				int kind = atom_kind_Toluene[serial];
				Positions_All.push_back({ posi_gromacs[i][0], posi_gromacs[i][1], posi_gromacs[i][2], kind, i });
				if (kind == 1) {
					Positions_1.push_back({ posi_gromacs[i][0], posi_gromacs[i][1], posi_gromacs[i][2], kind, i });
				}
				else if (kind == 2) {
					Positions_2.push_back({ posi_gromacs[i][0], posi_gromacs[i][1], posi_gromacs[i][2], kind, i });
				}
			}
			if (atom_kind == "1") {
				for (int j = 0; j < Num_1; j++) {
					trajectory[count][j][0] = Positions_1[j].x;
					trajectory[count][j][1] = Positions_1[j].y;
					trajectory[count][j][2] = Positions_1[j].z;
				}
			}
			else if (atom_kind == "2") {
				for (int j = 0; j < Num_2; j++) {
					trajectory[count][j][0] = Positions_2[j].x;
					trajectory[count][j][1] = Positions_2[j].y;
					trajectory[count][j][2] = Positions_2[j].z;
				}
			}
			else if (atom_kind == "C2") {
				for (int j = 0; j < N_mol; j++) {
					trajectory[count][j][0] = Positions_All[j * N_atom_mol + 1].x;
					trajectory[count][j][1] = Positions_All[j * N_atom_mol + 1].y;
					trajectory[count][j][2] = Positions_All[j * N_atom_mol + 1].z;
				}
			}

			else {
				for (int j = 0; j < N; j++) {
					trajectory[count][j][0] = Positions_All[j].x;
					trajectory[count][j][1] = Positions_All[j].y;
					trajectory[count][j][2] = Positions_All[j].z;
				}
			}
			count++;
		}
	}

	// 释放内存
	free(posi_gromacs);
	// 关闭 xtc 文件
	xdrfile_close(xtc);
#endif
}
inline void read_lammps_trajectory(std::ifstream& read_traj, std::vector<vector<double>>& L_traj, Trajectory& trajectory) {
	skip_lammps_trajectory(read_traj);
	int count = 0, kind, serial;
	for (int i_config = 0; i_config < total_config; i_config++) {
		if ((i_config % d_S == i_d_S)) {
			
			for (int j = 0; j < N_head; j++) {
				getline(read_traj, tmp);
			}
			double Lleft, Lright;
			double L[Dim] = { 0.0 };
			for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
				read_traj >> Lleft >> Lright;
				L[i_Dim] = Lright - Lleft;
				L_traj[count][i_Dim] = L[i_Dim];
			}
			std::cout << dir_read << count << ". L_x: " << L_traj[count][0] << endl;
			getline(read_traj, tmp);
			getline(read_traj, tmp);
			vector<Vec3> Positions;
			double posi[Dim] = { 0.0 };
			for (int j = 0; j < N; j++) {
				read_traj >> serial >> kind; //每个原子的编号和种类
				for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
					read_traj >> posi[i_Dim];
					posi[i_Dim] = posi[i_Dim] * L[i_Dim];
				}
				//std::cout << posi[0] << " " << posi[1] << " " << posi[2] << endl;
				PBC_r(posi, L);
				Positions.push_back({ posi[0], posi[1], posi[2], kind, serial }); //以盒子的中心为质心
			}
			getline(read_traj, tmp);
			std::sort(Positions.begin(), Positions.end(), compareByid); //粒子根据序号进行排序
			if (atom_kind == "1") {
				for (int j = 0; j < Num_1; j++) {
					trajectory[count][j][0] = Positions[j].x;
					trajectory[count][j][1] = Positions[j].y;
					trajectory[count][j][2] = Positions[j].z;
				}
			}
			else if (atom_kind == "2") {
				for (int j = Num_1; j < N; j++) {
					trajectory[count][j - Num_1][0] = Positions[j].x;
					trajectory[count][j - Num_1][1] = Positions[j].y;
					trajectory[count][j - Num_1][2] = Positions[j].z;
				}
			}
			else {
				for (int j = 0; j < N; j++) {
					trajectory[count][j][0] = Positions[j].x;
					trajectory[count][j][1] = Positions[j].y;
					trajectory[count][j][2] = Positions[j].z;
				}
			}
			count++;
		}
		else {
			for (int j = 0; j < N + N_infor; j++) {
				getline(read_traj, tmp);
			}
		}
	}
	read_traj.close();
}

inline void read_lammps_trajectory(std::ifstream& read_traj, vector<Vec3>& Positions, double(&L)[Dim]) {
	int serial, kind;
	for (int j = 0; j < N_head; j++) {
		getline(read_traj, tmp);
	}
	double Lleft, Lright;
	for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
		read_traj >> Lleft >> Lright;
		L[i_Dim] = Lright - Lleft;
	}
	getline(read_traj, tmp);
	getline(read_traj, tmp);
	double posi[Dim] = { 0.0 };
	for (int j = 0; j < N; j++) {
		read_traj >> serial >> kind; //每个原子的编号和种类
		for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
			read_traj >> posi[i_Dim];
			posi[i_Dim] = posi[i_Dim] * L[i_Dim];
		}
		PBC_r(posi, L);
		Positions.push_back({ posi[0], posi[1], posi[2], kind, serial }); //以盒子的中心为质心
	}
	getline(read_traj, tmp);
	std::sort(Positions.begin(), Positions.end(), compareByid); //粒子根据序号进行排序

}