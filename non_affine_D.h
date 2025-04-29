#pragma once
constexpr int N_Bins_D = 200;
constexpr double cut_D = 0.11;
//==============================================================================
// 构建邻居列表：基于 t0 时刻的粒子位置和盒子尺寸
//==============================================================================
void build_neighbors_SL(const vector<vector<double>>& positions, const vector<double>& L, vector<vector<Line>>& neighbors, const int& i) {
	if (i < Num_1) { //氧原子
		for (size_t j = 0; j < positions.size(); ++j) {
			if (i != j) {
				double dist2 = Distance_2(positions[i], positions[j], L);
				//cout << i << " " << j << " " << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << " " << positions[j][0] << " " << positions[j][1] << " " << positions[j][2] << " " << dist2 << " " << distance_cut << endl;
				if (dist2 <= distance_cut) {
					neighbors[i].push_back(Line(j, dist2));
				}
			}
		}
	}
	else { //硅原子
		for (size_t j = 0; j < positions.size(); ++j) {
			if (i != j) {
				double dist2 = Distance_2(positions[i], positions[j], L);
				if (dist2 <= Si_Si_d2) {
					neighbors[i].push_back(Line(j, dist2));
				}
			}
		}
	}
}
void build_neighbors_FL(const vector<vector<double>>& positions, const vector<double>& L, vector<vector<Line>>& neighbors, const double& cut, const int& i) {
	for (size_t j = 0; j < positions.size(); ++j) {
		if (i != j) {
			double dist2 = Distance_2(positions[i], positions[j], L);
			//cout << i << " " << j << " " << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << " " << positions[j][0] << " " << positions[j][1] << " " << positions[j][2] << " " << dist2 << " " << distance_cut << endl;
			if (dist2 <= cut) {
				neighbors[i].push_back(Line(j, dist2));
			}
		}
	}
}
//==============================================================================
// 计算局部变形梯度：给定 t0 和 t1 时刻的相对位移，计算 J = C⁻¹ * B（当 C 非奇异时）
//==============================================================================
Eigen::Matrix3d compute_deformation_gradient(const vector<vector<double>>& d_r0, const vector<vector<double>>& d_r) {
	Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
	for (size_t j = 0; j < d_r0.size(); ++j) {
		for (int a = 0; a < Dim; ++a) {
			for (int b = 0; b < Dim; ++b) {
				B(a, b) += d_r0[j][a] * d_r[j][b];
				C(a, b) += d_r0[j][a] * d_r0[j][b];
			}
		}
	}
	Eigen::Matrix3d C_inv;
	double epsilon = 1e-8;
	if (std::abs(C.determinant()) < epsilon) {
		// 使用伪逆或添加正则化项
		//C_inv = (C + epsilon * Eigen::Matrix3d::Identity()).inverse();
		Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix3d> cod(C);
		C_inv = cod.pseudoInverse();
	}
	else {
		C_inv = C.inverse();
	}
	return C_inv * B;
}

//==============================================================================
// 计算单个粒子的非仿射位移：给定 t0 和 t1 时刻的位置、邻居列表以及盒子尺寸, 返回值：非仿射位移 D 以及各方向分量 D_d（数组大小为 Dim）
//==============================================================================
void compute_non_affine_disp_for_particle(const vector<vector<double>>& pos_t0,
	const vector<vector<double>>& pos_t1,
	int i,
	const vector<vector<Line>>& neighbors,
	double& D,
	double D_d[Dim], 
	const vector<double>& L_t0,
	const vector<double>& L_t1)
{
	size_t neigh_num = neighbors[i].size();

	if (neigh_num <= 3) { //D_min = 0
		D = 0.0;
		for (int d = 0; d < Dim; ++d)
			D_d[d] = 0.0;
		return;
	}
	// 存储相对位移：d_r0：t0时刻，d_r：t1时刻，d_rt：仿射预测值
	vector<vector<double>> d_r0(neigh_num, vector<double>(Dim, 0.0));
	vector<vector<double>> d_r(neigh_num, vector<double>(Dim, 0.0));
	vector<vector<double>> d_rt(neigh_num, vector<double>(Dim, 0.0));
	double L_t1_list[Dim] = { L_t1[0], L_t1[1], L_t1[2] };
	double L_t0_list[Dim] = { L_t0[0], L_t0[1], L_t0[2] };
	double L_t1_h[Dim] = { L_t1[0] / 2.0, L_t1[1] / 2.0, L_t1[2] / 2.0 };
	double L_t0_h[Dim] = { L_t0[0] / 2.0, L_t0[1] / 2.0, L_t0[2] / 2.0 };
	// 计算相对位移
	for (size_t j = 0; j < neigh_num; ++j) {
		int nb_idx = neighbors[i][j].seq;
		double dr0[Dim] = { pos_t0[nb_idx][0] - pos_t0[i][0],
							pos_t0[nb_idx][1] - pos_t0[i][1],
							pos_t0[nb_idx][2] - pos_t0[i][2] };
		double dr[Dim] = { pos_t1[nb_idx][0] - pos_t1[i][0],
							pos_t1[nb_idx][1] - pos_t1[i][1],
							pos_t1[nb_idx][2] - pos_t1[i][2] };

		PBC(L_t0_h, L_t0_list, dr0);
		PBC(L_t1_h, L_t1_list, dr);
		for (int d = 0; d < Dim; ++d) {
			d_r0[j][d] = dr0[d];
			d_r[j][d] = dr[d];
		}
	}

	// 计算局部变形梯度 J
	Eigen::Matrix3d J = compute_deformation_gradient(d_r0, d_r);

	// 计算仿射预测的相对位移 d_rt = J * d_r0
	for (size_t j = 0; j < neigh_num; ++j) {
		for (int d = 0; d < Dim; ++d) {
			d_rt[j][d] = J(0, d) * d_r0[j][0] +
				J(1, d) * d_r0[j][1] +
				J(2, d) * d_r0[j][2];
		}
	}

	// 计算非仿射位移 D 及其各方向分量 D_d
	D = 0.0;
	for (int d = 0; d < Dim; ++d)
		D_d[d] = 0.0;
	for (size_t j = 0; j < neigh_num; ++j) {
		for (int d = 0; d < Dim; ++d) {
			double diff = d_rt[j][d] - d_r[j][d];
			D += diff * diff;
			D_d[d] += diff * diff;
		}
	}
	D /= neigh_num;
	for (int d = 0; d < Dim; ++d)
		D_d[d] /= neigh_num;

}
inline void change(double& D,
	double D_d[Dim], int i,
	vector<vector<double>>& D_kind) {
	D = std::max(D, 1e-6);
	D_d[0] = std::max(D_d[0], 1e-6);
	D_d[1] = std::max(D_d[1], 1e-6);
	D_d[2] = std::max(D_d[2], 1e-6);
	D_kind[0][i] = log(D);
	D_kind[1][i] = log(D_d[0]);
	D_kind[2][i] = log(D_d[1]);
	D_kind[3][i] = log(D_d[2]);
}
// 计算并保存统计数据的函数
inline void write_distribution(vector<vector<double>>& D_data, const string& outname, int N_Bins_D) {
	for (int i = 0; i < 4; ++i) {
		distribution(D_data[i], outname, N_Bins_D);
	}
}
//==============================================================================
// 主函数：根据 traj 计算非仿射位移
//==============================================================================
inline void get_non_affine_D(const Trajectory& traj, const vector<vector<double>>& L_traj) {
	string layer = "FL";
	if (!first_layer) {
		layer = "SL";
	}
	// 打开原子级和统计输出文件
	std::ofstream out_file(dir_out + to_string(Num_1) + "_dt" + to_string(dt) + "_" + layer + ".txt", ios::out);
	std::ofstream out_file_kind(dir_out + to_string(Num_1) + "_dt" + to_string(dt) + "_" + layer + "_kind.txt", ios::out);
	string outname1 = dir_out + "kind1_dt" + to_string(dt) + "_" + layer + "_sta.txt";
	string outname2 = dir_out + "kind2_dt" + to_string(dt) + "_" + layer + "_sta.txt";
	std::ofstream out_file1(outname1, ios::out);	
	std::ofstream out_file2(outname2, ios::out);
	out_file1.close();
	out_file2.close();
	if (dt > 0) {
		cout << "t0 is ergodic. Num_1:" << Num_1 << ", Num_2:" << Num_2 << ", N:" << traj[0].size() << endl;
		vector<vector<Line>> neighbors(traj[0].size());
		// 遍历所有 t₀ 时刻（确保 t₀ + delta_index 有效）
		for (size_t t0 = 0; t0 < traj.size() - dt; t0 += dt) {
			vector<vector<double>> D_1(4, vector<double>(Num_1, 0.0));
			vector<vector<double>> D_2(4, vector<double>(Num_2, 0.0));
			// 提取 t₀ 和 t₁ 时刻的粒子数据及盒子尺寸
			const auto& pos_t0 = traj[t0];
			const auto& pos_t1 = traj[t0 + dt];
			neighbors.clear(); 
			// 对 t₀ 时刻的每个粒子计算非仿射位移
			for (size_t i = Num_1; i < pos_t0.size(); ++i) {
				// 更新 t₀ 时刻的邻居列表
				if (first_layer) {
					build_neighbors_FL(pos_t0, L_traj[t0], neighbors, Si_O_d2, i);
				}else{
					build_neighbors_SL(pos_t0, L_traj[t0], neighbors, i);
				}
				double D = 0.0, D_d[Dim] = { 0.0 };
				compute_non_affine_disp_for_particle(pos_t0, pos_t1, i, neighbors, D, D_d, L_traj[t0], L_traj[t0 + dt]);
				out_file << D << " ";
				if (D > cut_D) {
					out_file_kind << "1 ";
				}
				else {
					out_file_kind << "0 ";
				}
				if (i < Num_1) { //氧原子
					change(D, D_d, i, D_1);
				}
				else { //硅原子
					change(D, D_d, i - Num_1, D_2);
				}
			}
			out_file_kind << endl;
			out_file << endl;
			// 写入总体统计数据
			write_distribution(D_1, outname1, N_Bins_D);
			write_distribution(D_2, outname2, N_Bins_D);
		}
	}
	else {
		cout << "t1 is ergodic. Num_1:" << Num_1 << ", Num_2:" << Num_2 << ", N:" << traj[0].size() << endl;
		const size_t t0 = 0;
		const auto& pos_t0 = traj[t0];
		vector<vector<Line>> neighbors(pos_t0.size());
		// 构建 t₀ 时刻的邻居列表
		for (size_t i = Num_1; i < pos_t0.size(); ++i) {			
			if (first_layer) {
				build_neighbors_FL(pos_t0, L_traj[t0], neighbors, Si_O_d2, i);
			}
			else {
				build_neighbors_SL(pos_t0, L_traj[t0], neighbors, i);
			}
		}
		// 遍历所有 t1时刻
		for (size_t t1 = 1; t1 < traj.size(); t1 += 1) {
			vector<vector<double>> D_1(4, vector<double>(Num_1, 0.0));
			vector<vector<double>> D_2(4, vector<double>(Num_2, 0.0));
			// 提取t₁ 时刻的粒子数据及盒子尺寸
			const auto& pos_t1 = traj[t1];
			// 对 t₀ 时刻的每个粒子计算非仿射位移
			for (size_t i = Num_1; i < pos_t0.size(); ++i) {
				double D = 0.0, D_d[Dim] = { 0.0 };
				compute_non_affine_disp_for_particle(pos_t0, pos_t1, i, neighbors, D, D_d, L_traj[t0], L_traj[t1]);
				//out_file << D_d[0] << " " << D_d[1] << " " << D_d[2] << endl;
				out_file << D << " ";
				if (D > cut_D) {
					out_file_kind << "1 ";
				}
				else {
					out_file_kind << "0 ";
				}
				if (i < Num_1) {
					change(D, D_d, i, D_1);
				}
				else {
					change(D, D_d, i - Num_1, D_2);
				}
			}
			out_file_kind << endl;
			out_file << endl;
			// 写入总体统计数据
			write_distribution(D_1, outname1, N_Bins_D);
			write_distribution(D_2, outname2, N_Bins_D);
		}
	}
	out_file_kind.close();
	out_file.close();
}

inline void get_global_non_affine_D(const std::vector<Vec3>& positions,
	const double(&L)[Dim],
	const std::vector<Vec3>& positions_NA_last, 
	const std::vector<Vec3>& positions_NA_new,
	const std::vector<Vec3>& positions_ini,
	const double(&L_ini)[Dim],
	const int i_config) {

	string outname_1 = dir_out + "atomic_D_NA_1.txt"; //储存1粒子的非仿射位移
	string outname_2 = dir_out + "atomic_D_NA_2.txt"; //储存2粒子的非仿射位移
	string outname_3 = dir_out + "atomic_p_hop_1.txt"; //储存1粒子的重排
	string outname_4 = dir_out + "atomic_p_hop_2.txt"; //储存2粒子的重排
	std::ofstream out_file_1(outname_1, ios::app);
	std::ofstream out_file_2(outname_2, ios::app);
	std::ofstream out_file_3(outname_3, ios::app);
	std::ofstream out_file_4(outname_4, ios::app);
	double D_NA = zero, p_hop = zero, dr[Dim] = { zero }, dr_ini[Dim] = { zero };
	for (int i = 0; i < positions.size(); ++i) { //减去仿射变换的偏移量
		dr[0] = positions_NA_new[i].x - positions_NA_last[i].x;
		dr[1] = positions_NA_new[i].y - positions_NA_last[i].y;
		dr[2] = positions_NA_new[i].z - positions_NA_last[i].z;
		p_hop = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
		dr_ini[0] = positions[i].x - ((L[0]) / L_ini[0]) * positions_ini[i].x;
		dr_ini[1] = positions[i].y - ((L[1]) / L_ini[1]) * positions_ini[i].y;
		dr_ini[2] = positions[i].z - ((L[2]) / L_ini[2]) * positions_ini[i].z;
		D_NA = dr_ini[0] * dr_ini[0] + dr_ini[1] * dr_ini[1] + dr_ini[2] * dr_ini[2];
		
		if (positions[i].kind == 1) {
			out_file_1 << D_NA << " ";
			out_file_3 << p_hop << " ";
		}
		else {
			out_file_2 << D_NA << " ";
			out_file_4 << p_hop << " ";
		}
	}
	out_file_1 << endl;
	out_file_2 << endl;
	out_file_3 << endl;
	out_file_4 << endl;
	out_file_1.close();
	out_file_2.close();
	out_file_3.close();
	out_file_4.close();
}
inline void get_p_hop(Trajectory& trajectory, const vector<vector<double>>& L_traj, const string& atom_kind) {
	unwrap_positions(trajectory, L_traj);
	size_t numConfigs = trajectory.size(); // 构型数量
	size_t numParticles = trajectory[0].size(); // 每个构型的粒子数量
	const int t_R = 20, t_R_h = t_R/2, t_R_data = t_R + 1, t_R_h_data = t_R_h + 1;
	string outname_1 = dir_out + "p_hop_" + atom_kind + ".txt";
	string outname_2 = dir_out + "D_NA_" + atom_kind + ".txt";
	std::ofstream outputFile_1(outname_1, ios::out);
	std::ofstream outputFile_2(outname_2, ios::out);
	for (size_t t = t_R_h + 1; t < numConfigs - t_R_h; ++t) { //间隔的构型数量
		Trajectory trajectory_NA(t_R_data, Configuration(numParticles, Position(3, 0.0))); //非仿射位移
		int t_last, t_new;
		for (size_t t_traj = 0; t_traj < t_R_data; ++t_traj) {
			t_last = t + t_traj - t_R_h - 1;
			t_new = t_last + 1;
			for (size_t i = 0; i < numParticles; ++i) {
				for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
					trajectory_NA[t_traj][i][i_dim] = trajectory[t_new][i][i_dim] - ((L_traj[t_new][i_dim]) / L_traj[0][i_dim]) * trajectory[0][i][i_dim];
				}
			}
		}
		for (size_t i = 0; i < numParticles; ++i) {
			double r_last_m[Dim] = { 0.0 }, r_new_m[Dim] = { 0.0 };
			for (size_t t_traj = 0; t_traj < t_R_h_data; ++t_traj) {
				for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
					r_last_m[i_dim] += trajectory_NA[t_traj][i][i_dim] / t_R_h_data;
					r_new_m[i_dim] += trajectory_NA[t_traj + t_R_h][i][i_dim] / t_R_h_data;
				}
			}
			double p_hop = 0.0, dr_last_sum = { 0.0 }, dr_new_sum = { 0.0 }, dr_last[Dim] = { 0.0 }, dr_new[Dim] = { 0.0 };
			for (size_t t_traj = 0; t_traj < t_R_h_data; ++t_traj) {
				for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
					dr_last[i_dim] = (trajectory_NA[t_traj][i][i_dim] - r_new_m[i_dim]);
					dr_new[i_dim] = (trajectory_NA[t_traj + t_R_h][i][i_dim] - r_last_m[i_dim]);
				}
				dr_last_sum += (dr_last[0] * dr_last[0] + dr_last[1] * dr_last[1] + dr_last[2] * dr_last[2]);
				dr_new_sum += (dr_new[0] * dr_new[0] + dr_new[1] * dr_new[1] + dr_new[2] * dr_new[2]);
			}
			p_hop = sqrt(dr_last_sum * dr_new_sum) / (t_R_h_data);
			outputFile_1 << p_hop << " ";
		}
		outputFile_1 << endl;
		double D_NA = 0.0;
		double dr_ini[Dim] = { 0.0 };
		for (size_t i = 0; i < numParticles; ++i) {
			for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
				dr_ini[i_dim] = trajectory[t][i][i_dim] - ((L_traj[t][i_dim]) / L_traj[0][i_dim]) * trajectory[0][i][i_dim];
			}
			D_NA = sqrt(dr_ini[0] * dr_ini[0] + dr_ini[1] * dr_ini[1] + dr_ini[2] * dr_ini[2]);
			outputFile_2 << D_NA << " ";
		}
		outputFile_2 << endl;
	}
	outputFile_1.close();
	outputFile_2.close();
}
inline void get_positions_NA(Trajectory& trajectory, const vector<vector<double>>& L_traj) {
	size_t T = trajectory.size();      // 构型数量
	size_t N = trajectory[0].size();  // 粒子数
	for (size_t i = 0; i < N; ++i) {
		for (size_t t = 0; t < T; ++t) {
			for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
				trajectory[t][i][i_dim] = trajectory[t][i][i_dim] - ((L_traj[t][i_dim]) / L_traj[0][i_dim]) * trajectory[0][i][i_dim];
			}
		}
	}
}
inline void get_p_u(Trajectory& trajectory, const vector<vector<double>>& L_traj, const string& atom_kind) {
	unwrap_positions(trajectory, L_traj);
	get_positions_NA(trajectory, L_traj);
	int N_Bins = 200;
	size_t numConfigs = trajectory.size(); // 构型数量
	size_t numParticles = trajectory[0].size(); // 每个构型的粒子数量
	string outname_1 = dir_out + "u_" + atom_kind + "_" + to_string(dt) + ".txt";
	string outname_2 = dir_out + "u_hist_" + atom_kind + "_" + to_string(dt) + ".txt";
	std::ofstream outputFile_1(outname_1, ios::out);
	reset_file(outname_2);
	for (size_t t = dt; t < numConfigs; ++t) { //间隔的构型数量
		vector<double> length_u;
		for (size_t i = 0; i < numParticles; ++i) {
			double u = 0.0;
			for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
				double diff = trajectory[t][i][i_dim] - trajectory[t - dt][i][i_dim];
				u += diff * diff;
			}
			u = sqrt(u);
			outputFile_1 << u << " ";
			length_u.emplace_back(u);
		}
		outputFile_1 << endl;
		distribution(length_u, outname_2, N_Bins);
	}
	outputFile_1.close();
}