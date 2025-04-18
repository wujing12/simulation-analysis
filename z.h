#pragma once
constexpr int N_Bins_z = 500; //
constexpr double Cut_d2_2 = 100.0;
constexpr double Cut_d2_1 = 50.0;
// 计算参数z
inline void get_z(std::vector<Vec3>& positions, const double(&L)[Dim], const std::vector<Vec3>& positions_ini, const double(&L_ini)[Dim]) {
	vector<double> zs;
	double N_Si_Si_m = 0.0, N_Si_O_m = 0.0, z;
	/*	for (int i = 0; i < positions.size(); ++i) { //减去仿射变换的偏移量
		positions[i].x = positions[i].x - ((L[0] - L_ini[0]) / L_ini[0]) * positions_ini[i].x;
		positions[i].y = positions[i].y - ((L[1] - L_ini[1]) / L_ini[1]) * positions_ini[i].y;
		positions[i].z = positions[i].z - ((L[2] - L_ini[2]) / L_ini[2]) * positions_ini[i].z;
	}*/
	string outname = dir_out + "z_atomic.txt";
	std::ofstream out_file(outname, ios::app);
	for (int i = 0; i < positions.size(); ++i) {
		if (positions[i].kind == 2) { //硅原子
			vector<double> z_1; //1邻居的距离平方
			vector<double> z_2; //2邻居的距离平方
			for (int j = 0; j < positions.size(); ++j) {
				if (positions[j].kind == 1) {
					double dist_2 = Distance_2(positions[i], positions[j], L);
					if (dist_2 < Cut_d2_1) {
						if (dist_2 <= Si_O_d2) {
							N_Si_O_m += 1.0;
						}
						z_1.emplace_back(dist_2);
					}
				}
				if ((j != i) && (positions[j].kind == 2)) {
					double dist_2 = Distance_2(positions[i], positions[j], L);
					if (dist_2 < Cut_d2_2) {
						if (dist_2 <= Si_Si_d2) { //Si邻居
							N_Si_Si_m += 1.0;
						}
						z_2.emplace_back(dist_2);
					}
				}
			}
			std::sort(z_1.begin(), z_1.end());
			std::sort(z_2.begin(), z_2.end());
			if (z_2.size() > 4 && z_1.size() > 3) {
				z = sqrt(z_2[4]) - sqrt(z_1[3]);
				zs.emplace_back(z);
			}
			else {
				cout << "error:" << z_1.size() << " " << z_2.size() << endl;
			}
			out_file << z << endl;
		}
	}
	out_file.close();

	double sum = 0.0;
	for (double value : zs) {
		sum += value;
	}
	double mean = sum / zs.size();
	double min_value = *std::min_element(zs.begin(), zs.end()) - 0.1;
	double max_value = *std::max_element(zs.begin(), zs.end()) + 0.1;

	// 计算区间大小
	double bin_size = (max_value - min_value) / N_Bins_z;

	// 创建一个计数器来统计每个区间的数量
	std::vector<int> histogram(N_Bins_z, 0);

	// 填充直方图
	for (double value : zs) {
		int bin_index = static_cast<int>((value - min_value) / bin_size);
		histogram[bin_index]++;
	}
	N_Si_Si_m = N_Si_Si_m / (Num_2);
	N_Si_O_m = N_Si_O_m / (Num_2);
	// 输出
	string outname_sta_z = dir_out + "z_sta.txt";
	std::ofstream outputFile(outname_sta_z, ios::app);
	outputFile << N_Si_Si_m << " " << N_Si_O_m << " " << mean << " " << min_value << " " << bin_size << " ";
	for (int i = 0; i < N_Bins_z; ++i) {
		outputFile << histogram[i] << " "; // 输出 g(r) 值
	}
	outputFile << endl;
	outputFile.close();
}