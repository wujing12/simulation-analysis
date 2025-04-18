#pragma once
constexpr int N_Bins_angle = 200;
constexpr double toangle = 180 / pi;

// 计算三种键角的分布
inline void get_angle(const std::vector<Vec3>& positions, const double(&L)[Dim]) {
	vector<double> Si_O_Si, Si_Si_Si, O_Si_O;
	double L_h[Dim] = { L[0] / 2.0, L[1] / 2.0, L[2] / 2.0 };
	// 计算所有粒子对之间的距离
	for (int i = 0; i < positions.size(); ++i) {
		if (positions[i].kind == 2) { //硅原子
			vector<std::vector<double>> O_d_r(3);
			vector<std::vector<double>> Si_d_r(3); //Si邻居的相对位置
			for (int j = 0; j < positions.size(); ++j) {
				if (j != i) {
					double dist_2 = Distance_2(positions[i], positions[j], L);
					if (((positions[j].kind == 1) && (dist_2 < Si_O_d2)) || ((positions[j].kind == 2) && (dist_2 < Si_Si_d2))) { //O,Si邻居
						double d_r[Dim] = { positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z };
						for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
							if (d_r[i_Dim] > L_h[i_Dim]) {
								d_r[i_Dim] = d_r[i_Dim] - L[i_Dim];
							}
							else if (d_r[i_Dim] < -L_h[i_Dim]) {
								d_r[i_Dim] = d_r[i_Dim] + L[i_Dim];
							}
						}
						if (positions[j].kind == 1) {
							for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
								O_d_r[i_Dim].emplace_back(d_r[i_Dim]);
							}
						}
						else {
							for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
								Si_d_r[i_Dim].emplace_back(d_r[i_Dim]);
							}
						}
					}
				}
			}
			double dot, modu, cosTheta;
			for (int i_n = 0; i_n < O_d_r[0].size(); i_n++) {
				for (int j_n = i_n + 1; j_n < O_d_r[0].size(); j_n++) {
					double d_r1[Dim] = { O_d_r[0][i_n], O_d_r[1][i_n], O_d_r[2][i_n] };
					double d_r2[Dim] = { O_d_r[0][j_n], O_d_r[1][j_n], O_d_r[2][j_n] };
					dot = d_r1[0] * d_r2[0] + d_r1[1] * d_r2[1] + d_r1[2] * d_r2[2];
					modu = dot / sqrt((d_r1[0] * d_r1[0] + d_r1[1] * d_r1[1] + d_r1[2] * d_r1[2]) * (d_r2[0] * d_r2[0] + d_r2[1] * d_r2[1] + d_r2[2] * d_r2[2]));
					cosTheta = std::max(-1.0, std::min(1.0, modu));
					O_Si_O.emplace_back(acos(cosTheta) * toangle);
				}
			}
			for (int i_n = 0; i_n < Si_d_r[0].size(); i_n++) {
				for (int j_n = i_n + 1; j_n < Si_d_r[0].size(); j_n++) {
					double d_r1[Dim] = { Si_d_r[0][i_n], Si_d_r[1][i_n], Si_d_r[2][i_n] };
					double d_r2[Dim] = { Si_d_r[0][j_n], Si_d_r[1][j_n], Si_d_r[2][j_n] };
					dot = d_r1[0] * d_r2[0] + d_r1[1] * d_r2[1] + d_r1[2] * d_r2[2];
					modu = dot / sqrt((d_r1[0] * d_r1[0] + d_r1[1] * d_r1[1] + d_r1[2] * d_r1[2]) * (d_r2[0] * d_r2[0] + d_r2[1] * d_r2[1] + d_r2[2] * d_r2[2]));
					cosTheta = std::max(-1.0, std::min(1.0, modu));
					Si_Si_Si.emplace_back(acos(cosTheta) * toangle);
				}
			}
		}
		else { //氧原子
			vector<std::vector<double>> Si_d_r(3); //Si邻居的位置
			for (int j = 0; j < positions.size(); ++j) {
				if (j != i) {
					double dist_2 = Distance_2(positions[i], positions[j], L);
					if (((positions[j].kind == 2) && (dist_2 < Si_O_d2))) { //O邻居
						double d_r[Dim] = { positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z };
						for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
							if (d_r[i_Dim] > L_h[i_Dim]) {
								d_r[i_Dim] = d_r[i_Dim] - L[i_Dim];
							}
							else if (d_r[i_Dim] < -L_h[i_Dim]) {
								d_r[i_Dim] = d_r[i_Dim] + L[i_Dim];
							}
						}
						for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
							Si_d_r[i_Dim].emplace_back(d_r[i_Dim]);
						}
					}
				}
			}
			double dot, modu, cosTheta;
			for (int i_n = 0; i_n < Si_d_r[0].size(); i_n++) {
				for (int j_n = i_n + 1; j_n < Si_d_r[0].size(); j_n++) {
					double d_r1[Dim] = { Si_d_r[0][i_n], Si_d_r[1][i_n], Si_d_r[2][i_n] };
					double d_r2[Dim] = { Si_d_r[0][j_n], Si_d_r[1][j_n], Si_d_r[2][j_n] };
					dot = d_r1[0] * d_r2[0] + d_r1[1] * d_r2[1] + d_r1[2] * d_r2[2];
					modu = dot / sqrt((d_r1[0] * d_r1[0] + d_r1[1] * d_r1[1] + d_r1[2] * d_r1[2]) * (d_r2[0] * d_r2[0] + d_r2[1] * d_r2[1] + d_r2[2] * d_r2[2]));
					cosTheta = std::max(-1.0, std::min(1.0, modu));
					Si_O_Si.emplace_back(acos(cosTheta) * toangle);
				}
			}
		}

	}
	string outname = dir_out + "angle.txt";
	std::ofstream outputFile(outname, ios::app);
	std::cout << outname << endl;
	if (outputFile.fail()) {
		std::cout << "File_angle opening failed! " << endl;
	}
	double min_value = *std::min_element(Si_O_Si.begin(), Si_O_Si.end()) - 0.00001;
	double max_value = *std::max_element(Si_O_Si.begin(), Si_O_Si.end()) + 0.00001;
	double bin_size = (max_value - min_value) / N_Bins_angle;
	std::vector<int> histogram(N_Bins_angle, 0);
	for (double value : Si_O_Si) {
		int bin_index = static_cast<int>((value - min_value) / bin_size);
		histogram[bin_index]++;
	}
	outputFile << min_value << " " << bin_size << " ";
	cout << min_value << " " << bin_size << " ";
	for (int i = 0; i < N_Bins_angle; ++i) {
		outputFile << histogram[i] << " "; // 输出 g(r) 值
		cout << histogram[i] << " ";
	}
	cout << endl;
	double min_value1 = *std::min_element(O_Si_O.begin(), O_Si_O.end()) - 0.00001;
	double max_value1 = *std::max_element(O_Si_O.begin(), O_Si_O.end()) + 0.00001;
	double bin_size1 = (max_value1 - min_value1) / N_Bins_angle;
	std::vector<int> histogram1(N_Bins_angle, 0);
	for (double value : O_Si_O) {
		int bin_index = static_cast<int>((value - min_value1) / bin_size1);
		histogram1[bin_index]++;
	}
	outputFile << min_value1 << " " << bin_size1 << " ";
	cout << min_value1 << " " << bin_size1 << " ";
	for (int i = 0; i < N_Bins_angle; ++i) {
		outputFile << histogram1[i] << " "; // 输出 g(r) 值
		cout << histogram1[i] << " ";
	}
	cout << endl;
	double min_value2 = *std::min_element(Si_Si_Si.begin(), Si_Si_Si.end()) - 0.00001;
	double max_value2 = *std::max_element(Si_Si_Si.begin(), Si_Si_Si.end()) + 0.00001;
	double bin_size2 = (max_value2 - min_value2) / N_Bins_angle;
	std::vector<int> histogram2(N_Bins_angle, 0);
	for (double value : Si_Si_Si) {
		int bin_index = static_cast<int>((value - min_value2) / bin_size2);
		histogram2[bin_index]++;
	}
	outputFile << min_value2 << " " << bin_size2 << " ";
	cout << min_value2 << " " << bin_size2 << " ";
	for (int i = 0; i < N_Bins_angle; ++i) {
		outputFile << histogram2[i] << " "; // 输出 g(r) 值
		cout << histogram2[i] << " ";
	}
	cout << endl;
	outputFile << endl;
	outputFile.close();
}