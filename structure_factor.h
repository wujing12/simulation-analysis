#pragma once
constexpr double k_max = 6.0;
constexpr double k_max_2 = k_max * k_max;
constexpr bool Debye_scattering_function = false;
constexpr int N_Bins_Sq = 600;
class Structure_Factor {
public:
	int Size;
	vector<double> S_11;  // 类型1粒子的结构因子
	vector<double> S_12;  // 类型1和2粒子间的结构因子
	vector<double> S_22;  // 类型2粒子的结构因子
	vector<double> S_All; // 所有粒子的结构因子
	vector<long double> Lengths;  // 波矢长度
};
// 波矢
inline std::vector<Vec3_wave_vector> read_wavevector() {
	std::vector<Vec3_wave_vector> wavevectors;
	ifstream fw;
	string name = dir_read + "wavevector-k.dat";
	fw.open(name, ios::in);
	int N_lines = count_lines(name);
	std::cout << N_lines << endl;
	string tmp;
	int N_copy;
	double x, y, z;
	long double length;
	getline(fw, tmp);
	for (int i_l = 0; i_l < N_lines - 1; i_l++) {
		fw >> length >> N_copy;
		for (int i = 0; i < N_copy; i++) {
			fw >> x >> y >> z;
			if (length < k_max) {
				wavevectors.push_back({ x, y, z, length });
			}
		}
	}
	fw.close();
	std::cout << wavevectors.size() << endl;
	if (wavevectors.size() == 0)
	{
		cerr << "There is no wavevector!"
			<< endl;
		exit(-1);
	}
	return wavevectors;
}

inline std::vector<Vec3_wave_vector> get_wavevectors(const double(&L)[Dim]) {
	//std::cout << "Start generate wavevectors! ";
	std::vector<Vec3_wave_vector> wavevectors;
	double deltaK[3] = { 2.0 * pi / L[0], 2.0 * pi / L[1], 2.0 * pi / L[2] };
	long double length, length_2;
	//int Nk[3] = { int((k_max) / (2.0 * pi / L[0])), int((k_max) / (2.0 * pi / L[1])) , int((k_max) / (2.0 * pi / L[2])) };
	int Nk[3] = { 20, 20, 20 };
	for (int nx = -Nk[0]; nx <= Nk[0]; ++nx) { //注意如果只取自然数有一定的误差。 -Nk[0] -Nk[1] -Nk[2]
		for (int ny = -Nk[1]; ny <= Nk[1]; ++ny) {
			for (int nz = -Nk[2]; nz <= Nk[2]; ++nz) {
				length_2 = (nx * deltaK[0] * nx * deltaK[0] + ny * deltaK[1] * ny * deltaK[1] + nz * deltaK[2] * nz * deltaK[2]);
				if ((!((nx == 0) && (ny == 0) && (nz == 0))) && (length_2 < k_max_2)) {
					length = sqrt(length_2);
					wavevectors.push_back({ nx * deltaK[0], ny * deltaK[1], nz * deltaK[2], length });
				}
			}
		}
	}
	std::sort(wavevectors.begin(), wavevectors.end(), compareByid_wave_vector); //从小到大排列
	string outname = dir_out + "wave_vector.txt";
	ofstream outputFile(outname, ios::out);
	outputFile << wavevectors.size() << " " << L[0] << " " << L[1] << " " << L[2] << endl;
	for (int j = 0; j < wavevectors.size(); ++j) {
		outputFile << std::setprecision(10) << wavevectors[j].length << " " << wavevectors[j].x << " " << wavevectors[j].y << " " << wavevectors[j].z << endl;
	}
	std::cout << deltaK[0] << " " << deltaK[1] << " " << deltaK[2] << " " << Nk[0] << " " << Nk[1] << " " << Nk[2] << " " << wavevectors.size() << endl;
	if (wavevectors.size() == 0)
	{
		cerr << "There is no wavevector!"
			<< endl;
		exit(-1);
	}

	return wavevectors;
}
inline void init_structure_factor(const vector<Vec3_wave_vector>& wavevectors, Structure_Factor& Sq) {
	Sq.Size = wavevectors.size();
	Sq.Lengths.resize(Sq.Size);
	Sq.S_11.resize(Sq.Size, 0.0);
	Sq.S_12.resize(Sq.Size, 0.0);
	Sq.S_22.resize(Sq.Size, 0.0);
	Sq.S_All.resize(Sq.Size, 0.0);

	for (int i = 0; i < Sq.Size; ++i) {
		Sq.Lengths[i] = wavevectors[i].length;
	}
}
// 计算结构因子
inline void get_structure_factor(const std::vector<Vec3>& positions, const double(&L)[Dim], const std::vector<Vec3_wave_vector>& wavevectors, Structure_Factor& Sq) {
	if (Debye_scattering_function) {
		// 计算所有粒子对之间的距离
		vector<double> S_All(N_Bins_Sq, 0.0); //22
		vector<double> S_11(N_Bins_Sq, 0.0);
		vector<double> S_12(N_Bins_Sq, 0.0);
		vector<double> S_22(N_Bins_Sq, 0.0);
		vector<double> q(N_Bins_Sq, 0.0);
		double r, r2, Sq;
		for (int i = 0; i < N_Bins_Sq; ++i) {
			q[i] = k_max * (i + 0.5) / N_Bins_Sq;
		}
		double R_c = (L[0] + L[1] + L[2]) / 6.0;// sqrt((L[0] * L[0] + L[1] * L[1] + L[2] * L[2]) / (2.0 * 2.0));
		cout << "R_c: " << R_c << endl;
		double pi_Rc = pi / R_c;
		for (int j = 0; j < positions.size(); ++j) {
			for (int l = j + 1; l < positions.size(); ++l) {
				r2 = Distance_2(positions[l], positions[j], L);
				r = sqrt(r2);
				if (r < R_c) {
					for (int i = 0; i < q.size(); ++i) {
						Sq = sin(r * q[i]) * sin(pi_Rc * r) / (q[i] * pi_Rc * r2);
						S_All[i] += Sq;
						
						if ((positions[j].kind == 1) && (positions[l].kind == 1)) {
							S_11[i] += Sq;
						}
						else if ((positions[j].kind == 2) && (positions[l].kind == 2)) {
							S_22[i] += Sq;
						}
						else {
							S_12[i] += Sq;
						}
					}
				}

			}
		}
		for (int i = 0; i < N_Bins_Sq; ++i) {
			S_All[i] = 1.0 + S_All[i] * 2.0 / N;
			S_11[i] = 1.0 + S_11[i] * 2.0 / sqrt(Num_1 * Num_1);
			S_12[i] = 1.0 + S_12[i] / sqrt(Num_1 * Num_2);
			S_22[i] = 1.0 + S_22[i] * 2.0 / sqrt(Num_2 * Num_2);
		}
		string outname = dir_out + "S_q.txt";
		ofstream outputFile(outname, ios::app);
		for (int i = 0; i < N_Bins_Sq; ++i) {
			outputFile << S_All[i] << " ";
		}
		outputFile << endl;
		outputFile.close();
	}
	else {
		for (int k_idx = 0; k_idx < wavevectors.size(); ++k_idx) {
			const auto& k = wavevectors[k_idx];
			double rho_All_x = 0.0, rho_All_y = 0.0;
			double rho_1_x = 0.0, rho_1_y = 0.0;
			double rho_2_x = 0.0, rho_2_y = 0.0;

			for (int i = 0; i < N; ++i) {
				double phi = positions[i].x * k.x + positions[i].y * k.y + positions[i].z * k.z;
				double v_cos = cos(phi);
				double v_sin = sin(phi);

				rho_All_x += v_cos;
				rho_All_y += v_sin;

				if (positions[i].kind == 1) {
					rho_1_x += v_cos;
					rho_1_y += v_sin;
				}
				else if (positions[i].kind == 2) {
					rho_2_x += v_cos;
					rho_2_y += v_sin;
				}
			}

			Sq.S_All[k_idx] += (rho_All_x * rho_All_x + rho_All_y * rho_All_y) / N;
			Sq.S_11[k_idx] += (rho_1_x * rho_1_x + rho_1_y * rho_1_y) / Num_1;
			Sq.S_12[k_idx] += (rho_1_x * rho_2_x + rho_1_y * rho_2_y) / sqrt(Num_1 * Num_2);
			Sq.S_22[k_idx] += (rho_2_x * rho_2_x + rho_2_y * rho_2_y) / Num_2;
		}
	}
}
// 对相同波长的结构因子进行平均处理
void average_same_wavelengths(Structure_Factor& sf)
{
	// 1. 创建波长到索引的映射
	std::map<long double, std::vector<size_t>> wavelength_groups;
	int precision = 10;
	// 计算舍入因子
	long double factor = std::pow(10, precision);

	// 2. 按舍入后的波长分组
	for (size_t i = 0; i < sf.Lengths.size(); ++i) {
		// 舍入波长以消除浮点误差
		long double rounded_length = std::round(sf.Lengths[i] * factor) / factor;
		wavelength_groups[rounded_length].push_back(i);
	}

	// 3. 准备新的存储空间
	std::vector<double> new_S_11, new_S_12, new_S_22, new_S_All;
	std::vector<long double> new_Lengths;

	// 预留空间以避免多次分配
	new_S_11.reserve(wavelength_groups.size());
	new_S_12.reserve(wavelength_groups.size());
	new_S_22.reserve(wavelength_groups.size());
	new_S_All.reserve(wavelength_groups.size());
	new_Lengths.reserve(wavelength_groups.size());

	// 4. 对每个波长组进行平均
	for (const auto& group : wavelength_groups) {
		const auto& indices = group.second;

		// 计算平均值
		double avg_S_11 = std::accumulate(indices.begin(), indices.end(), 0.0,
			[&sf](double sum, size_t idx) { return sum + sf.S_11[idx]; }) / indices.size();

		double avg_S_12 = std::accumulate(indices.begin(), indices.end(), 0.0,
			[&sf](double sum, size_t idx) { return sum + sf.S_12[idx]; }) / indices.size();

		double avg_S_22 = std::accumulate(indices.begin(), indices.end(), 0.0,
			[&sf](double sum, size_t idx) { return sum + sf.S_22[idx]; }) / indices.size();

		double avg_S_All = std::accumulate(indices.begin(), indices.end(), 0.0,
			[&sf](double sum, size_t idx) { return sum + sf.S_All[idx]; }) / indices.size();

		// 存储结果
		new_Lengths.push_back(group.first); // 使用舍入后的波长
		new_S_11.push_back(avg_S_11);
		new_S_12.push_back(avg_S_12);
		new_S_22.push_back(avg_S_22);
		new_S_All.push_back(avg_S_All);
	}

	// 5. 更新结构因子数据
	sf.Size = wavelength_groups.size();
	sf.S_11 = std::move(new_S_11);
	sf.S_12 = std::move(new_S_12);
	sf.S_22 = std::move(new_S_22);
	sf.S_All = std::move(new_S_All);
	sf.Lengths = std::move(new_Lengths);
}
inline void write_structure_factor(Structure_Factor& Sq, const int& number_sample) {
	if (!Debye_scattering_function) {
		average_same_wavelengths(Sq);
		string outname = dir_out + "sq.txt";
		ofstream outputFile(outname, ios::out);
		for (int i = 0; i < Sq.Size; ++i) {
			Sq.S_11[i] /= (number_sample);
			Sq.S_12[i] /= (number_sample);
			Sq.S_22[i] /= (number_sample);
			Sq.S_All[i] /= (number_sample);
			outputFile << std::setprecision(10) << Sq.Lengths[i] << " " << Sq.S_11[i] << " " << Sq.S_12[i] << " " << Sq.S_22[i] << " " << Sq.S_All[i]  << endl;
		}
		outputFile.close();
	}
}
