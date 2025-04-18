#pragma once
const double T_b = 2.728;
const double T_k = -0.0001027;
double S_k_max = 2.6295;//AlSm 2.67 SiO2 1.6275 CuZr: 2.75778, 2.76258, 2.76497, 2.76737
constexpr double d_k = 0.005;

inline std::vector<Vec3_wave_vector> max_wavevectors(const double(&L)[Dim]) {
	cout << "S_k_max: " << S_k_max << endl;
	//std::cout << "Start generate wavevectors! ";
	std::vector<Vec3_wave_vector> wavevectors;
	long double length;
	double deltaK[3] = { 2.0 * pi / L[0], 2.0 * pi / L[1], 2.0 * pi / L[2] };
	int Nk[3] = { int((S_k_max + d_k) / (2.0 * pi / L[0])) + 1, int((S_k_max + d_k) / (2.0 * pi / L[1])) + 1, int((S_k_max + d_k) / (2.0 * pi / L[2])) + 1 };
	for (int nx = -Nk[0]; nx <= Nk[0]; ++nx) {
		for (int ny = -Nk[1]; ny <= Nk[1]; ++ny) {
			for (int nz = -Nk[2]; nz <= Nk[2]; ++nz) {
				length = sqrt(nx * deltaK[0] * nx * deltaK[0] + ny * deltaK[1] * ny * deltaK[1] + nz * deltaK[2] * nz * deltaK[2]);
				if (abs(length - S_k_max) <= d_k) {
					wavevectors.push_back({ nx * deltaK[0], ny * deltaK[1], nz * deltaK[2], length });
				}
			}
		}
	}
	std::sort(wavevectors.begin(), wavevectors.end(), compareByid_wave_vector);
	string outname = dir_out + "wavevectors_max.txt";
	ofstream outputFile(outname, ios::out);
	outputFile << wavevectors.size() << " " << L[0] << " " << L[1] << " " << L[2] << endl;
	for (int j = 0; j < wavevectors.size(); ++j) {
		outputFile << std::setprecision(10) << wavevectors[j].length << " " << wavevectors[j].x << " " << wavevectors[j].y << " " << wavevectors[j].z << endl;
	}

	outputFile.close();
	std::cout << deltaK[0] << " " << deltaK[1] << " " << deltaK[2] << " " << wavevectors.size() << endl;
	if (wavevectors.size() == 0)
	{
		cerr << "There is no wavevector!"
			<< endl;
		exit(-1);
	}
	return wavevectors;
}
// 计算中间散射函数
inline void get_ISF(const std::vector<Vec3>& Positions, vector<double>& F_t_x, vector<double>& F_t_y, const std::vector<Vec3_wave_vector>& wavevectors) {

	double phi, value_x, value_y;
	for (int i = 0; i < N; ++i) {
		for (const auto& k : wavevectors) {
			phi = Positions[i].x * k.x + Positions[i].y * k.y + Positions[i].z * k.z;
			value_x = cos(phi);
			//value_y = sin(phi);
			F_t_x[2] += value_x;
			//F_t_y[2] += value_y;
			if (Positions[i].kind == 1) {
				F_t_x[0] += value_x;
				//F_t_y[0] += value_y;
			}
			else if (Positions[i].kind == 2) {
				F_t_x[1] += value_x;
				//F_t_y[1] += value_y;
			}
		}
	}
	for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
		F_t_x[i_Dim] = F_t_x[i_Dim] / (wavevectors.size());
		//F_t_y[i_Dim] = F_t_y[i_Dim] / (wavevectors.size());
	}
}
// 计算自中间散射函数
inline vector<double> get_SISF(Trajectory& trajectory, const vector<vector<double>>& L_traj, const std::vector<Vec3_wave_vector>& wavevectors) {
	unwrap_positions(trajectory, L_traj);
	size_t numConfigs = trajectory.size(); // 构型数量
	size_t numParticles = trajectory[0].size(); // 每个构型的粒子数量
	vector<double> F_t(numConfigs - 2, 0.0);// F_t
	// 遍历所有时间原点 t_k
	double dr[Dim] = { 0.0 };
	double phi;
	for (size_t t = 1; t < numConfigs - 1; ++t) { //间隔的构型数量
		cout << t << endl;
		for (size_t k = 0; k < numConfigs - t; ++k) { //初始第k个构型
			// 对所有粒子求和
			for (size_t i = 0; i < numParticles; ++i) {
				dr[0] = trajectory[k + t][i][0] - trajectory[k][i][0];
				dr[1] = trajectory[k + t][i][1] - trajectory[k][i][1];
				dr[2] = trajectory[k + t][i][2] - trajectory[k][i][2];
				for (const auto& waveVec : wavevectors) {
					phi = dr[0] * waveVec.x + dr[1] * waveVec.y + dr[2] * waveVec.z;
					F_t[t - 1] += cos(phi);
				}
			}
			// 平均所有时间原点和粒子
			
		}
		F_t[t - 1] /= (numParticles * (numConfigs - t) * wavevectors.size());
	}
	return F_t;
}