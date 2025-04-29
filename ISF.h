#pragma once
double S_k_max = 2.6295;//AlSm 2.67 SiO2 1.6275 CuZr: 2.75778, 2.76258, 2.76497, 2.76737
constexpr double d_k = 0.005;
inline void get_index_particle_first(size_t t, size_t i,  size_t d, size_t& index, size_t numConfigs)
{
	index = (i * numConfigs * Dim) + (t * Dim) + d;
}
inline void get_index_time_first(size_t t, size_t i, size_t d, size_t& index, size_t numParticles)
{
	index = (t * numParticles * Dim) + (i * Dim) + d;
}
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
	// 找到最接近 S_k_max 的波矢长度
	double min_diff = std::numeric_limits<double>::max();
	for (const auto& waveVec : wavevectors) {
		double diff = abs(waveVec.length - S_k_max);
		if (diff < min_diff) {
			min_diff = diff;
		}
	}

	// 只保留那些与 S_k_max 长度相差最小的波矢
	std::vector<Vec3_wave_vector> closest_wavevectors;
	for (const auto& waveVec : wavevectors) {
		if (abs(abs(waveVec.length - S_k_max) - min_diff) <= error) {
			closest_wavevectors.push_back(waveVec);
		}
	}
	wavevectors = closest_wavevectors;
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
inline vector<long double> get_SISF(Trajectory& trajectory, const vector<vector<double>>& L_traj, const std::vector<Vec3_wave_vector>& wavevectors) {
	unwrap_positions(trajectory, L_traj);
	size_t numConfigs = trajectory.size(); // 构型数量
	size_t numParticles = trajectory[0].size(); // 每个构型的粒子数量
	if (atom_kind == "com") {
		get_centroid(trajectory, numParticles);
	}

	vector<long double> F_t(numConfigs - 1, 0.0);// F_t
	double dr[Dim] = { 0.0 };
	double phi;
	if (parallel) {
		// 将三维 trajectory 展开为一维数组
		vector<double> trajectory_flat(numConfigs * numParticles * Dim, 0.0); //(t * numParticles * Dim) + (i * Dim) + d;
		// 将 trajectory 数据填充到一维数组 trajectory_flat 中
		for (size_t t = 0; t < numConfigs; ++t) {
			for (size_t p = 0; p < numParticles; ++p) {
				for (size_t d = 0; d < Dim; ++d) {
					trajectory_flat[t * numParticles * Dim + p * Dim + d] = trajectory[t][p][d];
				}
			}
		}
		Trajectory empty;
		trajectory.swap(empty); //清空内存
#pragma omp parallel for schedule(guided, 1)
		for (size_t t = 1; t < numConfigs; ++t) {  // Loop over time steps
			double F_t_t = 0.0;
			for (size_t k = 0; k < numConfigs - t; ++k) {  // Loop over initial configurations
				double F_t_k = 0.0;
				const size_t t1 = (t + k) * numParticles * Dim;
				const size_t t2 = k * numParticles * Dim;
				for (size_t p = 0; p < numParticles; ++p) {  // Loop over particles
					// Loop over time steps
					const double dx = trajectory_flat[t1 + p * Dim] - trajectory_flat[t2 + p * Dim];
					const double dy = trajectory_flat[t1 + p * Dim + 1] - trajectory_flat[t2 + p * Dim + 1];
					const double dz = trajectory_flat[t1 + p * Dim + 2] - trajectory_flat[t2 + p * Dim + 2];
					// Calculate the scattering physics contribution for each wavevector
					for (const auto& waveVec : wavevectors) {
						F_t_k += cos(dx * waveVec.x + dy * waveVec.y + dz * waveVec.z);  // Accumulate result in thread-local vector
					}
				}
				F_t_t += F_t_k / (numParticles * wavevectors.size());
			}
			F_t[t - 1] = F_t_t / (numConfigs - t);
		}
	}
	else {
		// 遍历所有时间原点 t_k
		for (size_t t = 1; t < numConfigs; ++t) { //间隔的构型数量
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
	}
	return F_t;
}