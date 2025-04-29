#pragma once
constexpr double maxDistance = 40.0; //埃 L/2
constexpr int N_Bins_r = 4000;//区间间隔为0.01	
constexpr double binWidth = maxDistance / N_Bins_r;

// 计算径向分布函数 g(r)
inline void get_rdf(const std::vector<Vec3>& positions, const double(&L)[Dim], vector<double>& rdf_11, vector<double>& rdf_12, vector<double>& rdf_22, vector<double>& rdf_All) {
	vector<double> i_rdf_11(N_Bins_r, 0.0); //11
	vector<double> i_rdf_12(N_Bins_r, 0.0); //12 21
	vector<double> i_rdf_22(N_Bins_r, 0.0); //22
	vector<double> i_rdf_All(N_Bins_r, 0.0); //22
	//Sect RDF(zero, maxDistance, N_Bins_r);
	double V = L[0] * L[1] * L[2];

	double density = positions.size() / V;

	// 计算所有粒子对之间的距离
	for (int i = 0; i < positions.size(); ++i) {
		for (int j = i + 1; j < positions.size(); ++j) {
			double distance = sqrt(Distance_2(positions[i], positions[j], L));
			//std::cout << distance << endl;
			// 将距离放入适当的bin
			if (distance < maxDistance) {
				int binIndex = static_cast<int>(distance / binWidth);
				i_rdf_All[binIndex] = i_rdf_All[binIndex] + 2.0;
				if ((positions[i].kind == 1) && (positions[j].kind == 1)) {
					i_rdf_11[binIndex] = i_rdf_11[binIndex] + 2.0;
				}
				else if (((positions[i].kind == 1) && (positions[j].kind == 2)) || ((positions[i].kind == 2) && (positions[j].kind == 1))) {
					i_rdf_12[binIndex] = i_rdf_12[binIndex] + 2.0;
				}
				else if ((positions[i].kind == 2) && (positions[j].kind == 2)) {
					i_rdf_22[binIndex] = i_rdf_22[binIndex] + 2.0;
				}
				else {
					std::cout << positions[i].kind << ", error " << positions[j].kind << endl;
				}
			}
		}
	}
	double p_1 = (double)Num_1 / (double)positions.size();
	double p_2 = (double)Num_2 / (double)positions.size();
	cout << "p_1: " << p_1 << ". p_2: " << p_2 << endl;
	// 归一化 RDF
	for (int i = 0; i < N_Bins_r; ++i) {
		double r1 = i * binWidth;
		double r2 = (i + 1) * binWidth;
		double shellVolume = (4.0 / 3.0) * pi * (pow(r2, 3) - pow(r1, 3)); // 壳层体积
		double idealGasCount = shellVolume * density; // 理想气体中该距离范围内的粒子对数
		rdf_11[i] += i_rdf_11[i] / (p_1 * p_1 * positions.size() * idealGasCount); // 归一化为 g(r)
		rdf_12[i] += i_rdf_12[i] / (2.0 * p_1 * p_2 * positions.size() * idealGasCount);
		rdf_22[i] += i_rdf_22[i] / (p_2 * p_2 * positions.size() * idealGasCount);
		rdf_All[i] += i_rdf_All[i] / (positions.size() * idealGasCount);
	}
}

// 输出 RDF
inline void write_rdf(const vector<double>& rdf_11, const vector<double>& rdf_12, const vector<double>& rdf_22, const vector<double>& rdf_All, const int& number_sample) {
	string outname = dir_out + "rdf.txt";
	std::ofstream outputFile(outname, ios::out);
	for (int i = 0; i < N_Bins_r; ++i) {
		outputFile << (i + 0.5) * maxDistance / N_Bins_r << " " << rdf_11[i] / number_sample << " " << rdf_12[i] / number_sample << " " << rdf_22[i] / number_sample << " " << rdf_All[i] / number_sample << endl;
	}
	outputFile.close();
}

// 计算O的配位数
inline void get_O_Si(const std::vector<Vec3>& positions, const double(&L)[Dim]) {
	double N_O_Si = 0.0, N_O_Si_2 = 0.0;

	// 计算所有粒子对之间的距离
	for (int i = 0; i < positions.size(); ++i) {
		if (positions[i].kind == 1) { //氧原子
			double N_Si = 0.0;
			for (int j = 0; j < positions.size(); ++j) {
				if (positions[j].kind == 2) {//Si邻居
					double dist_2 = Distance_2(positions[i], positions[j], L);
					if (dist_2 <= Si_O_d2) {
						N_Si += 1.0;
					}
				}
			}
			N_O_Si_2 += N_Si * N_Si;
			N_O_Si += N_Si;
		}
	}
	N_O_Si = N_O_Si / (Num_1);
	N_O_Si_2 = N_O_Si_2 / (Num_2);
	// 输出
	string outname = dir_out + "O_Si.txt";
	std::ofstream outputFile(outname, ios::app);
	outputFile << N_O_Si << " " << sqrt(N_O_Si_2 - N_O_Si * N_O_Si) << " " << N_O_Si_2 << endl;
	outputFile.close();
}