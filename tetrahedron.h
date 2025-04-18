#pragma once
inline double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline double magnitude_vec(const std::vector<double>& vec) {
	return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
// 计算向量夹角（返回弧度制角度）
inline double angleBetweenVectors(const std::vector<double>& a, const std::vector<double>& b) {
	double dot = dotProduct(a, b);
	double magA = magnitude_vec(a);
	double magB = magnitude_vec(b);
	return std::acos(dot / (magA * magB)); // 结果为弧度
}
inline void get_tetrahedron(const std::vector<Vec3>& positions, const double(&L)[Dim]) {
	string outname = dir_out + "vector_atomic.txt";
	std::ofstream out_file(outname, ios::app);
	double L_h[Dim] = { L[0] / 2.0, L[1] / 2.0, L[2] / 2.0 };
	vector<vector<Line> > Link_2; //邻居是Si
	Link_2.resize(Num_2);
	int count = 0; //Si的编号
	double dist_2, module;
	double dr[Dim] = { 0.0 };
	vector<vector<double> > vector_tetra(Num_2, std::vector<double>(Dim, zero));; //四面体的矢量
	for (int i = 0; i < positions.size(); ++i) { // 计算两种邻居
		if (positions[i].kind == 2) {
			//cout << count << endl;
			for (int j = 0; j < positions.size(); ++j) {
				if (i != j) {
					dist_2 = Distance_2(positions[i], positions[j], L);
					if ((positions[j].kind == 1) && (dist_2 < Si_O_d2)) { //氧原子						
						dr[0] = positions[j].x - positions[i].x;
						dr[1] = positions[j].y - positions[i].y;
						dr[2] = positions[j].z - positions[i].z;
						PBC(L_h, L, dr);
						for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
							vector_tetra[count][i_Dim] += dr[i_Dim];
						}
					}
					else if ((positions[j].kind == 2) && (dist_2 < Si_Si_d2)) { //硅原子
						Link_2[count].push_back(Line(j - Num_1, dist_2));
					} 
				}

			}
			count++;
		}
	}
	count = 0;
	for (int i = 0; i < positions.size(); ++i) {
		if (positions[i].kind == 2) { //硅原子
			module = magnitude_vec(vector_tetra[count]);
			out_file << module << " " << vector_tetra[count][2];
			vector<double> angle_tetra(6, -1.0); //四面体间的角度
			for (int j = 0; j < Link_2[count].size() && j < angle_tetra.size(); ++j) {
				int serial = Link_2[count][j].seq;
				double angleRad = angleBetweenVectors(vector_tetra[count], vector_tetra[serial]);
				angle_tetra[j] = angleRad * 180.0 / pi;
			}
			for (int j = 0; j < angle_tetra.size(); ++j) {
				out_file << " " << angle_tetra[j];
			}
			out_file << endl;
			count++;
		}		
	}
	out_file.close();
}
