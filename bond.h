#pragma once
inline void get_bond(const std::vector<Vec3>& positions, const double(&L)[Dim], const bool& initial) {
	string Si_O_bond = dir_out + "Si_O_bond.txt";
	string length_Si_O_sta = dir_out + "length_Si_O_sta.txt"; //储存非仿射位移的统计
	string dr_Si_O_sta = dir_out + "dr_Si_O_sta.txt";
	string dr2_Si_O_sta = dir_out + "dr2_Si_O_sta.txt";
	string force_Si_O_sta = dir_out + "force_Si_O_sta.txt";
	string force_r_sta = dir_out + "force_r_sta.txt";
	string projection_Si_O_sta = dir_out + "projection_Si_O_sta.txt";
	string projection_Si_Si_sta = dir_out + "projection_Si_Si_sta.txt";
	if (initial) {
		reset_file(Si_O_bond);
		reset_file(length_Si_O_sta);
		reset_file(dr_Si_O_sta);
		reset_file(dr2_Si_O_sta);
		reset_file(force_Si_O_sta);
		reset_file(force_r_sta);
		reset_file(projection_Si_O_sta);
		reset_file(projection_Si_Si_sta);
	}
	
	std::ofstream out_file(Si_O_bond, ios::app);
	double L_h[Dim] = { L[0] / 2.0, L[1] / 2.0, L[2] / 2.0 };
	int i_r, N_Bins = 200, serial; 
	const int N_bond = 6;
	double dist_2, r, force_stretch, dr_2;
	vector<double> length_Si_O;
	vector<double> dr_Si_O;
	vector<double> force_Si_O;
	vector<double> force_r;
	vector<double> dr2_Si_O;
	vector<double> projection_Si_O;
	vector<double> projection_Si_Si;
	for (int i = 0; i < positions.size(); ++i) { // 计算两种邻居
		if (positions[i].kind == 2) { //硅原子
			vector<Line> Link_1; //邻居是O
			for (int j = 0; j < positions.size(); ++j) {
				if (positions[j].kind == 1) {
					dist_2 = Distance_2(positions[i], positions[j], L);
					if (dist_2 < Si_O_d2) { //氧原子	
						r = sqrt(dist_2);
						Link_1.push_back(Line(j, r));
					}
				}
				else if (j > i) { //硅邻居	
					dist_2 = Distance_2(positions[i], positions[j], L);
					if (dist_2 < Si_Si_d2) { 
						r = sqrt(dist_2);
						double dr[Dim] = { (positions[i].x - positions[j].x), (positions[i].y - positions[j].y), (positions[i].z - positions[j].z) };
						PBC(L_h, L, dr);
						projection_Si_Si.emplace_back(dr[2] / r);
					}
				}

			}
			vector<double> length_Si_O_bond(N_bond, zero); //
			vector<double> dr_Si_O_bond(N_bond, zero);
			vector<double> force_Si_O_bond(N_bond, zero);
			for (int j = 0; j < Link_1.size() && j < N_bond; ++j) {
				serial = Link_1[j].seq;
				r = Link_1[j].r;
				double dr[Dim] = { (positions[i].x - positions[serial].x), (positions[i].y - positions[serial].y), (positions[i].z - positions[serial].z) };
				PBC(L_h, L, dr);
				i_r = (int)(r / interval);
				dr_2 = dr[2] * dr[2];
				force_stretch = forces_Si_O[i_r] * dr_2 / r; //拉伸方向的力
				length_Si_O_bond[j] = r;
				dr_Si_O_bond[j] = dr_2;
				force_Si_O_bond[j] = force_stretch;
				length_Si_O.emplace_back(r);
				dr_Si_O.emplace_back(dr[2]);
				dr2_Si_O.emplace_back(dr_2);
				force_Si_O.emplace_back(force_stretch);
				force_r.emplace_back(forces_Si_O[i_r]);
				projection_Si_O.emplace_back(dr[2]/r);
			}
			for (int j = 0; j < N_bond; ++j) {
				out_file << length_Si_O_bond[j] << " ";
			}
			for (int j = 0; j < N_bond; ++j) {
				out_file << dr_Si_O_bond[j] << " ";
			}
			for (int j = 0; j < N_bond; ++j) {
				out_file << force_Si_O_bond[j] << " ";
			}
			out_file << endl;
		}
	}

	out_file.close();
	distribution(length_Si_O, length_Si_O_sta, N_Bins);
	distribution(dr_Si_O, dr_Si_O_sta, N_Bins);
	distribution(dr2_Si_O, dr2_Si_O_sta, N_Bins);
	distribution(force_Si_O, force_Si_O_sta, N_Bins);
	distribution(force_r, force_r_sta, N_Bins);
	distribution(projection_Si_O, projection_Si_O_sta, N_Bins);
	distribution(projection_Si_Si, projection_Si_Si_sta, N_Bins);
}