#pragma once
inline double distance3D(const vector<double>& a, const vector<double>& b) {
    double dr2 = 0.0;
    for (int d = 0; d < Dim; ++d) {
        double diff = a[d] - b[d];
        dr2 += diff * diff;
    }
    return sqrt(dr2);
}
inline void get_chi4(Trajectory& traj, const vector<vector<double>>& L_traj, const string& atom_kind) {
    unwrap_positions(traj, L_traj);
    size_t T = traj.size();      // 总时间帧数
    size_t N = traj[0].size();   // 粒子数

    string outname_1 = dir_out + "chi4_" + atom_kind + ".txt";
    std::ofstream outputFile_1(outname_1, ios::out);
    for (size_t dt = 1; dt < T; ++dt) {
        size_t num_t0 = T - dt;
        vector<double> Q_list;

        for (size_t t0 = 0; t0 < num_t0; ++t0) {
            size_t t1 = t0 + dt;
            double Q = 0.0;
            for (size_t i = 0; i < N; ++i) {
                double dr = distance3D(traj[t1][i], traj[t0][i]);
                if (dr < Al_Al_d) {
                    Q += 1.0;
                }
            }
            Q_list.push_back(Q);
        }

        // 计算 <Q> 和 <Q^2>
        double Q_avg = accumulate(Q_list.begin(), Q_list.end(), 0.0) / Q_list.size();

        double Q2_avg = 0.0;
        for (double q : Q_list) {
            Q2_avg += q * q;
        }
        Q2_avg /= Q_list.size();

        double chi4 = (Q2_avg - Q_avg * Q_avg) / N;
        outputFile_1 << chi4 << endl;
    }
    outputFile_1.close();
}