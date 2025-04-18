#pragma once
inline void get_centroid(Trajectory& trajectory, size_t& numParticles) {
    // Precompute normalized masses
    double total_mass = 0.0;
    vector<double> norm_masses(N_atom_mol);
    for (size_t serial = 0; serial < N_atom_mol; ++serial) {
        total_mass += atom_masses_Toluene[serial];
    }
    for (size_t serial = 0; serial < N_atom_mol; ++serial) {
        norm_masses[serial] = atom_masses_Toluene[serial] / total_mass;
    }

    // Convert trajectory to COM in-place
    for (size_t t = 0; t < trajectory.size(); ++t) {
        for (size_t j = 0; j < N_mol; ++j) {
            double com[3] = { 0.0, 0.0, 0.0 };
            for (size_t serial = 0; serial < N_atom_mol; ++serial) {
                com[0] += norm_masses[serial] * trajectory[t][j * N_atom_mol + serial][0];
                com[1] += norm_masses[serial] * trajectory[t][j * N_atom_mol + serial][1];
                com[2] += norm_masses[serial] * trajectory[t][j * N_atom_mol + serial][2];
            }
            // Store COM in first atom's position
            trajectory[t][j][0] = com[0];
            trajectory[t][j][1] = com[1];
            trajectory[t][j][2] = com[2];
        }
    }
    numParticles = N_mol;
}
inline void unwrap_positions(Trajectory& trajectory, const vector<vector<double>>& L_traj) {
    size_t T = trajectory.size();      // ��������
    size_t N = trajectory[0].size();  // ������
    
    for (size_t i = 0; i < N; ++i) {
        Configuration dr(T, Position(3, 0.0)); //�����ӵ�dr
        for (size_t t = 1; t < T; ++t) {
            for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
                double L_h = L_traj[t][i_dim] / 2.0;
                double delta = trajectory[t][i][i_dim] - trajectory[t - 1][i][i_dim];
                // ������߽��λ��
                if (delta > L_h) { //�Ӻ��ӵ���߳�ȥ
                    delta -= L_traj[t][i_dim];
                }
                else if (delta < -L_h) { //�Ӻ��ӵ��ұ߳�ȥ
                    delta += L_traj[t][i_dim];
                }
                // ��������
                dr[t][i_dim] = delta;
            }
        }

        for (size_t t = 1; t < T; ++t) {
            for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
                // �ۼ�����
                trajectory[t][i][i_dim] = trajectory[t - 1][i][i_dim] + dr[t][i_dim];
            }
        }
    }
}

// �������λ�� (MSD)
inline void get_MSD(Trajectory& trajectory, const vector<vector<double>>& L_traj, const string& atom_kind) {
    unwrap_positions(trajectory, L_traj);
    size_t numConfigs = trajectory.size(); // ��������
    size_t numParticles = trajectory[0].size(); // ÿ�����͵���������

    if (atom_kind == "com") {
        get_centroid(trajectory, numParticles);
    }
    vector<double> msd(numConfigs - 2, 0.0);// MSD �ۼ�
    // ��������ʱ��ԭ�� t_k
    double dr[Dim] = { 0.0 };
    for (size_t t = 1; t < numConfigs - 1; ++t) { //����Ĺ�������
        for (size_t k = 0; k < numConfigs - t; ++k) { //��ʼ��k������
            // �������������
            for (size_t i = 0; i < numParticles; ++i) {
                dr[0] = trajectory[k + t][i][0] - trajectory[k][i][0];
                dr[1] = trajectory[k + t][i][1] - trajectory[k][i][1];
                dr[2] = trajectory[k + t][i][2] - trajectory[k][i][2];
                // �ۼ�����λ��ƽ��
                msd[t - 1] += dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            }
        }
        // ƽ������ʱ��ԭ�������
        msd[t - 1] /= (numParticles * (numConfigs - t));
    }
    string outname = dir_out + "MSD_" + atom_kind + ".txt";
    std::ofstream outputFile(outname, ios::out);
    for (int j = 0; j < msd.size(); ++j) {
        outputFile << msd[j] << " ";
    }
    outputFile << endl;
    outputFile.close();
}
inline void get_momentum(Trajectory& trajectory, const vector<vector<double>>& L_traj, const string& atom_kind) {
    unwrap_positions(trajectory, L_traj);
    size_t numConfigs = trajectory.size(); // ��������
    size_t numParticles = trajectory[0].size(); // ÿ�����͵���������

    string outname = dir_out + atom_kind + ".txt";
    std::ofstream outputFile(outname, ios::out);
    for (size_t t = 1; t < numConfigs; ++t) { //����Ĺ�������
        double dr[Dim] = { 0.0 };
        for (size_t i = 0; i < numParticles; ++i) {
            for (size_t i_dim = 0; i_dim < Dim; ++i_dim) {
                dr[i_dim] += trajectory[t][i][i_dim] - trajectory[t - 1][i][i_dim];
            }
        }
        outputFile << (dr[0] / numParticles) << " " << (dr[1] / numParticles) << " " << (dr[2] / numParticles) << endl;
    }
    outputFile.close();
}
