#pragma once
const double Al_mass = 26.98;
const double Sm_mass = 150.36;
constexpr double matrix_m = 1.0;
constexpr int Dim_grid = 6;
constexpr int time_interval = 100;

inline std::vector<std::array<double, Dim>> generate_grid(const int& grid_size, const vector<double>& grid_w) { //格点的位置x, y, z

    std::vector<std::array<double, Dim>> grid;
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            for (int k = 0; k < grid_size; ++k) {
                grid.push_back({ k * grid_w[0], j * grid_w[1], i * grid_w[2] });
            }
        }
    }

    return grid;
}
inline std::vector<std::array<double, Dim>> generate_velocity(const Configuration& pos1, const Configuration& pos2) {
    std::vector<std::array<double, Dim>> velocity;
    double sum_v[3] = { 0.0 }, dr[3] = { 0.0 };
    for (size_t i = 0; i < pos1.size(); ++i) {
        if (i < Num_1) {
            for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                dr[i_Dim] = Al_mass * (pos2[i][i_Dim] - pos1[i][i_Dim]);
            }
        }
        else {
            for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                dr[i_Dim] = Sm_mass * (pos2[i][i_Dim] - pos1[i][i_Dim]);
                //cout << dr[i_Dim] << " ";
            }
            //cout << endl;
        }
//        sum_v[0] += dr[0];
//        sum_v[1] += dr[1];
//        sum_v[2] += dr[2];
        velocity.push_back({ dr[0], dr[1], dr[2] });
    }
//    sum_v[0] /= pos1.size();
//    sum_v[1] /= pos1.size();
//    sum_v[2] /= pos1.size();
//    for (size_t i = 0; i < velocity.size(); ++i) {
//        velocity[i][0] -= sum_v[0];
//        velocity[i][1] -= sum_v[1];
//        velocity[i][2] -= sum_v[2];
//    }

    return velocity;
}
inline double Distance_square(const std::vector<double>& particle1, const std::array<double, Dim>& particle2, const vector<double>& L)
{
    double r_2 = zero, diff = zero;
    for (int i = 0; i < Dim; i++) {
        diff = particle1[i] - particle2[i];
        diff = diff - L[i] * round(diff / L[i]);
        r_2 += diff * diff;
    }
    return r_2;
}
inline double Gaussian(const double& radius_square, const double& coefficient1, const double& coefficient2) {
    return exp(-radius_square * coefficient1) * coefficient2;
}
inline std::vector<std::array<double, Dim_grid>> calculate_grid_velocity(const std::vector<std::array<double, Dim>> velocity, const Configuration& posi,
    const std::vector<std::array<double, Dim>>& grid, const double& radius_square, const double& sigma, const vector<double>& L, const int& grid_size) {
    std::vector<std::array<double, Dim_grid>> grid_velocity(grid.size(), { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    double dist_square, C_ik;
    double coefficient1 = 1.0 / (2.0 * sigma * sigma);
    double coefficient2 = 1.0 / (sigma * sqrt(2.0 * pi));
    double grid_length[Dim] = { L[0] / grid_size, L[1] / grid_size, L[2] / grid_size};  // 每个格点的边长
    int i_grid;
    // 遍历每个格点

    if (smoothing == "Gaussian") {
        for (size_t i = 0; i < grid.size(); ++i) {
            for (size_t k = 0; k < Num_1; ++k) {
                dist_square = Distance_square(posi[k], grid[i], L);
                if (dist_square < radius_square) {  // 粒子落在这个格点的影响范围内
                    C_ik = Gaussian(dist_square, coefficient1, coefficient2);
                    grid_velocity[i][0] += velocity[k][0] * C_ik;
                    grid_velocity[i][1] += velocity[k][1] * C_ik;
                    grid_velocity[i][2] += velocity[k][2] * C_ik;
                }
            }
            for (size_t k = Num_1; k < posi.size(); ++k) {
                dist_square = Distance_square(posi[k], grid[i], L);
                if (dist_square < radius_square) {  // 粒子落在这个格点的影响范围内
                    C_ik = Gaussian(dist_square, coefficient1, coefficient2);
                    grid_velocity[i][3] += velocity[k][0] * C_ik;
                    grid_velocity[i][4] += velocity[k][1] * C_ik;
                    grid_velocity[i][5] += velocity[k][2] * C_ik;
                }
            }
        }
    }
    else {
        double v[Dim] = { zero };
        int i_v[Dim] = { 0 };
        for (size_t k = 0; k < Num_1; ++k) {
            for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                v[i_Dim] = posi[k][i_Dim] - L[i_Dim] * (floor(posi[k][i_Dim] / L[i_Dim]));
                i_v[i_Dim] = (int)(v[i_Dim] / grid_length[i_Dim]);
                if (i_v[i_Dim] == grid_size) {
                    i_v[i_Dim] = i_v[i_Dim] - 1;
                }
            }
            if (i_v[0] >= 0 && i_v[0] < grid_size && i_v[1] >= 0 && i_v[1] < grid_size && i_v[2] >= 0 && i_v[2] < grid_size) {
                i_grid = i_v[2] * grid_size * grid_size + i_v[1] * grid_size + i_v[0];
                grid_velocity[i_grid][0] += velocity[k][0];
                grid_velocity[i_grid][1] += velocity[k][1];
                grid_velocity[i_grid][2] += velocity[k][2];
            }
            else {
                cout << "error!" << " ix: " << i_v[0] << " iy: " << i_v[1] << " iz: " << i_v[2] << endl;
                exit(-1);
            }
        }
        for (size_t k = Num_1; k < posi.size(); ++k) {
            for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                v[i_Dim] = posi[k][i_Dim] - L[i_Dim] * (floor(posi[k][i_Dim] / L[i_Dim]));
                i_v[i_Dim] = (int)(v[i_Dim] / grid_length[i_Dim]);
                if (i_v[i_Dim] == grid_size) {
                    i_v[i_Dim] = i_v[i_Dim] - 1;
                }
            }
            if (i_v[0] >= 0 && i_v[0] < grid_size && i_v[1] >= 0 && i_v[1] < grid_size && i_v[2] >= 0 && i_v[2] < grid_size) {
                i_grid = i_v[2] * grid_size * grid_size + i_v[1] * grid_size + i_v[0];
                grid_velocity[i_grid][3] += velocity[k][0];
                grid_velocity[i_grid][4] += velocity[k][1];
                grid_velocity[i_grid][5] += velocity[k][2];
            }
            else {
                cout << "error!" << " ix: " << i_v[0] << " iy: " << i_v[1] << " iz: " << i_v[2] << endl;
                exit(-1);
            }
        }
    }

    return grid_velocity;
}

inline void write_matrix(Trajectory& trajectory, const vector<vector<double>>& L_all, const string& atom_kind) {
	unwrap_positions(trajectory, L_all);
	size_t numConfigs = trajectory.size(); // 构型数量
	size_t numParticles = trajectory[0].size(); // 每个构型的粒子数量
	int grid_size = static_cast<int>(std::floor(pow(numParticles, 1.0/3.0)));
	cout << "grid_size: " << grid_size << endl;
    std::ofstream output_Matrix(matrix_name, ios::out);

    for (size_t t = 0; t < numConfigs - time_interval; ++t) {
        cout << t << endl;
        double sigma[Dim] = { L_all[t][0] / (grid_size * matrix_m * 2.355), L_all[t][1] / (grid_size * matrix_m * 2.355), L_all[t][2] / (grid_size * matrix_m * 2.355) };
        double radius_square = 9.0 * (sigma[0] * sigma[0] + sigma[1] * sigma[1] + sigma[2] * sigma[2]);
        double sigma_length = sqrt(sigma[0] * sigma[0] + sigma[1] * sigma[1] + sigma[2] * sigma[2]);
        vector<double> grid_w = { L_all[t][0] / (grid_size), L_all[t][1] / (grid_size), L_all[t][2] / (grid_size) };
        std::vector<std::array<double, Dim>> velocity = generate_velocity(trajectory[t], trajectory[t + time_interval]);
        std::vector<std::array<double, Dim>> grid = generate_grid(grid_size, grid_w); 
        std::vector<std::array<double, Dim_grid>> grid_velocity = calculate_grid_velocity(velocity, trajectory[t], grid, radius_square, sigma_length, L_all[t], grid_size);
        
        for (int i = 0; i < grid_velocity.size(); ++i) {
            for (int j = 0; j < Dim_grid; ++j) {
                output_Matrix << grid_velocity[i][j] << " ";
            }
        }
        output_Matrix << std::endl;
        
    }
    output_Matrix.close();
}