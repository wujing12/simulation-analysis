#pragma once

inline void get_matrix_water(const vector<Vec3>& Positions, const float L, const int GRID_SIZE) {
    std::vector<std::vector<std::vector<int>>> grid_count(GRID_SIZE, std::vector<std::vector<int>>(GRID_SIZE, std::vector<int>(GRID_SIZE, 0)));

    double grid_length = L / GRID_SIZE;  // 每个格点的边长

    //统计每个格点的粒子数
    for (int i = 0; i < Positions.size(); i += 3) { //统计所有的氧原子
        int ix = (int)(Positions[i].x / grid_length);
        int iy = (int)(Positions[i].y / grid_length);
        int iz = (int)(Positions[i].z / grid_length);

        // 确保索引在有效范围内
        if (ix >= 0 && ix < GRID_SIZE && iy >= 0 && iy < GRID_SIZE && iz >= 0 && iz < GRID_SIZE) {
            grid_count[ix][iy][iz]++;
        }
        else {
            printf("Error: Particle %d is outside the box! (%.2f, %.2f, %.2f)\n", i, Positions[i].x, Positions[i].y, Positions[i].z);
            exit(-1);
        }
    }
    std::ofstream output_Matrix(matrix_name, ios::app);
    output_Matrix << L * L * L << " ";
    for (int k = 0; k < GRID_SIZE; k++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int i = 0; i < GRID_SIZE; i++) {
                output_Matrix << grid_count[i][j][k] << " ";
            }
        }
    }
    output_Matrix << endl;
    output_Matrix.close();
}