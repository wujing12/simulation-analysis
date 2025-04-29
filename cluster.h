#pragma once
// 迭代DFS找团簇
void iterative_dfs(int start_si_idx, const vector<vector<int>>& Link_Si_Si,
    const vector<int>& kind_Si, vector<bool>& visited,
    vector<int>& current_cluster) {
    stack<int> stk;
    stk.push(start_si_idx);
    visited[start_si_idx] = true;

    while (!stk.empty()) {
        int si_idx = stk.top();
        stk.pop();
        current_cluster.push_back(si_idx);

        for (int neighbor_idx : Link_Si_Si[si_idx]) {  // Updated to access the integer index directly
            if (!visited[neighbor_idx] && kind_Si[neighbor_idx] == 1) {
                visited[neighbor_idx] = true;
                stk.push(neighbor_idx);
            }
        }
    }
}

inline void get_cluster(const std::vector<Vec3>& positions, const double(&L)[Dim], int lines) {
    if (lines > 0) { //第一帧没有团簇
        vector<vector<int>> Link_Si_O; //Si的O邻居
        vector<vector<int>> Link_Si_Si; // Si的Si邻居
        Link_Si_O.resize(Num_2);
        Link_Si_Si.resize(Num_2);
        double dist_2;

        for (int i = Num_1; i < positions.size(); ++i) { // Si的编号
            for (int j = 0; j < Num_1; ++j) { // O的编号
                dist_2 = Distance_2(positions[i], positions[j], L);
                if (dist_2 <= Si_O_d2) {
                    Link_Si_O[i - Num_1].push_back(j);  // Store the index directly
                }
            }
        }

        // 预先存储每个Si的O邻居到哈希表
        vector<unordered_set<int>> Si_O_neighbors(Num_2);
        for (int i = 0; i < Num_2; ++i) {
            for (int o : Link_Si_O[i]) {
                Si_O_neighbors[i].insert(o);
            }
        }

        // 检查Si(i)和Si(j)是否有共同O邻居
        for (int i = 0; i < Num_2; ++i) {
            for (int j = i + 1; j < Num_2; ++j) {
                // 遍历Si(j)的O邻居，检查是否在Si(i)的O邻居集合中
                for (int o : Link_Si_O[j]) {
                    if (Si_O_neighbors[i].count(o)) {
                        Link_Si_Si[i].push_back(j); // Store the index directly (now integers)
                        Link_Si_Si[j].push_back(i); // 双向添加
                        break;  // 找到一个共同O即可停止
                    }
                }
            }
        }

        // Clear Link_Si_O to free up memory after processing
        Link_Si_O.clear(); Link_Si_O.shrink_to_fit();
        Si_O_neighbors.clear(); Si_O_neighbors.shrink_to_fit();

        // File reading for Si kind
        string Si_kind_filename = dir_out + to_string(Num_1) + "_dt" + to_string(dt) + "_FL_kind.txt";
        ifstream read_kind(Si_kind_filename);
        if (!read_kind.is_open()) {
            cerr << "Error opening file " << Si_kind_filename << endl;
            exit(-1);
        }

        for (int i = 0; i < lines - 1; i++) { //从第二帧开始储存数据
            getline(read_kind, tmp);
        }
        vector<int> kind_Si(Num_2);
        for (int i = 0; i < Num_2; i++) {
            read_kind >> kind_Si[i];
        }
        read_kind.close();

        // 团簇分析
        vector<bool> visited(Num_2, false);
        vector<vector<int>> clusters; // 所有团簇
        vector<int> cluster_sizes;   // 每个团簇的尺寸

        for (int i = 0; i < Num_2; i++) {
            if (!visited[i] && kind_Si[i] == 1) {
                vector<int> current_cluster;
                iterative_dfs(i, Link_Si_Si, kind_Si, visited, current_cluster);
                clusters.push_back(current_cluster);
                cluster_sizes.push_back(current_cluster.size());
            }
        }

        // 最大团簇尺寸和团簇数目
        ofstream out_max(dir_out + "max_cluster.txt", ios::app);
        int max_size = cluster_sizes.empty() ? 0 : *max_element(cluster_sizes.begin(), cluster_sizes.end());
        out_max << max_size << " " << clusters.size() << "\n";
        out_max.close();

        // 每个团簇的尺寸
        ofstream out_sizes(dir_out + "cluster_sizes.txt", ios::app);
        for (int size : cluster_sizes) out_sizes << size << " ";
        out_sizes << "\n";
        out_sizes.close();

        // 所有团簇的粒子编号
        ofstream out_clusters(dir_out + "clusters.txt", ios::app);
        for (const auto& cluster : clusters) {
            for (int si_idx : cluster) out_clusters << si_idx << " ";
        }
        out_clusters << "\n";
        out_clusters.close();
    }

}
