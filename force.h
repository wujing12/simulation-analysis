#pragma once
// 定义常量和参数
std::map<std::string, double> q;
std::map<std::string, double> A, B, C, D;
const double r_cut_Coul = 10.0;
const double r_cut_LJ = 8.0;
const double gamma_val = 0.2;
const double r_min = 0.0;
const double r_max = 2.2;
const double interval = 0.001;
const int num_elements = 10000;
std::vector<double> forces_Si_O(num_elements);
std::vector<double> forces_Si_Si(num_elements);
std::vector<double> forces_O_O(num_elements);
std::vector<double> forces_O_Si(num_elements);
// 势能函数及其导数

double U_Coul(double r, const std::string& pair) {
    if (r < r_cut_Coul) {
        return q[pair.substr(0, 2)] * q[pair.substr(3, 2)] * (1.0 / r - 1.0 / r_cut_Coul + (r - r_cut_Coul) / (r_cut_Coul * r_cut_Coul));
    }
    else {
        return 0.0;
    }
}

double G_Coul(double r, const std::string& pair) {
    return exp(-gamma_val * gamma_val / ((r - r_cut_Coul) * (r - r_cut_Coul)));
}

double U_LJ(double r, const std::string& pair) {
    if (r < r_cut_LJ) {
        return A[pair] * exp(-B[pair] * r) - C[pair] / pow(r, 6) + D[pair] / pow(r, 24);
    }
    else {
        return 0.0;
    }
}

double G_LJ(double r, const std::string& pair) {
    return exp(-gamma_val * gamma_val / ((r - r_cut_LJ) * (r - r_cut_LJ)));
}

double U(double r, const std::string& pair) {
    return U_Coul(r, pair) * G_Coul(r, pair) + U_LJ(r, pair) * G_LJ(r, pair);
}

double dU_dr(double r, const std::string& pair) {
    // 对 U(r) 关于 r 求导数
    double dU_Coul_dr = 0.0;
    if (r < r_cut_Coul) {
        dU_Coul_dr = -q[pair.substr(0, 2)] * q[pair.substr(3, 2)] * (1.0 / (r * r) - 1.0 / (r_cut_Coul * r_cut_Coul));
    }
    double dG_Coul_dr = G_Coul(r, pair) * (2 * gamma_val * gamma_val / ((r - r_cut_Coul) * (r - r_cut_Coul) * (r - r_cut_Coul)));

    double dU_LJ_dr = 0.0;
    if (r < r_cut_LJ) {
        dU_LJ_dr = -A[pair] * B[pair] * exp(-B[pair] * r) + 6 * C[pair] / pow(r, 7) - 24 * D[pair] / pow(r, 25);
    }
    double dG_LJ_dr = G_LJ(r, pair) * (2 * gamma_val * gamma_val / ((r - r_cut_LJ) * (r - r_cut_LJ) * (r - r_cut_LJ)));

    return dU_Coul_dr * G_Coul(r, pair) + U_Coul(r, pair) * dG_Coul_dr + dU_LJ_dr * G_LJ(r, pair) + U_LJ(r, pair) * dG_LJ_dr;
}

double force(double r, const std::string& pair) {
    // 力的形式：F(r) = -dU/dr
    return -dU_dr(r, pair);
}
//初始化力
inline void initialization_force() {
    std::string pair_Si_O = "Si-O_";
    std::string pair_Si_Si = "Si-Si";
    std::string pair_O_O = "O_-O_";
    std::string pair_O_Si = "O_-Si";
    q["Si"] = 6.737574454;
    q["O_"] = -3.368787227;

    A["Si-Si"] = 2797.9;
    A["Si-O_"] = 23107.8;
    A["O_-Si"] = 23107.8;
    A["O_-O_"] = 1120.5;

    B["Si-Si"] = 4.407;
    B["Si-O_"] = 5.098;
    B["O_-Si"] = 5.098;
    B["O_-O_"] = 2.893;

    C["Si-Si"] = 0.0;
    C["Si-O_"] = 139.7;
    C["O_-Si"] = 139.7;
    C["O_-O_"] = 26.1;

    D["Si-Si"] = 3423204.0;
    D["Si-O_"] = 66.0;
    D["O_-Si"] = 66.0;
    D["O_-O_"] = 16800.0;
    cout << pair_Si_O.substr(0, 2) << " " << pair_Si_O.substr(3, 2) << endl;
    cout << pair_Si_Si.substr(0, 2) << " " << pair_Si_Si.substr(3, 2) << endl;
    cout << pair_O_O.substr(0, 2) << " " << pair_O_O.substr(3, 2) << endl;
    cout << pair_O_Si.substr(0, 2) << " " << pair_O_Si.substr(3, 2) << endl;
    for (int i = 0; i < num_elements; ++i) {
        double r = r_min + (i + 0.5) * interval;
        forces_Si_O[i] = force(r, pair_Si_O);
        forces_Si_Si[i] = force(r, pair_Si_Si);
        forces_O_O[i] = force(r, pair_O_O);
        forces_O_Si[i] = force(r, pair_O_Si);
    }
    string outname = dir_out + "force.txt";
    std::ofstream outputFile(outname, ios::out);
    for (int i = 1000; i < num_elements; i += 1) {
        double r = r_min + (i + 0.5) * interval;
        outputFile << r << " " << forces_Si_O[i] << " " << forces_Si_Si[i] << " " << forces_O_O[i] << " " << forces_O_Si[i] << endl;
    }
    outputFile.close();
}
class Line //邻居结构
{
public:
    Line(const int& seq, const double& r) : seq(seq), r(r) {}
    Line() : seq(0), r(2.0) {}
    int seq; //唯一编号
    double r; //距离
    // 重载操作符
    bool operator<(const Line& rhs) {
        return (*this).seq < rhs.seq;
    };
    bool operator>(const Line& rhs) {
        return (*this).seq > rhs.seq;
    };
    bool operator==(const Line& rhs) {
        return (*this).seq == rhs.seq;
    }
};

// 计算力链网络
inline void compute_force_chain(const std::vector<Vec3>& positions, const double(&L)[Dim]) {
    string force_chain_txt = "force_atomic.txt"; //储存三个维度的位力
    string force_stat_txt = "stress_sta.txt"; //储存三维合力
    string position_txt = "position_atomic.txt"; //储存盒子信息、三维位置
    string filename_force_chain = dir_out + force_chain_txt;
    string filename_force_stat = dir_out + force_stat_txt;
    string filename_position = dir_out + position_txt;
    std::ofstream out_force_chain(filename_force_chain, ios::app);
    std::ofstream out_force_stat(filename_force_stat, ios::app);
    std::ofstream out_position(filename_position, ios::app);
    double r, dist_2, r_F[Dim] = {zero};
    int i_r, i_link;
    vector<vector<Line> > Link_1; //邻居是O
    vector<vector<Line> > Link_2; //邻居是Si
    Link_1.resize(positions.size());
    Link_2.resize(positions.size());
    for (int i = 0; i < positions.size(); ++i) { // 计算两种邻居
        for (int j = i + 1; j < positions.size(); ++j) {
            dist_2 = Distance_2(positions[i], positions[j], L);
            r = sqrt(dist_2);
            if (r <= r_cut_Coul) {
                if ((positions[i].kind == 1) && (positions[j].kind == 1)) {
                    Link_1[i].push_back(Line(j, r));
                    Link_1[j].push_back(Line(i, r));
                }
                else if ((positions[i].kind == 1) && (positions[j].kind == 2)) {
                    Link_2[i].push_back(Line(j, r));
                    Link_1[j].push_back(Line(i, r));
                }
                else if ((positions[i].kind == 2) && (positions[j].kind == 1)) {
                    Link_1[i].push_back(Line(j, r));
                    Link_2[j].push_back(Line(i, r));
                }
                else {
                    Link_2[i].push_back(Line(j, r));
                    Link_2[j].push_back(Line(i, r));
                }
            }
        }
    }
    //out_position << L[0] << " " << L[1] << " " << L[2] << endl;
    double L_h[Dim] = { L[0] / 2.0, L[1] / 2.0, L[2] / 2.0 };
    for (int i = 0; i < positions.size(); ++i) {
        double force[Dim] = { zero }; //三个方向的合力
        if (positions[i].kind == 1) {
            for (size_t j = 0; j < Link_1[i].size(); j++) { //O-O
                i_link = Link_1[i][j].seq;
                r = Link_1[i][j].r;
                double dr[Dim] = { (positions[i].x - positions[i_link].x), (positions[i].y - positions[i_link].y), (positions[i].z - positions[i_link].z) };
                PBC(L_h, L, dr);
                i_r = (int)(r / interval);
                for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                    force[i_Dim] += forces_O_O[i_r] * dr[i_Dim] * dr[i_Dim] / r;
                }
            }
            for (size_t j = 0; j < Link_2[i].size(); j++) { //O-Si
                i_link = Link_2[i][j].seq;
                r = Link_2[i][j].r;
                double dr[Dim] = { (positions[i].x - positions[i_link].x), (positions[i].y - positions[i_link].y), (positions[i].z - positions[i_link].z) };
                PBC(L_h, L, dr);
                i_r = (int)(r / interval);
                for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                    force[i_Dim] += forces_O_Si[i_r] * dr[i_Dim] * dr[i_Dim] / r;
                }
            }
        } else {
            for (size_t j = 0; j < Link_1[i].size(); j++) { //Si-O
                i_link = Link_1[i][j].seq;
                r = Link_1[i][j].r;
                double dr[Dim] = { (positions[i].x - positions[i_link].x), (positions[i].y - positions[i_link].y), (positions[i].z - positions[i_link].z) };
                PBC(L_h, L, dr);
                i_r = (int)(r / interval);
                for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                    force[i_Dim] += forces_Si_O[i_r] * dr[i_Dim] * dr[i_Dim] / r;
                }
            }
            for (size_t j = 0; j < Link_2[i].size(); j++) { //Si-Si
                i_link = Link_2[i][j].seq;
                r = Link_2[i][j].r;
                double dr[Dim] = { (positions[i].x - positions[i_link].x), (positions[i].y - positions[i_link].y), (positions[i].z - positions[i_link].z) };
                PBC(L_h, L, dr);
                i_r = (int)(r / interval);
                for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
                    force[i_Dim] += forces_Si_Si[i_r] * dr[i_Dim] * dr[i_Dim] / r;
                }
            }
        }
        r_F[0] += force[0];
        r_F[1] += force[1];
        r_F[2] += force[2];
        out_position << positions[i].x << " " << positions[i].y << " " << positions[i].z << endl;
        out_force_chain << force[0]<< " " << force[1] << " " << force[2] << endl;
    }
    out_force_stat << r_F[0] << " "
        << r_F[1] << " "
        << r_F[2] << " "
        << L[0] << " "
        << L[1] << " "
        << L[2] << endl;
    out_force_chain.close();
    out_force_stat.close();
    out_position.close();
}