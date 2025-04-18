#pragma once
//#define gromacs_traj
#define gen_cluster
//#define gen_tetrahedron
//#define gen_non_affine_D
//#define gen_sq
//#define gen_rdf
//#define gen_MSD
//#define gen_SISF
//#define gen_position
//#define gen_trajectory
//#define gen_z
//#define gen_O_Si
//#define gen_force_chain
//#define gen_angle
//#define gen_matrix
//#define gen_bond
//#define gen_eigen
//#define gen_D_NA
//#define gen_p_hop
//#define gen_momentum
//#define gen_p_u
//#define gen_chi4
//#define gen_matrix_water
//#define read_wavevectors
//#define gen_ISF

int N = 10000;// ������
int N_mol = N;// ������
int N_atom_mol = 1;// ÿ�����ӵ�ԭ����
double mass_1 = 63.55; //����1������
double mass_2 = 91.22; //����2������
int ini_config = 0;
int total_config = 1000; 
int d_S = 1;
constexpr int i_d_S = 0;
constexpr double zero = 0.0;
constexpr int Dim = 3;
constexpr double pi = 3.1415926535898;
constexpr int N_head = 5;
constexpr int N_infor = 9;
constexpr double Si_Si_d2 = 3.45 * 3.45; //����ĵ�λΪ��
constexpr double Si_Si_S_d2 = 6.17 * 6.17;
constexpr double Si_O_d2 = 2.19 * 2.19; 
constexpr double error = 1e-6;
double Al_Al_d = 2.8, O_O_d2 = 0.0;
bool first_layer;
int Num_1 = 0, Num_2 = 0;
int kind, serial, M, N_micro, dt, min_serial;
string tmp, traj_filename = "traj.atom", micro_form, matrix_name, eigen_name, decompose_way = "SVD", if_relative, matrix_kind, finite_size, atom_kind, part, functions, matrix_label = "K_out.dat", eigen_label = "K_eigen.dat";
string dir_read, dir_out, dir_matrix, dir_eigen, dir_system, dir_data = "data", N_str, Pressure, Temperature, Simu_case, out_label, smoothing; // CuZr AlSm SiO2
vector<int> atom_kind_Toluene = { 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2 }; //�ױ�
vector<double> atom_masses_Toluene = { 12.01,  12.01,  12.01,  12.01,  12.01,  12.01, 1.008, 1.008, 1.008, 1.008, 1.008,  12.01, 1.008, 1.008, 1.008 };
vector<string> Temp_list;

inline void reset_file(const string& outname) {
	std::ofstream outputFile(outname, ios::out);
	outputFile.close();
}
inline void distribution(const std::vector<double>& D, string outname, const int& N_Bins) {
	double sum = 0.0;
	for (double value : D) {
		sum += value;
	}
	double mean = sum / D.size();
	double min_value = *std::min_element(D.begin(), D.end()) - error;
	double max_value = *std::max_element(D.begin(), D.end()) + error;
	double bin_size = (max_value - min_value) / N_Bins;
	std::vector<int> histogram(N_Bins, 0);

	// ���ֱ��ͼ
	for (double value : D) {
		int bin_index = static_cast<int>((value - min_value) / bin_size);
		if (bin_index >= N_Bins) bin_index = N_Bins - 1;
		histogram[bin_index]++;
	}
	// ���
	std::ofstream outputFile(outname, ios::app);
	outputFile << mean << " " << min_value << " " << bin_size << " ";
	for (int i = 0; i < N_Bins; ++i) {
		outputFile << histogram[i] << " "; // ��� g(r) ֵ
	}
	outputFile << endl;
	outputFile.close();
}

inline void PBC(const double(&L_h)[Dim], const double(&L)[Dim], double(&dr)[Dim])
{

	for (int i = 0; i < Dim; i++) {
		if (dr[i] > L_h[i]) {
			dr[i] = dr[i] - L[i];
		}
		else if (dr[i] < -L_h[i]) {
			dr[i] = dr[i] + L[i];
		}

	}
	// Check if the vector is now inside the box       
}
inline void PBC_r(double(&r)[Dim], const double(&L)[Dim])
{
	for (int i = 0; i < Dim; i++) {
		if (r[i] > L[i]) {
			r[i] = r[i] - L[i];
		}
		else if (r[i] < zero) {
			r[i] = r[i] + L[i];
		}
		if (r[i] > L[i] || r[i] < zero)
		{
			cerr << i << " " << r[i] << " is out of simulation box."
				<< endl;
			exit(-1);
		}
	}
}

// ����һ��3D����������
struct Vec3 {
	double x, y, z;
	int kind; //ԭ�ӵ�����
	int serial; //ԭ�ӵ����
	// ���������ļӼ���������
	Vec3 operator-(const Vec3& other) const {
		return { x - other.x, y - other.y, z - other.z };
	}

	double dot(const Vec3& other) const {
		return x * other.x + y * other.y + z * other.z;
	}
	// ���������ĳ���
	double magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}
};
bool compareByid(const Vec3& a, const Vec3& b) {
	return a.serial < b.serial;
}
inline double Distance_2(const Vec3& A, const Vec3& B, const double(&L)[Dim]) //[-L/2, L)
{
	double r_2 = zero, diff = zero;
	diff = A.x - B.x;
	diff = diff - L[0] * round(diff / L[0]);
	r_2 += diff * diff;
	diff = A.y - B.y;
	diff = diff - L[1] * round(diff / L[1]);
	r_2 += diff * diff;
	diff = A.z - B.z;
	diff = diff - L[2] * round(diff / L[2]);
	r_2 += diff * diff;
	return r_2;
}
inline double Distance_2(const vector<double>& A, const vector<double>& B, const vector<double>& L) //[-L/2, L)
{
	double r_2 = zero, diff = zero;
	for (int i_Dim = 0; i_Dim < Dim; i_Dim++) {
		diff = A[i_Dim] - B[i_Dim];
		diff = diff - L[i_Dim] * round(diff / L[i_Dim]);
		r_2 += diff * diff;
	}
	return r_2;
}
// ����һ��3D��ʸ������
struct Vec3_wave_vector {
	double x, y, z;
	long double length; //��ʸ�ĳ���
	// ���������ļӼ���������
	Vec3_wave_vector operator-(const Vec3_wave_vector& other) const {
		return { x - other.x, y - other.y, z - other.z };
	}

	double dot(const Vec3_wave_vector& other) const {
		return x * other.x + y * other.y + z * other.z;
	}
	// ���������ĳ���
	double magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}
};
bool compareByid_wave_vector(const Vec3_wave_vector& a, const Vec3_wave_vector& b) {
	return a.length < b.length;
}


inline int CountLines(char* filename)
{
	ifstream read_traj;
	int n = 0;
	string tmp;
	read_traj.open(filename, ios::in);//ios::in ��ʾ��ֻ���ķ�ʽ��ȡ�ļ�
	if (read_traj.fail()) { //�ļ���ʧ��:����0
		return 0;
	}
	else { //�ļ�����
		while (getline(read_traj, tmp)) {
			n++;
		}
		std::cout << "num of lines: " << n << endl;
		return n;
	}
	read_traj.close();
}

class Sect //����ͳ�ƺ͹�������
{
public:
	//~Sect();
	int nint;
	double max;
	double min;
	double leng;
	vector<double> mid;
	vector<double> count;
	vector<double> val;
	vector<double> valstd;
	int Numcount;
	Sect(const double& Min, const double& Max, const int& Nint);
	void set(const double& Int, const double& Val);
	void equ(const double& Length);
	double sum();
};
Sect::Sect(const double& Min, const double& Max, const int& Nint)
{
	nint = Nint;
	mid.resize(nint);
	count.resize(nint);
	val.resize(nint);
	valstd.resize(nint);
	min = Min;
	max = Max;
	leng = (max - min) / nint;
	Numcount = 0;
	for (int i = 0; i < nint; i++) {
		mid[i] = min + i * leng + leng / 2.0;
		val[i] = zero;
		valstd[i] = zero;
		count[i] = zero;
	}
}
void Sect::set(const double& Int, const double& Val)
{
	if ((Int >= min) && (Int < max)) {
		int fract = (int)((Int - min) / leng); // get the integer part
		if ((fract >= 0) && (fract < nint)) {
			val[fract] += Val;
			valstd[fract] += Val * Val;
			count[fract] += 1.0;
			//printf("Int: %f fract: %d pri: %f val[fract]: %f\n", Int, fract, pri, val[fract]);
		}
	}
}
void Sect::equ(const double& Length)
{
	for (int i = 0; i < nint; i++) {
		if (count[i] > 0.5) {
			val[i] = val[i] / count[i];
			valstd[i] = sqrt(valstd[i] / count[i] - val[i] * val[i]);
		}
		if (Length > 0.5) {
			count[i] = count[i] / Length;
		}
	}
}
double Sect::sum()
{
	double sum = zero;
	for (int i = 0; i < nint; i++) {
		sum += count[i];
	}
	Numcount += (int)sum;
	return sum;
}
inline int count_lines_eigen(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Unable to open file: " << filename << std::endl;
		return -1; // ���ش�����
	}

	int line_count = 0;
	std::string line;
	while (std::getline(file, line)) {
		line_count++;
	}

	file.close(); // �ر��ļ�
	return line_count;
}
inline int count_lines(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Unable to open file: " << filename << std::endl;
		return -1; // ���ش�����
	}
	
	int line_count = 0;
	const std::streamsize buffer_size = 4096; // ÿ�ζ�ȡ���ֽ���
	char buffer[buffer_size];

	while (file.read(buffer, buffer_size) || file.gcount() > 0) {
		for (std::streamsize i = 0; i < file.gcount(); ++i) {
			if (buffer[i] == '\n') {
				line_count++;
			}
		}
	}

	file.close(); // �ر��ļ�
	return line_count;
}
int createDir(const std::string& dirName) //��������ļ���
{
#ifdef __unix
	char order[200] = "mkdir -p ";
	strcat(order, dirName.c_str());
	system(order);
	return 0;
#else
	int m = 0, n;
	string str1, str2;
	str1 = dirName;
	if (str1.substr(0, 2) == ".\\") {
		str2 = str1.substr(0, 1);
		str1 = str1.substr(2, str1.size());
	}
	else {
		str2 = str1.substr(0, 2);
		str1 = str1.substr(3, str1.size());
	}
	while (m >= 0)
	{
		m = str1.find('\\');
		str2 += '\\' + str1.substr(0, m);
		n = _access(str2.c_str(), 0); //�жϸ�Ŀ¼�Ƿ����
		if (n == -1)
		{
			if (_mkdir(str2.c_str()) != 0)     //����Ŀ¼
			{
				return 1;
			}
		}
		str1 = str1.substr(m + 1, str1.size());
	}
	return 0;
#endif
}
