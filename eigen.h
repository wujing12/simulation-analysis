#pragma once
constexpr int max_cols = 10;
inline bool containsNaNInf(const Eigen::MatrixXf& mat) {
	for (int i = 0; i < mat.rows(); ++i) {
		for (int j = 0; j < mat.cols(); ++j) {
			if (std::isnan(mat(i, j)) || std::isinf(mat(i, j))) {
				return true;  // 如果找到 NaN 或者 Inf 元素，则返回 true
			}
		}
	}
	return false;  
}
inline void Matrix_processing(Eigen::MatrixXf& Matrix_run, vector<double>& mean_class, vector<double>& mean_class_2, const int M) {
	if (project == "super_critical") {
		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化			
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				mean_class[0] += Matrix_run(j, i);
			}
		}
		mean_class[0] /= (Matrix_run.cols() * N_micro);
		std::cout << mean_class[0] << endl;
		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				Matrix_run(j, i) = Matrix_run(j, i) - mean_class[0];
			}
		}
	}
	else if (project == "AlSm") {	
		int num_kind = stoi(micro_form);
		cout << "num_kind: " << num_kind << endl;
		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				int micro_kind = j % num_kind;
				if (micro_kind < 3) {
					mean_class_2[0] += Matrix_run(j, i) * Matrix_run(j, i);
				}
				else {
					mean_class_2[1] += Matrix_run(j, i) * Matrix_run(j, i);
				}
			}
		}
		mean_class_2[0] = sqrt(mean_class_2[0] / (M * 9000.0)); //Al
		mean_class_2[1] = sqrt(mean_class_2[1] / (M * 1000.0)); //Sm
		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				int micro_kind = j % num_kind;
				if (micro_kind < 3) {
					Matrix_run(j, i) = Matrix_run(j, i) / mean_class_2[0];
				}
				else {
					Matrix_run(j, i) = Matrix_run(j, i) / mean_class_2[1];
				}
			}
		}
	}
	else {
		int num_kind = stoi(micro_form);  // 将 string 转为 int

		int N_one = N_micro / num_kind;
		cout << "num_kind: " << num_kind << ". N_pro: " << N_one << endl;
		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				int micro_kind = j % num_kind;
				mean_class_2[micro_kind] += Matrix_run(j, i) * Matrix_run(j, i);
				mean_class[micro_kind] += Matrix_run(j, i);
			}
		}
		for (size_t j = 0; j < mean_class.size(); j++) {
			mean_class[j] = mean_class[j] / (Matrix_run.cols() * N_one);
			mean_class_2[j] = mean_class_2[j] / (Matrix_run.cols() * N_one);
			std::cout << mean_class[j] << " " << mean_class_2[j] << endl;
		}

		for (int i = 0; i < Matrix_run.cols(); i++) { //对列进行标准化
			for (size_t j = 0; j < Matrix_run.rows(); j++) {
				int micro_kind = j % num_kind;
				Matrix_run(j, i) = (Matrix_run(j, i) - mean_class[micro_kind]) /
					sqrt(mean_class_2[micro_kind] - mean_class[micro_kind] * mean_class[micro_kind]);
			}
		}
	}
	if (matrix_kind != "Full" && project == "supercooled") {
		std::vector<int> keepRows;
		if (matrix_kind == "3") {
			for (int row = 0; row < Matrix_run.rows(); ++row) {
				int remainder = row % 7;
				if (remainder == 0 || remainder == 1 || remainder == 2) {
					keepRows.push_back(row);
				}
			}
		}
		else if (matrix_kind == "6") {
			for (int row = 0; row < Matrix_run.rows(); ++row) {
				int remainder = row % 7;
				if (remainder == 3 || remainder == 4 || remainder == 5) {
					keepRows.push_back(row);
				}
			}
		}
		else if (matrix_kind == "7") {
			for (int row = 0; row < Matrix_run.rows(); ++row) {
				int remainder = row % 7;
				if (remainder == 6) {
					keepRows.push_back(row);
				}
			}
		}
		Eigen::MatrixXf newMatrix(keepRows.size(), Matrix_run.cols());
		for (int i = 0; i < static_cast<int>(keepRows.size()); ++i) {
			//std::cout << keepRows[i] << " ";
			newMatrix.row(i) = Matrix_run.row(keepRows[i]);
		}
		//std::cout << endl;
		Matrix_run.resize(keepRows.size(), Matrix_run.cols());
		Matrix_run = newMatrix;
		newMatrix.resize(0, 0);
	}
}
inline void get_eigen(const string matrix_name, const string eigen_name) {

	const int M = count_lines_eigen(matrix_name);
	std::cout << "Name: " << matrix_name << ", N_micro: " << N_micro << ", M: " << M << endl;
	Eigen::MatrixXf Matrix_run;
	Matrix_run.resize(N_micro, M); //每列是一个构型的微观态

	ifstream Matrixdata(matrix_name);
	double V = zero;

	for (int i = 0; i < Matrix_run.cols(); i++) { //每行是一个构型的微观态
		if (project == "super_critical") {
			Matrixdata >> V; //格点体积
			//cout << V << " ";
		}
		for (size_t j = 0; j < Matrix_run.rows(); j++) {
			Matrixdata >> Matrix_run(j, i);
			if (project == "super_critical") {
				Matrix_run(j, i) = Matrix_run(j, i) / V;
			}
		}
	}
	//cout << endl;
	Matrixdata.close();
	if (containsNaNInf(Matrix_run)) {
		std::cout << "Matrix contains NaN or Inf elements!" << std::endl;
		exit(-1);
	}
	else {
		std::cout << "Matrix does not contain NaN or Inf elements." << std::endl;
	}

	vector<double> mean_class(7, zero);
	vector<double> mean_class_2(7, zero);
	Matrix_processing(Matrix_run, mean_class, mean_class_2, M);

	if (decompose_way == "Eigen") {
		double norm = Matrix_run.norm();
		std::cout << Matrix_run.rows() << " " << Matrix_run.cols() << " " << norm << endl;
		Matrix_run.normalize(); //除以矩阵的Frobenius范数
		Eigen::MatrixXf C = Matrix_run * Matrix_run.transpose(); //关联矩阵
		Matrix_run.resize(0, 0);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(C);
		if (eigensolver.info() != Eigen::Success) {
			std::cerr << "Eigen decomposition failed.\n";
			exit(-1);
		}
		Eigen::VectorXf eigenvalues = eigensolver.eigenvalues();
		const Eigen::MatrixXf& eigenvectors = eigensolver.eigenvectors();
		
		std::ofstream output_eigen(eigen_name, ios::out);
		for (int i = 0; i < eigenvectors.rows(); ++i) {
			for (int j = 0; j < max_cols; ++j) {
				output_eigen << sqrt(static_cast<double>(N_micro)) * eigenvectors(i, eigenvectors.rows() - j - 1) << " ";
			}
			output_eigen << pow(eigenvalues(eigenvectors.rows() - i - 1), 0.5) << std::endl;
		}
	}
	else if (decompose_way == "SVD") {
		double norm = Matrix_run.norm();
		std::cout << Matrix_run.rows() << " " << Matrix_run.cols() << " " << norm << endl;
		Matrix_run.normalize(); //除以矩阵的Frobenius范数
		int K = Matrix_run.rows();
		if (Matrix_run.cols() < K) {
			K = Matrix_run.cols();
		}
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(
			Matrix_run,
			Eigen::ComputeThinU | Eigen::ComputeThinV
		); // ComputeThinU | ComputeThinV Eigen::ComputeFullU | Eigen::ComputeFullV BDCSVD的精度更低，但更快

		std::cout << "SVD Completed!" << endl;
		Matrix_run.resize(0, 0);
		Eigen::VectorXf singular_values = svd.singularValues(); //列向量
		const Eigen::MatrixXf& left_singular_vectors = svd.matrixU(); //double类型矩阵
		const Eigen::MatrixXf& right_singular_vectors = svd.matrixV();
		//std::cout.precision(2);
		//std::cout << svd.matrixV() << endl;

		std::ofstream output_eigen(eigen_name, ios::out);
		// 输出左奇异值矩阵和本征值
		for (int i = 0; i < left_singular_vectors.rows(); ++i) {
			for (int j = 0; j < max_cols; ++j) {
				output_eigen << sqrt(static_cast<double>(N_micro)) * left_singular_vectors(i, j) << " ";
			}
			if (i < K) {
				output_eigen << singular_values(i) << std::endl;
			}
			else {
				output_eigen << "0" << std::endl;
			}
		}
		// 输出右奇异值矩阵
		for (int i = 0; i < right_singular_vectors.rows(); ++i) {
			for (int j = 0; j < max_cols; ++j) {
				output_eigen << sqrt(static_cast<double>(M)) * right_singular_vectors(i, j) << " ";
			}
			output_eigen << "0" << std::endl;
		}

		/*	Eigen::MatrixXf Sigma = Eigen::MatrixXf::Zero(left_singular_vectors.cols(), right_singular_vectors.cols());
		for (int i = 0; i < singular_values.size(); ++i) {
			Sigma(i, i) = singular_values(i);
		}
		Eigen::MatrixXf A_reconstructed = left_singular_vectors * Sigma * right_singular_vectors.transpose();
		//std::cout << A_reconstructed;
		// 计算重构误差
		double reconstruction_error = (Matrix_run - A_reconstructed).norm();
		std::cout << "Error of Frobenius norm" << reconstruction_error << std::endl;*/
		// 重构矩阵 A_reconstructed = U * Σ * V^T


		//std::cout << "Matrix U (Left singular vectors):\n" << left_singular_vectors << std::endl;
		//std::cout << "Matrix V (Right singular vectors):\n" << right_singular_vectors << std::endl;
		output_eigen.close();
	}
	else {
		std::ofstream output_Matrix(matrix_name, ios::out);
		int num_kind = stoi(micro_form);
		for (int i = 0; i < Matrix_run.rows(); ++i) {
			for (int j = 0; j < Matrix_run.cols(); ++j) {
				output_Matrix << Matrix_run(i, j) << " ";
			}
			output_Matrix << std::endl;
		}
		output_Matrix.close();
	}
}