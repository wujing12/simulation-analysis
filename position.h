#pragma once
const int write_time = 1000;
void get_trajectory(Trajectory& trajectory, const vector<vector<double>>& L_traj) {
	string outname_posi = dir_out + "traj.txt";
	unwrap_positions(trajectory, L_traj);

	std::ofstream output_posi(outname_posi, ios::out);
	for (int i = 0; i < trajectory.size(); i++) {
		output_posi << L_traj[i][0] << " " << L_traj[i][1] << " " << L_traj[i][2] << endl;
	}
	for (int i = 0; i < trajectory[0].size(); ++i) {
		for (int j = 0; j < Dim - 1; ++j) {
			output_posi << trajectory[write_time][i][j] << " ";
		}
		output_posi << trajectory[write_time][i][Dim - 1] << endl;
	}

	std::vector<std::array<double, Dim>> velocity = generate_velocity(trajectory[write_time], trajectory[write_time + time_interval]);
	for (int i = 0; i < trajectory[0].size(); ++i) {
		for (int j = 0; j < Dim - 1; ++j) {
			output_posi << velocity[i][j] << " ";
		}
		output_posi << velocity[i][Dim - 1] << endl;
	}
	output_posi.close();
}
void get_position(const std::vector<Vec3>& positions, const double(&L)[Dim]) {
	string outname_posi = dir_out + "posi.txt";
	std::ofstream output_posi(outname_posi, ios::app);
	output_posi << L[0] << " " << L[1] << " " << L[2] << endl;

	for (int i = min_serial; i < positions.size(); ++i) {
		output_posi << positions[i].x << " " << positions[i].y << " " << positions[i].z << endl;
	}
	output_posi.close();
}
