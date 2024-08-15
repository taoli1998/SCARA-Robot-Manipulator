#include "function.h"
#include "stdafx.h"
#include <stdio.h>
#include <conio.h>
#include "ensc-488.h"
#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <thread>
#include <vector>
#include <algorithm>
#include <Windows.h>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace std::chrono;

const int resolution = 10;
const double L1 = 195;
const double L2 = 142;
bool eofoutput = false;
double num_solutions = 2;
double station_angle = 0;
static bool Check_joint_limit(double* joint);
static double Correct_Angle(double angle);

double** ForwKin(double* q) {
	 //Inverse Base wrt Station
	double** invT_BS = new double* [4];
	for (int i = 0; i < 4; i++) {
		invT_BS[i] = new double[4];
	}
	invT_BS[0][0] = cos(station_angle);
	invT_BS[0][1] = sin(station_angle);
	invT_BS[0][2] = 0;
	invT_BS[0][3] = 0;

	invT_BS[1][0] = -sin(station_angle);
	invT_BS[1][1] = cos(station_angle);
	invT_BS[1][2] = 0;
	invT_BS[1][3] = 0;

	invT_BS[2][0] = 0;
	invT_BS[2][1] = 0;
	invT_BS[2][2] = 1;
	invT_BS[2][3] = 0;

	invT_BS[3][0] = 0;
	invT_BS[3][1] = 0;
	invT_BS[3][2] = 0;
	invT_BS[3][3] = 1;

	return WHERE(q, invT_BS);
}

double** WHERE(double* q, double** invT_BS)
{
	double num;
	double** T_ST = new double* [4];
	for (int i = 0; i < 4; i++) {
		T_ST[i] = new double[4];
	}

	double** T_BT = KIN(q);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			num = 0;
			for (int k = 0; k < 4; k++) {
				num += invT_BS[i][k] * T_BT[k][j];
			}
			T_ST[i][j] = num;
		}
	}

	return T_ST;

}

double** KIN(double* q)
{
	double theta2 = DEG2RAD(q[0]);
	double theta3 = DEG2RAD(q[1]);
	double theta5 = DEG2RAD(q[3]);
	double d4 = q[2];

	double** T_BT = new double* [4];
	for (int i = 0; i < 4; i++) {
		T_BT[i] = new double[4];
	}
	T_BT[0][0] = cos(theta2 + theta3 - theta5);
	T_BT[0][1] = sin(theta2 + theta3 - theta5);
	T_BT[0][2] = 0;
	T_BT[0][3] = L1 * cos(theta2) + L2 * cos(theta2 + theta3);

	T_BT[1][0] = sin(theta2 + theta3 - theta5);
	T_BT[1][1] = -cos(theta2 + theta3 - theta5);
	T_BT[1][2] = 0;
	T_BT[1][3] = L1 * sin(theta2) + L2 * sin(theta2 + theta3);

	T_BT[2][0] = 0;
	T_BT[2][1] = 0;
	T_BT[2][2] = -1;
	T_BT[2][3] = -d4 - 75;

	T_BT[3][0] = 0;
	T_BT[3][1] = 0;
	T_BT[3][2] = 0;
	T_BT[3][3] = 1;

	return T_BT;
}
 
//double** WHERE(double* joint)
//{
//	double theta2 = DEG2RAD(joint[0]);
//	double theta3 = DEG2RAD(joint[1]);
//	double theta5 = DEG2RAD(joint[3]);
//	double d4 = joint[2];
//
//	double** Matrix = new double* [4];
//	for (int i = 0; i < 4; i++) {
//		Matrix[i] = new double[4];
//	}
//	Matrix[0][0] = cos(theta2 + theta3 - theta5);
//	Matrix[0][1] = sin(theta2 + theta3 - theta5);
//	Matrix[0][2] = 0;
//	Matrix[0][3] = L1 * cos(theta2) + L2 * cos(theta2 + theta3);
//
//	Matrix[1][0] = sin(theta2 + theta3 - theta5);
//	Matrix[1][1] = -cos(theta2 + theta3 - theta5);
//	Matrix[1][2] = 0;
//	Matrix[1][3] = L1 * sin(theta2) + L2 * sin(theta2 + theta3);
//
//	Matrix[2][0] = 0;
//	Matrix[2][1] = 0;
//	Matrix[2][2] = -1;
//	Matrix[2][3] = -d4 - 75;
//
//	Matrix[3][0] = 0;
//	Matrix[3][1] = 0;
//	Matrix[3][2] = 0;
//	Matrix[3][3] = 1;
//
//	return Matrix;
//}

double* SOLVE(double** matrix) {

	double** sol = new double* [2];
	for (int i = 0; i < 2; i++) {
		sol[i] = new double[4];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 4; j++) {
			sol[i][j] = 0;
		}
	}
	
	double Px = matrix[0][3];
	double Py = matrix[1][3];
	double Pz = matrix[2][3];
	double joint1[2];

	double C2 = (Px*Px + Py*Py - (L1*L1 + L2*L2)) / (2 * L1 * L2);
	if (C2 > 1) {
		cout << C2 << endl;
		cout << "There is no solution! for Theta2" << endl;
		num_solutions = -1;
	}
	double S2 = sqrt(1 - C2 * C2);
	double joint2[2] = { atan2(S2, C2), atan2(-S2, C2) };

	for (int i = 0; i < 2; i++) {
		double a = L1 + L2 * cos(joint2[i]);
		double b = L2 * sin(joint2[i]);
		double r = sqrt(a * a + b * b);
		if (r == 0) {
			//cout << "There are infinite solutions!" << endl;
			num_solutions = -999;
		}
		else {
			joint1[i] = atan2(a * Py - b * Px, a * Px + b * Py);
		}
	}

	double joint3 = -Pz - 75;
	sol[0][2] = joint3;
	sol[1][2] = joint3;
	double temp = atan2(matrix[1][0], matrix[0][0]);

	double joint4[2] = { joint1[0] + joint2[0] - temp, joint1[1] + joint2[1] - temp };

	for (int i = 0; i < 2; i++) {
		joint1[i] = RAD2DEG(joint1[i]);
		joint2[i] = RAD2DEG(joint2[i]);
		joint4[i] = RAD2DEG(joint4[i]);
	}

	for (int i = 0; i < 2; i++) {
		joint1[i] = Correct_Angle(joint1[i]);
		joint2[i] = Correct_Angle(joint2[i]);
		joint4[i] = Correct_Angle(joint4[i]);
	}

	sol[0][0] = joint1[0];
	sol[1][0] = joint1[1];
	sol[0][1] = joint2[0];
	sol[1][1] = joint2[1];
	sol[0][3] = joint4[0];
	sol[1][3] = joint4[1];

	for (int i = 0;i < 2;i++) {
		for (int j = 0;j < 4;j++) {
			sol[i][j] = round(sol[i][j]);
			//cout << sol[i][j] << "\t";
		}
		//cout << endl;
	}

	if (Check_joint_limit(sol[1]) && Check_joint_limit(sol[0])) {		
		num_solutions = -1;
	}
	else if (Check_joint_limit(sol[0])) {
		num_solutions = 1;
	}
	else if (Check_joint_limit(sol[1])) {
		num_solutions = 0;
	}

	double* q1 = new double[4];
	JOINT q0 = { 0,0,0,0 };
	bool a = GetConfiguration(q0);

	if (num_solutions == -1)
	{
		//cout << "No Valid Solution found" << endl;
	}
	else if (num_solutions == 0) {
		//cout << "First solution is valid" << endl;
		q1[0] = sol[0][0];
		q1[1] = sol[0][1];
		q1[2] = sol[0][2];
		q1[3] = sol[0][3];
	}
	else if (num_solutions == 1) {
		//cout << "Second solution is valid" << endl;
		q1[0] = sol[1][0];
		q1[1] = sol[1][1];
		q1[2] = sol[1][2];
		q1[3] = sol[1][3];
	}
	else if (num_solutions == 2) {
		//cout << "We have two solution" << endl;
		double dist[2] = {0, 0};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 4; j++) {
				dist[i] += abs(q0[j] - sol[i][j]);
			}
		}
		if (dist[0] < dist[1]) {
			//cout << "Choose first solution" << endl;
			q1[0] = sol[0][0];
			q1[1] = sol[0][1];
			q1[2] = sol[0][2];
			q1[3] = sol[0][3];
		}
		else if (dist[0] > dist[1]) {
			//cout << "Choose second solution" << endl;
			q1[0] = sol[1][0];
			q1[1] = sol[1][1];
			q1[2] = sol[1][2];
			q1[3] = sol[1][3];
		}
		else {
			//cout << "They have the same distance, choose first solution" << endl;
			q1[0] = sol[1][0];
			q1[1] = sol[1][1];
			q1[2] = sol[1][2];
			q1[3] = sol[1][3];
		}
	}
	return q1;
}





void Is_Grasp(void) {
	cout << endl;
	cout << "Press 1 to close the gripper or 2 to open the gripper" << endl;
	char gripper;
	gripper = _getch();
	if (gripper == '1') {
		Grasp(true);
	}
	else if (gripper == '2') {
		Grasp(false);
	}
	else {
		return;
	}
}

double** Joint_To_Tuple(double* joint) {

	double** user_tuple = new double* [4];
	for (int i = 0; i < 4; i++) {
		user_tuple[i] = new double[4];
	}
	user_tuple[0][0] = cos(joint[3]);
	user_tuple[0][1] = sin(joint[3]);
	user_tuple[0][2] = 0;
	user_tuple[0][3] = joint[0];

	user_tuple[1][0] = sin(joint[3]);
	user_tuple[1][1] = -cos(joint[3]);
	user_tuple[1][2] = 0;
	user_tuple[1][3] = joint[1];

	user_tuple[2][0] = 0;
	user_tuple[2][1] = 0;
	user_tuple[2][2] = -1;
	user_tuple[2][3] = joint[2];

	user_tuple[3][0] = 0;
	user_tuple[3][1] = 0;
	user_tuple[3][2] = 0;
	user_tuple[3][3] = 1;

	return user_tuple;
}

void printMatrix(double** matrix) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << matrix[i][j] << "\t";
		}
		cout <<endl;
	}
}


void printMatrix5(double** matrix) {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 4; j++) {
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void printJoint(double* joint) {
	for (int i = 0; i < 4; i++) {
		cout << joint[i] << "\t";
	}
	cout << endl;
}

static bool Check_joint_limit(double *joint) {
	if ((joint[0] > 150) || (joint[0] < -150)){
		//cout << "ERROR: Joint 1 is out of limit" << endl;
		return true;
	}
	else if ((joint[1] > 100) || (joint[1] < -100)) {
		//cout << "ERROR: Joint 2 is out of limit" << endl;
		return true;
	}
	else if ((joint[2] > -100) || (joint[2] < -200)) {
		//cout << "ERROR: Joint 3 is out of limit" << endl;
		return true;
	}
	else if ((joint[3] > 160) || (joint[3] < -160)) {
		//cout << "ERROR: Joint 4 is out of limit" << endl;
		return true;
	}
	else {
		return false;
	}
}

static double Correct_Angle(double angle) {
	while (angle <= -180) {
		angle += 360;
	}
	while (angle > 180) {
		angle -= 360;
	}
	return angle;
}

//double* INKIN(double** matrix) {
//
//	double** T_BS = new double* [4];
//	for (int i = 0; i < 4; i++) {
//		T_BS[i] = new double[4];
//	}
//	T_BS[0][0] = cos(station_angle);
//	T_BS[0][1] = sin(station_angle);
//	T_BS[0][2] = 0;
//	T_BS[0][3] = 0;
//
//	T_BS[1][0] = -sin(station_angle);
//	T_BS[1][1] = cos(station_angle);
//	T_BS[1][2] = 0;
//	T_BS[1][3] = 0;
//
//	T_BS[2][0] = 0;
//	T_BS[2][1] = 0;
//	T_BS[2][2] = 1;
//	T_BS[2][3] = 0;
//
//	T_BS[3][0] = 0;
//	T_BS[3][1] = 0;
//	T_BS[3][2] = 0;
//	T_BS[3][3] = 1;
//
//
//	double** invT_WT = new double* [4];
//	for (int i = 0; i < 4; i++) {
//		invT_WT[i] = new double[4];
//	}
//	invT_WT[0][0] = 1;
//	invT_WT[0][1] = 0;
//	invT_WT[0][2] = 0;
//	invT_WT[0][3] = 0;
//
//	invT_WT[1][0] = 0;
//	invT_WT[1][1] = 1;
//	invT_WT[1][2] = 0;
//	invT_WT[1][3] = 0;
//
//	invT_WT[2][0] = 0;
//	invT_WT[2][1] = 0;
//	invT_WT[2][2] = 1;
//	invT_WT[2][3] = -60;
//
//	invT_WT[3][0] = 0;
//	invT_WT[3][1] = 0;
//	invT_WT[3][2] = 0;
//	invT_WT[3][3] = 1;
//
//	double* q1 = new double[4];
//	return q1;
//
//}
double get_distance(double* joint_vec1, double* joint_vec2) {
	double j0 = joint_vec2[0]-joint_vec1[0];
	double j1 = joint_vec2[1] - joint_vec1[1];
	double j2 = joint_vec2[2] - joint_vec1[2];
	double j3 = joint_vec2[3] - joint_vec1[3];

	return sqrt(pow(j0,2.0)+pow(j1,2.0)+pow(j2,2.0)+pow(j3,2.0));

}


bool Check_Workspace(double a, double b, double c, double d) {
	if (sqrt(a * a + b * b) > 337 || sqrt(a * a + b * b) < 53) {
		return false;
	}
	return true;
}

int num_sample(double timef)
{
	//here the resolution is define to be 10, can change later
	int samples = (int)ceil(timef * resolution);
	return samples;
}

double intervals(double time, int samples)
{
	if (samples == 0) {
		return 0;
	}
	//here the resolution is define to be 10, can change later
	double interval = time / double(samples);
	return interval;
}

void export_csv(double** matrix) {
	ofstream myFile("C:/Users/alons/Desktop/Amirali's 488 Project/txt.csv");
	myFile << "0,337\n";
	for (int i = 0; i < 4; i++) {
		myFile << matrix[i][1];
		myFile << ",";
		myFile << matrix[i][2];
		myFile << "\n";
	}
	myFile.close();
}

void Trajectory_Planner(double** traj_joint, double timef) {
	
	cout << endl << "* Inside Trajectory Planner *" << endl;

	bool jlim_invalid = false;
	for (int i = 0; i < 4; i++){
		jlim_invalid = Check_joint_limit(traj_joint[i]);
		if (jlim_invalid) {
			cout << "Joint angle/angles in frame " << i+1 << " exceed limits! Please try different values." << endl;
			return;
		}
	}
	

	///////////////// First, we will find the time required between each frame based on joint limits /////////////////
	// Initialize variable for time based on limits
	double** tf_min = new double* [4]; // Create a matrix for tf
	for (int i = 0; i < 4; i++) {
		tf_min[i] = new double[4];
	}

	double* time = new double [4]; // Create a matrix for time in each spline


	// Initialize variable for number of samples required
	int* samples = new int [4]; // Create a matrix for samples

	// Initialize variable for the time interval between each sample
	double* interv = new double [4]; // Create a matrix for intervals
	
	
	// Holder values for calculating tf_max vel and acc
	double tf_vel;
	double tf_acc;
	
	// Iterate through each spline and calculate the tf for each joint based on joint limits
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {

			if (j != 2) {
				tf_vel = abs((double(1) / double(100)) * (traj_joint[i+1][j] - traj_joint[i][j]));
				tf_acc = sqrt(abs((double(1) / double(100)) * (traj_joint[i+1][j] - traj_joint[i][j])));
			}
			else {
				tf_vel = abs((double(3) / double(100)) * (traj_joint[i+1][j] - traj_joint[i][j]));
				tf_acc = sqrt(abs(((double(3) / double(100)) * (traj_joint[i+1][j] - traj_joint[i][j]))));
			}

			if (tf_vel >= tf_acc) {
				tf_min[i][j] = tf_vel;
			}
			else {
				tf_min[i][j] = tf_acc;
			}

			//samples[i][j] = num_sample(tf_min[i][j]);
			//interv[i][j] = intervals(tf_min[i][j],samples[i][j]);
		}
	}

	double temp_tf[4] = {0, 0, 0, 0};
	bool timeflaw = false;
	// Check each joint to see if input tf is unreasonable
	for (int j = 0; j < 4; j++) {

		for (int i = 0; i < 4; i++) {
			temp_tf[j] = temp_tf[j] + tf_min[i][j];
		}
	}

	double maxtime_min = 0;
	for (int i = 0; i < 4; i++) {
		maxtime_min = max(maxtime_min, temp_tf[i]);
	}
		
	if (maxtime_min > timef) {
		// If there is a time flaw in the input, then return a warning and stop planning
		cout << " WARNING: Time input violates joint velocity/acceleration limits! Try another input." << endl;
		cout << " Make sure that the input time is greater than " << maxtime_min << "." << endl;
		return;
	}

	cout << "Minimum time: " << maxtime_min << endl;

	double tot_dist = get_distance(traj_joint[0], traj_joint[1]) + get_distance(traj_joint[1], traj_joint[2]) + get_distance(traj_joint[2], traj_joint[3]) + get_distance(traj_joint[3], traj_joint[4]);
	double time_ratio = 0.0;
	double total_segment_time = 0.0;

	for (int i = 0; i < 4; i++) {

		time_ratio = get_distance(traj_joint[i], traj_joint[i+1]) / tot_dist;

		time[i] = (timef - maxtime_min) * time_ratio + maxtime_min;

		total_segment_time += time[i];

	}
	
	for (int i = 0; i < 4; i++) {
		time[i] *= timef / total_segment_time;

		samples[i] = num_sample(time[i]);

		interv[i] = intervals(time[i], samples[i]);

		cout << "Segment " << i + 1 << ": time of " << time[i] << ", " << samples[i] << " samples, for " << interv[i] << " intervals." << endl;
		if (time[i] < maxtime_min) {
			// If there is a time flaw in the input, then return a warning and stop planning
			cout << " WARNING: The calculated segment time is less than the required time." << endl;
			cout << " Setting segmet time to minimum time of " << maxtime_min << ". This may increase total time of travel." << endl;
			time[i] = maxtime_min;
			//return;
		}
	}

	

	// Now, for each spline, we have a time vector that determines the minimum time needed for each joint and samples per spline

	///////////////// Second, we will find the trajectory between each frame for each joint /////////////////

	// Initialize matrix for trajectory functions - each will hold the coefficients
	double*** traj_coeffs = new double** [4]; // Create a matrix for tf
	for (int i = 0; i < 4; i++) {
		traj_coeffs[i] = new double*[4];
		for (int j = 0; j < 4; j++) {
			traj_coeffs[i][j] = new double[4];
		}
	}

	// Find function coefficients for position and velocity
	for (int i = 0; i < 4; i++) { // Each trajectory frame to frame

		for (int j = 0; j < 4; j++) { // Each spline for each joint

			if (i == 0) {
				traj_coeffs[i][j] = Find_Spline(traj_joint[i][j], traj_joint[i + 1][j], time[i]);
			}
			else if (i < 2) {
				traj_coeffs[i][j] = Find_Spline(traj_joint[i][j], traj_joint[i + 1][j], time[i]);
			}
			else {
				traj_coeffs[i][j] = Find_Spline(traj_joint[i][j], traj_joint[i + 1][j], time[i]);
			}

		}
	}

	// Find max size within each sample set
	/*int max_array[4] = {0,0,0,0};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			max_array[i] = max(max_array[i], samples[i]);
		}
	}*/

	// Find max interval within each sample set
	/*double max_interval[4] = { 0,0,0,0 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			max_interval[i] = max(max_interval[i], interv[i]);
		}
	}*/

	// Find position, velocity, and acceleration of each spline for each joint
	//double*** pos = new double** [4];
	//double*** vel = new double** [4];
	//double*** acc = new double** [4];
	//double time_interval;

	//for (int i = 0; i < 4; i++) {
	//	
	//	pos[i] = new double* [4];
	//	vel[i] = new double* [4];
	//	acc[i] = new double* [4];

	//	for (int j = 0; j < 4; j++) {

	//		pos[i][j] = new double[samples[i][j]];
	//		vel[i][j] = new double[samples[i][j]];
	//		acc[i][j] = new double[samples[i][j]];

	//		time_interval = interv[i][j];

	//		for (int k = 0; k < samples[i][j]; k++) {
	//			
	//			pos[i][j][k] = traj_coeffs[i][j][0] + (traj_coeffs[i][j][1]) * (time_interval)+(traj_coeffs[i][j][2]) * (pow(2.0, time_interval)) + (traj_coeffs[i][j][3]) * (pow(3.0, time_interval));
	//			vel[i][j][k] = traj_coeffs[i][j][1] + 2*(traj_coeffs[i][j][2]) * (time_interval) + 3*(traj_coeffs[i][j][3]) * (pow(2.0, time_interval));
	//			acc[i][j][k] = 2*(traj_coeffs[i][j][2])  + 6*(traj_coeffs[i][j][3]) * (time_interval);

	//			time_interval = time_interval + time_interval;
	//		}
	//	}

	//}

	vector<vector<vector<double>>> pos(4, vector<vector<double>>(4));
	vector<vector<vector<double>>> vel(4, vector<vector<double>>(4));
	vector<vector<vector<double>>> acc(4, vector<vector<double>>(4));
	double time_interval;
	//double temp_pos;
	//double temp_vel;
	//double temp_acc;

	// Find max sample within each joint set and basically fill the entire vector with the same amount of samples
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {

			time_interval = interv[i];
			//cout << max_array[i] << endl;
			for (int k = 0; k < samples[i]; k++) {
				
				if (k < samples[i]) { // If samples are less than the max sample, then calculate the proper value for each joint
					pos[i][j].push_back(floor(traj_coeffs[i][j][0] + (traj_coeffs[i][j][1]) * (time_interval)+(traj_coeffs[i][j][2]) * (pow(time_interval, 2.0)) + (traj_coeffs[i][j][3]) * (pow(time_interval, 3.0))));
					vel[i][j].push_back(floor(traj_coeffs[i][j][1] + 2 * (traj_coeffs[i][j][2]) * (time_interval)+3 * (traj_coeffs[i][j][3]) * (pow(time_interval, 2.0))));
					acc[i][j].push_back(floor(2 * (traj_coeffs[i][j][2]) + 6 * (traj_coeffs[i][j][3]) * (time_interval)));
				}
				else { // Else, push back in the final values that were reached

					pos[i][j].push_back(traj_joint[i][j]);
					vel[i][j].push_back(0);
					acc[i][j].push_back(0);
				}

				//for testing
				if (j == 2) {
					if (vel[i][j][k] > 50.0000000000 || vel[i][j][k] < -50.0000000000) {
						cout << "Joint " << j << ":" << endl;
						cout << std::fixed << std::setprecision(32) << vel[i][j][k] << " at frame " << i + 1 << ", at sample point " << k << endl;
					}
				}


				time_interval = time_interval + interv[i];

			}
		}
	}


	// Now that we know the trajectory of each spline, we can iterate moving from one frame to another
	
	JOINT joint_p;
	JOINT joint_v;
	JOINT joint_a;

	//JOINT joint_v_end = {0, 0, 0, 0};
	//JOINT joint_a_end = {0, 0, 0, 0};
	
	double time_sleep = 0;

	// Define a clock for elapsed time
	auto start = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < 4; i++) {
		for (int k = 0; k < samples[i]; k++) {

			time_sleep = interv[i];
			/*std::chrono::duration<double>(time_sleep);*/

			for (int n = 0; n < 4; n++) {
				joint_p[n] = pos[i][n][k];
				joint_v[n] = vel[i][n][k];
				joint_a[n] = acc[i][n][k];
			}


			/*start_time = steady_clock::now();*/

			MoveWithConfVelAcc(joint_p, joint_v, joint_a);

			//std::this_thread::sleep_for(std::chrono::duration<double>(time_sleep));
			Sleep(time_sleep*1000);
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
	std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds." << std::endl;

	StopRobot();

	// For plotting position
	int num_of_samples = samples[0] + samples[1] + samples[2] + samples[3];
	//double time_interval_plot = 0;
	//
	//JOINT temp_pos_joint;
	//double** cart_arr = new double*[4];

	//ofstream myFile("C:/Users/alons/Desktop/Amirali's 488 Project/position.csv");
	//myFile << "0,337\n";

	//for (int i = 0; i < 4; i++) {
	//	cart_arr[i] = new double[4];
	//}

	//if (myFile.is_open()) {

	//	for (int i = 0; i < 4; i++) {

	//		for (int k = 0; k < samples[i]; k++) {

	//			for (int n = 0; n < 4; n++) {
	//				temp_pos_joint[n] = pos[i][n][k];
	//				/*joint_v[n] = vel[i][n][k];
	//				joint_a[n] = acc[i][n][k];*/
	//			}

	//			cart_arr = ForwKin(temp_pos_joint);

	//			myFile << time_interval_plot;
	//			myFile << ",";
	//			myFile << pos[i][0][k];
	//			myFile << ",";
	//			myFile << pos[i][1][k];
	//			myFile << ",";
	//			myFile << pos[i][2][k];
	//			myFile << ",";
	//			myFile << pos[i][3][k];
	//			myFile << ",";
	//			myFile << cart_arr[0][3];
	//			myFile << ",";
	//			myFile << cart_arr[1][3];
	//			myFile << ",";
	//			myFile << "\n";

	//			time_interval_plot = time_interval_plot + interv[i];
	//			
	//		}
	//		
	//	}
	//	myFile.close();

	//}
	//else {
	//	cout << "Issues opening file. Please try again!" << endl;
	//	myFile.close();
	//}

	//// Velocity plot
	//time_interval_plot = 0;

	//ofstream myVelFile("C:/Users/alons/Desktop/Amirali's 488 Project/velocity.csv");
	//myVelFile << "0,337\n";

	//if (myVelFile.is_open()) {

	//	for (int i = 0; i < 4; i++) {

	//		for (int k = 0; k < samples[i]; k++) {

	//			myVelFile << time_interval_plot;
	//			myVelFile << ",";
	//			myVelFile << vel[i][0][k];
	//			myVelFile << ",";
	//			myVelFile << vel[i][1][k];
	//			myVelFile << ",";
	//			myVelFile << vel[i][2][k];
	//			myVelFile << ",";
	//			myVelFile << vel[i][3][k];
	//			myVelFile << ",";
	//			myVelFile << "\n";

	//			time_interval_plot = time_interval_plot + interv[i];

	//		}

	//	}
	//	myVelFile.close();

	//}
	//else {
	//	cout << "Issues opening file. Please try again!" << endl;
	//	myVelFile.close();
	//}
	//	
	//
	//// Acceleration plot
	//time_interval_plot = 0;

	//ofstream myAccFile("C:/Users/alons/Desktop/Amirali's 488 Project/acceleration.csv");
	//myAccFile << "0,337\n";

	//if (myAccFile.is_open()) {

	//	for (int i = 0; i < 4; i++) {

	//		for (int k = 0; k < samples[i]; k++) {

	//			myAccFile << time_interval_plot;
	//			myAccFile << ",";
	//			myAccFile << acc[i][0][k];
	//			myAccFile << ",";
	//			myAccFile << acc[i][1][k];
	//			myAccFile << ",";
	//			myAccFile << acc[i][2][k];
	//			myAccFile << ",";
	//			myAccFile << acc[i][3][k];
	//			myAccFile << ",";
	//			myAccFile << "\n";

	//			time_interval_plot = time_interval_plot + interv[i];

	//		}

	//	}
	//	myAccFile.close();

	//}
	//else {
	//	cout << "Issues opening file. Please try again!" << endl;
	//	myAccFile.close();
	//}


	cout << "Total samples = " << num_of_samples << endl;

}

double* Find_Spline(double init_joint, double dest_joint, double timef) {

	double t_init = 0;
	double* a_i = new double[4];
	
	if (timef != 0) {
		a_i[0] = init_joint;
		a_i[1] = 0;
		a_i[2] = (double(3) / (pow(timef, 2.0))) * (dest_joint - init_joint);
		a_i[3] = -(double(2) / (pow(timef, 3.0))) * (dest_joint - init_joint);
	}
	else {
		a_i[0] = 0;
		a_i[1] = 0;
		a_i[2] = 0;
		a_i[3] = 0;
	}
	

	return a_i;
}