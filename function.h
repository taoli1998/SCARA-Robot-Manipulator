#ifndef _function
#define _function

double** ForwKin(double* q);
double** WHERE(double* q, double** invT_BS);
double** KIN(double* q);

double* SOLVE(double** matrix);
void Is_Grasp(void);
void printMatrix(double** matrix);
void printJoint(double* joint);
double** Joint_To_Tuple(double* joint);

bool Check_Workspace(double a, double b, double c, double d);
void printMatrix5(double** matrix);
int num_sample(double timef);
void Trajectory_Planner(double** traj_joint, double timef);
double* Find_Spline(double init_joint, double dest_joint, double timef);
double intervals(double time, int samples);
double get_distance(double* joint_vec1,double* joint_vec2);
void export_csv(double** matrix);
#endif 

