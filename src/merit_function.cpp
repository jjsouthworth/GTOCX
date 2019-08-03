#define _USE_MATH_DEFINES

#include <cmath>
#include "common.h"
#include "events.h"

using namespace std;


void export_solution_file();
void parse_solution_file();
void compute_merit_function();

void merit_function(Solution inputs, double merit)
{
	/* Computes the merit function. The task of GTOCX is to design trajectories
		 of the ships to maximize the merit function described in the problem
		 statement.
	 * inputs 'inputs' are inputs //TODO: populate these
	 * Output 'merit' is a float containing the merit function value
   */
	double J, B, E_r, E_theta, dV_max, dV_used;
	int N;
	double r, theta;
	double R_k, Theta_k;
	double f_r, f_theta, g_r, g_theta;
	double t_start, t_end, t_submission;

	// What are 'inputs'?
	// TODO: need some sort of settled stars object/struct
	// TODO: dV values
	// TODO: time values


	// Radial and angular error functions
	E_r = 0;
	E_theta = 0;
	for (int k=0, k<=32, k++){

		// Ingredients for summations
		R_k = k + 2.0;  // kpc
		Theta_k = -M_PI + 2.0 * M_PI * k / 32.0;  // rad
		f_and_g(R_k, Theta_k, f_r, f_theta, g_r, g_theta);

		// Summations have different limits
		if (k <= 30){
			E_r += pow(f_r / g_r - 1.0, 2);
		}
		E_theta += pow(f_theta / g_theta) - 1.0, 2);
	}

	// Time bonus factor (for submitting early)
	B = 1.0 + pow((t_end - t_submission) / (t_end - t_start), 4);  // TODO: is pow() best?

	// Top level merit definition
	J = B * (N / (1.0e4 * N * (E_r + E_theta))) * (dV_max / dV_used);

	// Return merit value
	double merit = J;
}


void f_and_g(double N, double r, double theta, double f_r, double f_theta, double g_r, double, g_theta)
{
	/* f and g are defined in the problem statement // TODO: complete this doc string
	 */
  double f_r, f_theta, g_r, g_theta;

	double alpha, beta;

	double r_i, theta_i;
	double s_r, s_theta;

	// Scaling functions are relatved piecewise
	// TODO: the operators on these inequalities seem odd (should there be more gt/lt?)
	if (r == 2.0){
		alpha = 0.5833;
	}
	else if (r == 32.0){
		alpha = 0.4948;
	}
	else{
		alpha = 1.0;
	}
	if (theta == -M_PI){
		beta = 0.5;
	}
	else if (theta == M_PI){
		beta = 0.5;
	}
	else{
		beta = 1.0;
	}

	// f functions
	f_r = 0;
	f_theta = 0;
	for (int i=0, i<=N, i++){
		r_i = ;  // orbital radius of ith settled star, TODO: how to get?
		theta_i = ; // final polar angle (theta_f) of the ith settled star (in stars.txt file) TODO: how to get this?
		s_r = 1.0;  // kpc
		s_theta = 2.0 * M_PI / 32.0;  // rad
		f_r += f(r, r_i, s_r);
		f_theta += f(theta, theta_i, s_theta);
	}
	f_r /= N;
	f_theta /= N;

	// g functions
	R_max = ;
	R_min = ;
	g_r = alpha * 2.0 * r / (pow(R_max, 2) - pow(R_min, 2));
	g_theta = beta * 1.0 / (2.0 * M_PI);

}


void f(double x, double mu, double s, double f)
{
	/* TODO: write this doc string
	*/
	double f;
	if (abs(x-mu) >= s){
		f = 0.0;
	}
	else{
		f = 1.0 / s - abs(x-mu) / pow(s, 2);
	}
}
