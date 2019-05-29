#include <cmath>
#include "common.h"
#include "propagate.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void eqns_of_motion(const vector<double> &x, vector<double> &dxdt, double t)
{
	/* Central-Force equations of motion for galactic motion.
	 * State vector 'x' is 6-elements containing position (kpc) and velocity (kpc/Myr)
	 * Output 'dxdt' is 6-elements containing velocity (kpc/Myr) and acceleration (kpc/Myr^2) */
	double r, vc, fr;


	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);  // kpc

	vc = KM_PER_SEC_TO_KPC_PER_MYR / _kr(r);  // kpc/Myr
	fr = vc * vc / r;  // kpc/Myr^2

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	dxdt[3] = -x[0] / r * fr;
	dxdt[4] = -x[1] / r * fr;
	dxdt[5] = -x[2] / r * fr;
}

void eqns_of_motion_stm(const vector<double> &x, vector<double> &dxdt, double t)
{
	/* Central-Force equations of motion including a State-Transition Matrix (STM).
	 * State vector 'x' is 42-elements contiaining position (kpc), velocity (kpc/Myr), 
	 * and 36 elements representing the state transition matrix. Output 'dxdt' contains
	 * velocity, acceleration (kpc/Myr^2) and the derivative of the elements of the STM. */
	
	double r, r2, vc, kr, dvcdx, dqdx;

	// Populate the velocity and acceleration terms
	eqns_of_motion(x, dxdt, t);

	// Create a matrix to store the STM.
	vector<double>::const_iterator first = x.begin() + 6;
	vector<double>::const_iterator last = x.begin() + 42;
	vector<double> phi_vec(first, last)
	MatrixXd phi(phi_vec.data()); jacobian(6,6), phi_dot(6,6);
	phi.resize(6, 6)  // Resize the input STM to a 6x6 matrix.

	// Build the Jacobian of the equations of motion.
	MatrixXd jacobian = MatrixXd::Zero(6, 6);
	jacobian(0, 3) = 1.0;
	jacobian(1, 4) = 1.0;
	jacobian(2, 5) = 1.0;
	
	// Analytic expression for derivative of acceleration w.r.t. position.
	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);  // kpc
	r2 = r*r;
	kr = _kr(r);
	vc = KM_PER_SEC_TO_KPC_PER_MYR / kr;

	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
		{
			dvcdx = -KM_PER_SEC_TO_KPC_PER_MYR * _dkdx(r, x[j]) / (kr*kr)
			dqdx = (2.0*vc / r2) * (dvcdx - vc*x[j] / r2);
			jacobian(i+3, j) = -x[i]*dqdx;
			if (i == j)
				jacobian(i+3, j) -= vc*vc / r2;
		}

	// Compute the derivative of the STM.
	MatrixXd phi_dot = jacobian * phi;	
	
	// Insert the elements of the STM derivative into the state vector.
	phi_dot.resize(36);
	for (int i=0; i<36; i++)
		dxdt[i+6] = phi_dot(i);		
}


double _kr(double r)
{
	double val = K8;
	for (int n = 7; n >= 0; n--)
		val = val * r + K[n];  // (km/s)^-1

	return val;
}


double _dkdx(double r, double x)
{
	double val = 0.0;
	for (int n = 1; n < 9; n++)
		val += n*K[n]*pow(r, n-2);

	val *= x;
	return val;
}
