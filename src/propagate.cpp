#include "propagate.h"

using namespace std;

const double k[9] = { 0.00287729, 0.0023821, -0.0010625, 0.000198502, -1.88428e-05, 9.70521e-07, -2.70559e-08,  3.7516e-10, -1.94316e-12 };

void eqns_of_motion(const vector<double> &x, vector<double> &dxdt, double t)
{
	/* Central-Force equations of motion for galactic motion.
	State vector 'x' is 6-elements containing position (kpc) and velocity (kpc/Myr)
	Output 'dxdt' is 6-elements containing velocity (kpc/Myr) and acceleration (kpc/Myr^2)
	*/
	double vel_units = 31557600e6 / 30856775814671900.0;  // conversion between km/s and kpc/Myr
	double r, vc, denom, fr;


	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);  // kpc

	denom = k[8];
	for (int i = 7; i >= 0; i--)
		denom = denom * r + k[i];  // (km/s)^-1

	vc = 1.0 / denom;  // kpc/Myr
	fr = vc * vc / r;  // kpc/Myr^2

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	dxdt[3] = -x[0] / r * fr;
	dxdt[4] = -x[1] / r * fr;
	dxdt[5] = -x[2] / r * fr;
}
