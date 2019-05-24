#include <cmath>
#include "common.h"
#include "propagate.h"

using namespace std;


void eqns_of_motion(const vector<double> &x, vector<double> &dxdt, double t)
{
	/* Central-Force equations of motion for galactic motion.
	State vector 'x' is 6-elements containing position (kpc) and velocity (kpc/Myr)
	Output 'dxdt' is 6-elements containing velocity (kpc/Myr) and acceleration (kpc/Myr^2)
	*/
	double r, vc, denom, fr;


	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);  // kpc

	denom = K[8];
	for (int i = 7; i >= 0; i--)
		denom = denom * r + K[i];  // (km/s)^-1

	vc = KM_PER_SEC_TO_KPC_PER_MYR / denom;  // kpc/Myr
	fr = vc * vc / r;  // kpc/Myr^2

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	dxdt[3] = -x[0] / r * fr;
	dxdt[4] = -x[1] / r * fr;
	dxdt[5] = -x[2] / r * fr;
}
