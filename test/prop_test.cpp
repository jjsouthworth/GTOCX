#include "propagate.h"
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

int main() {
	vector<double> x = { 8.34, 0.0, 0.0, 0.0, 0.262773, 0.0 };

	size_t steps = integrate(eqns_of_motion, x, 0.0, 90.0, 0.1);
	printf("N steps: %d\n", (int)steps);
	printf("Final State: %.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n\n", x[0], x[1], x[2], x[3], x[4], x[5]);
 
	vector<double> x1 = {4.17, 7.11292356, 1.25420033, 0.2275681, 0.12939045, 0.02281503,
			     1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

	cout << "x( 0.0) = " << x1[42] << endl;
	steps = integrate(eqns_of_motion_stm, x1, 0.0, 90.0, 0.1);
	cout << "x(90.0) = " << x1[42] << endl;

	return 0;
}

