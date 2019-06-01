#include "propagate.h"
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace boost::numeric::odeint;

int main() {

	// Define types and steppers for the integrator.
	typedef vector<double> vec_type;
	typedef runge_kutta_fehlberg78 <vec_type> error_stepper;
	controlled_runge_kutta < error_stepper > controlled_stepper;

	vector<double> x = { 8.34, 0.0, 0.0, 0.0, 0.262773, 0.0 };
	printf("\nEQUATIONS OF MOTION:\n");
	printf("  X(t= 0): %.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n", x[0], x[1], x[2], x[3], x[4], x[5]);

	// Propagate the state.
	size_t steps = integrate_adaptive(controlled_stepper, eqns_of_motion, x, 0.0, 90.0, 0.1);

	// Print results.
	printf("  X(t=90): %.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n", x[0], x[1], x[2], x[3], x[4], x[5]);
	printf("  N steps: %d\n", (int)steps);

	// Propagate an STM.
	vector<double> x_stm = {
		4.17, 7.11292356, 1.25420033, 0.2275681, 0.12939045, 0.02281503,
			     1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 1.0
				 };
  printf("\nEQUATIONS OF MOTION WITH STM:\n");
	printf("  X(t= 0): %.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n", x_stm[0], x_stm[1], x_stm[2], x_stm[3], x_stm[4], x_stm[5]);

	// Propagate the trajectory with an STM.
	steps = integrate_adaptive(controlled_stepper, eqns_of_motion_stm, x_stm, 0.0, 90.0, 0.1);

	 // Create a matrix to store the STM for printing purposes.
 	vector<double>::const_iterator first = x_stm.begin() + 6;
 	vector<double>::const_iterator last = x_stm.begin() + 42;
 	vector<double> phi_vec(first, last);
 	Eigen::Map<Eigen::MatrixXd> phi(phi_vec.data(), 6, 6);

	// Print results.
	printf("  X(t=90): %.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n", x_stm[0], x_stm[1], x_stm[2], x_stm[3], x_stm[4], x_stm[5]);
	printf("  N steps: %d\n", (int)steps);
	cout << "\n  STM(t=90):\n" << phi << endl;

	return 0;
}
