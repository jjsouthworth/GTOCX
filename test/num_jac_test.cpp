#include "common.h"
#include "transfer.h"
#include "star.h"
#include "propagate.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int main(){
  // Initialize variables.
  string star_file = "../data/stars.txt";
  double t1 = 0.0;
  double t2 = 15.0;
  int i1, i2;
  uint i, j;

  // Set the random number generator seed to the current time.
  srand(time(NULL));

  // Load all stars.
  Galaxy galaxy;
  Star star1, star2;
  galaxy.loadStars(star_file);
  cout << "Loaded: " << galaxy.starVec.size() << " stars." << endl;

  // Build transfers from Sol to all other stars.
  cout <<"\nBuilding transfers from t=" << t1 << " Myrs to t=" << t2 << " MYrs." << endl;
  // Create a random transfer.
  // i1 = rand() % (NUM_STARS - 1) + 1;
  // i2 = rand() % (NUM_STARS - 1) + 1;
  i1 = 0;
  i2 = 62245;
  star1 = galaxy[i1];
  star2 = galaxy[i2];
  printf("====================================\n");
  printf("Star %d -> Star %d\n", i1, i2);
  printf("====================================\n");

  vec_type ic(12, 0.0), ic_perturbed(12, 0.0);
  vec_type analytic_grad(12), numerical_grad(12), temp(12);

  const vec_type xi = star1.getState(t1).as_posvel_vector();
  const vec_type xf = star2.getState(t2).as_posvel_vector();
  vec_type deriv(192);

  // Define the transfer data needed by the optimizer.
  transfer_data data = { t1, t2, &xi, &xf, SETTLER_SHIP };

  // Create initial conditions. Start with a two-impulse transfer to get the initial delta-v.
  vec_type dv1(3, 0.0), dv2(3, 0.0);
  ::two_impulse_transfer(&star1, t1, &star2, t2, &dv1, &dv2);

  // Set times for the interior maneuvers at evenly spaced intervals.
  double delta_t = t2 - t1;
  ic[0] = t1 + delta_t*0.25;
  ic[4] = t1 + delta_t*0.5;
  ic[8] = t1 + delta_t*0.75;

  // Propagate the trajectory to get the positions at the interior points.
  vec_type x = xi;
  for (uint i=0; i<3; i++)
    x[i+3] += dv1[i];
  // Propagate to each interior point and store the position.
  propagate_eom(&x, t1, ic[0]);
  ic[1] = x[0];
  ic[2] = x[1];
  ic[3] = x[2];

  propagate_eom(&x, ic[0], ic[4]);
  ic[5] = x[0];
  ic[6] = x[1];
  ic[7] = x[2];

  propagate_eom(&x, ic[4], ic[8]);
  ic[9] = x[0];
  ic[10] = x[1];
  ic[11] = x[2];

  // Compute the numeric cost gradient.
  double cost1 = nlopt_cost(ic, analytic_grad, &data);
  double delta_x = 0.01;

  for (uint i=0; i<12; i++){
  // Perturb the initial conditions for finding the numeric gradient.
    ic_perturbed = ic;
    ic_perturbed[i] += delta_x;

    // Compute new cost and gradient.
    double cost2 = nlopt_cost(ic_perturbed, temp, &data);
    numerical_grad[i] = (cost2 - cost1) / delta_x;
  }

  // Now compare the numerical and analytic gradient.
  printf("\nCOST\nAnalytic - Numeric\n");
  for (i=0; i<12; i++){
    printf("%f - %f = %f\n", analytic_grad[i], numerical_grad[i], analytic_grad[i] - numerical_grad[i]);
  }

  // Repeat for constraints.
  double constraints1[9], constraints2[9];
  double icc[12];
  double ana_grad[84], num_grad[84], tmp[84];

  // Compute analytic graident.
  for (j=0; j<12; j++)
    icc[j] = ic[j];
  nlopt_constraints(7, constraints1, 12, icc, ana_grad, &data);

  for (i=0; i<12; i++){
    for (j=0; j<12; j++)
      icc[j] = ic[j];
    icc[i] += delta_x;

    nlopt_constraints(7, constraints2, 12, icc, tmp, &data);
    for (j=0; j<7; j++){
      num_grad[j*12+i] = (constraints2[j] - constraints1[j]) / delta_x;
    }
  }

  // Now compare the numerical and analytic gradient.
  printf("\nCONSTRAINTS\nAnalytic - Numeric\n");
  for (i=0; i<84; i++){
    printf("%d: %f - %f = %f\n", i, ana_grad[i], num_grad[i], ana_grad[i] - num_grad[i]);
  }

  return 0;
}
