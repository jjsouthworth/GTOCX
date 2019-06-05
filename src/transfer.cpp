#include "common.h"
#include "star.h"
#include "vessels.h"
#include "propagate.h"
#include <nlopt.hpp>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>


using namespace std;
typedef std::vector<double> vec_type;
typedef Eigen::Matrix<double, 6, 6> STM;
typedef Eigen::Matrix<double, 3, 3> SubSTM;

int eom_shooter(vec_type *x1, vec_type *x2, double dt, vec_type *dv1, vec_type *dv2, double tol = 1e-13){
  /*Compute the two delta-v's for transfering from R1 to R2 with
  transfer time dt.*/

  //Initialize matrices.
  int i, iter = 0;
  int max_iter = 100;
  double err;
  STM Phi;
  SubSTM B;
  vec_type x;

  vec_type &x1v = *x1;
  vec_type &x2v = *x2;

  // Initialize position and velocity error vectors.
  Eigen::Vector3d dR2 = Eigen::Vector3d::Zero(3);
  Eigen::Vector3d dV1 = Eigen::Vector3d::Zero(3);

  // Add an STM to the end of the state vector.
  vec_type xi = { x1v[0], x1v[1], x1v[2], x1v[3], x1v[4], x1v[5],
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 1.0
                  };

  // Add rectilinear approximation for initial delta-v.
  for (i=0; i<3; i++)
    xi[i+3] += (x2v[i] - x1v[i]) / dt;

  // Iterate up to some max number of iterations. Break early if converged.
  for (iter=0; iter < max_iter; iter++){
    // Copy the initial state and add the previous iteration's delta-v.
    x.assign(xi.begin(), xi.end());

    // Propagate the state forward by dt.
    propagate_eom_stm(&x, 0.0, dt);

    // Compute the final position error.
    for (i=0; i<3; i++)
      dR2[i] = x2v[i] - x[i];

    // Check if the position error is within toleratnce. If so, break the loop.
    err = sqrt(dR2[0]*dR2[0] + dR2[1]*dR2[1] + dR2[2]*dR2[2]);
    if (err < tol)
      break;

    // Convert the trailing 36 elements into a matrix and invert it.
    for (int i=0; i<36; i++)
      Phi(i) = x[i+6];

    // Extract the 3x3 submatrix and compute the velocity correction.
    B = Phi.block<3,3>(0,3);
    dV1 = B.colPivHouseholderQr().solve(dR2);

    // Compute the corrections to initial velocity and apply it to the initial conditions.
    for (i=0; i<3; i++)
      xi[i+3] += dV1[i];
  }

  // Compute and store the initial and final delta-v.
  for (i=0; i<3; i++){
    (*dv1)[i] = xi[i+3] - x1v[i+3];
    (*dv2)[i] = x2v[i+3] - x[i+3];
  }

  // Return 0 if we converged. Return 1 if not.
  if (iter == max_iter - 1)
    return 1;
  else
    return 0;
}


void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, vec_type *dv1, vec_type *dv2){
  vec_type x1_vec = x1->as_posvel_vector();
  vec_type x2_vec = x2->as_posvel_vector();
  eom_shooter(&x1_vec, &x2_vec, t2 - t1, dv1, dv2);
}

void two_impulse_transfer(vec_type *x1, double t1, vec_type *x2, double t2, vec_type *dv1, vec_type *dv2){
    eom_shooter(x1, x2, t2 - t1, dv1, dv2);
}

void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2){
    StateVec x1_sv = star1->getState(t1);
    StateVec x2_sv = star2->getState(t1);
    two_impulse_transfer(&x1_sv, t1, &x2_sv, t2, dv1, dv2);
}

// A struct that contains extra data needed to compute the objective function.
typedef struct {
  double ti, tf;
  vec_type *xi, *xf;
} transfer_data;


double deltav_cost(unsigned n, const double *x, double *grad, void *data){
  /*The objective function, which returns the total delta-v for a 5-impulse
   transfer.*/
  int i;
  vec_type dv4(3, 0.0), dv5(3, 0.0);
  double cost;

  vec_type dv1 = { x[0], x[1], x[2] };
  vec_type dv2 = { x[3], x[4], x[5] };
  vec_type dv3 = { x[6], x[7], x[8] };
  double t2 = x[9];
  double t3 = x[10];
  double t4 = x[11];

  // Type cast transfer data to it's actual type.
  transfer_data *d = (transfer_data *) data;

  // Propagate the first leg of the trajectory.
  vec_type x_vec((*d->xi));  // Create a vector with initial state.
  for (i=0; i<3; i++)
    x_vec[i+3] += dv1[i];
  propagate_eom(&x_vec, d->ti, t2);

  // Repeat for the second and third legs.
  for (i=0; i<3; i++)
    x_vec[i+3] += dv2[i];
  propagate_eom(&x_vec, t2, t3);

  for (i=0; i<3; i++)
    x_vec[i+3] += dv3[i];
  propagate_eom(&x_vec, t3, t4);

  // Now build a two-impulse transfer from time of fourth maneuver to final position.
  two_impulse_transfer(&x_vec, t4, d->xf, d->tf, &dv4, &dv5);

  // TODO: Compute Gradient!!!

  vector<vec_type> all_dvs{ dv1, dv2, dv3, dv4, dv5 };
  for(auto dv : all_dvs)
    cost += dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];

  return cost;

}

void build_settler_ship_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2){
  /*Build an optimal transfer between two stars for a Settler ship while
  meeting all of it's constraints.*/

  // Create an instance of a Settler Ship so we have access to the constraints.
  // StateVec sv;
  // SettlerShip ship = SettlerShip(sv);

  // Create a two-impulse transfer from Star 1 to Star 2. Ignore constraints.
  two_impulse_transfer(star1, t1, star2, t2, dv1, dv2);


}
