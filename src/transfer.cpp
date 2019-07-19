#include "common.h"
#include "transfer.h"
#include "star.h"
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

// Define functions that are used to build transfers.
int eom_shooter(vec_type *x1, vec_type *x2, double dt, vec_type *dv1, vec_type *dv2,
  double tol, vec_type *deriv){
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
    // Copy the initial state which has the previous iteration's delta-v added.
    x = xi;

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

  // Use the STM to construct the derivative of the velocities with respect to
  // the initial and final positions and times.
  if (deriv != NULL){
    // Declare new variables needed.
    vec_type &dd = *deriv;
    SubSTM A, C, D, Binv;
    Eigen::Vector3d dv1_dt1, dv1_dt2, dv2_dt1, dv2_dt2;
    SubSTM dv1_dr1, dv1_dr2, dv2_dr1, dv2_dr2;
    Eigen::Vector3d v1 = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d v2 = Eigen::Vector3d::Zero(3);
    Eigen::Vector3d a2 = Eigen::Vector3d::Zero(3);
    vec_type dxdt(6);

    // Store the initial and final velocities.
    for (i=0; i<3; i++){
      v1[i] = xi[i+3];
      v2[i] = x[i+3];
    }
    // Get the final acceleration as well.
    eqns_of_motion(x, dxdt);
    for (i=0; i<3; i++)
      a2[i] = dxdt[i+3];

    // Extract sub-matrices from the STM.
    A = Phi.block<3,3>(0, 0);
    B = Phi.block<3,3>(0,3);
    C = Phi.block<3,3>(3, 0);
    D = Phi.block<3,3>(3, 3);
    Binv = B.inverse(); // Eigen docuemntation says this is best for small (<= 4x4) matrices.

    // Compute derivatives.
    dv1_dt1 = Binv*v2;
    dv1_dr1 = (-Binv*A);
    dv1_dt2 = -Binv*v2;
    dv1_dr2 = Binv;

    dv2_dt1 = D*Binv*v2 - a2;
    dv2_dr1 = (C - D*Binv*A);
    dv2_dt2 = a2 - D*Binv*v2;
    dv2_dr2 = (D*Binv);

    // Resize the derivative vector and store derivative components.
    dd.resize(48);
    for (i=0; i<3; i++)
      dd[i] = dv1_dt1(i);
    for (i=0; i<9; i++)
      dd[i+3] = dv1_dr1(i);
    for (i=0; i<3; i++)
      dd[i+12] = dv1_dt2(i);
    for (i=0; i<9; i++)
      dd[i+15] = dv1_dr2(i);

    for (i=0; i<3; i++)
      dd[i+24] = dv2_dt1(i);
    for (i=0; i<9; i++)
      dd[i+27] = dv2_dr1(i);
    for (i=0; i<3; i++)
      dd[i+36] = dv2_dt2(i);
    for (i=0; i<9; i++)
      dd[i+39] = dv2_dr2(i);
  }

  // Return 0 if we converged. Return 1 if not.
  if (iter == max_iter - 1)
    return 1;
  else
    return 0;
}


void two_impulse_transfer(vec_type *x1, double t1, vec_type *x2, double t2, vec_type *dv1, vec_type *dv2){
    eom_shooter(x1, x2, t2 - t1, dv1, dv2);
}

void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, DeltaV *dv1, DeltaV *dv2){
  vec_type x1_vec = x1->as_posvel_vector();
  vec_type x2_vec = x2->as_posvel_vector();
  vec_type dv1_vec = dv1->get_vector();
  vec_type dv2_vec = dv2->get_vector();
  // Run the shooter, then re-assign discovered values.
  eom_shooter(&x1_vec, &x2_vec, t2 - t1, &dv1_vec, &dv2_vec);
  dv1->set_from_vector(dv1_vec);
  dv1->set_time(t1);
  dv2->set_from_vector(dv2_vec);
  dv2->set_time(t2);
}

void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, vec_type *dv1, vec_type *dv2){
  vec_type x1_vec = x1->as_posvel_vector();
  vec_type x2_vec = x2->as_posvel_vector();
  // Run the shooter, then re-assign discovered values.
  eom_shooter(&x1_vec, &x2_vec, t2 - t1, dv1, dv2);
}

void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2){
  StateVec x1_sv = star1->getState(t1);
  StateVec x2_sv = star2->getState(t2);
  two_impulse_transfer(&x1_sv, t1, &x2_sv, t2, dv1, dv2);
}

void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, DeltaV *dv1, DeltaV *dv2){
    StateVec x1_sv = star1->getState(t1);
    StateVec x2_sv = star2->getState(t2);
    vec_type dv1_vec = dv1->get_vector();
    vec_type dv2_vec = dv2->get_vector();

    // Compute transfer and reassign variables.
    two_impulse_transfer(&x1_sv, t1, &x2_sv, t2, &dv1_vec, &dv2_vec);
    dv1->set_from_vector(dv1_vec);
    dv1->set_time(t1);
    dv2->set_from_vector(dv2_vec);
    dv2->set_time(t2);
}

const vector<vec_type> nlopt_build_transfer(const vec_type *state_vec, transfer_data *data, vec_type *deriv){
  uint i;
  double tol = 1e-10;
  vec_type dv1, dv2, dv3, dv4, dv5;
  vec_type dv1_tmp(3, 0.0), dv2_tmp(3, 0.0);
  const vec_type& x = *state_vec;
  vec_type shooter_deriv;

  double t2 = x[0];
  vec_type R2 = { x[1], x[2], x[3] };
  double t3 = x[4];
  vec_type R3 = { x[5], x[6], x[7] };
  double t4 = x[8];
  vec_type R4 = { x[9], x[10], x[11] };

  // Ensure the derivative vector is the correct size.
  vec_type &global_deriv = *deriv;
  global_deriv.resize(144);

  // Propagate the first leg of the trajectory.
  double t1 = data->ti;
  vec_type x1 = (*data->xi);  // Create a vector with initial state.
  vec_type x2 = { R2[0], R2[1], R2[2], 0.0, 0.0, 0.0 }; // Initialize the velocity as 0, then use the solved "delta-v" as the velocity solution.
  ::eom_shooter(&x1, &x2, t2 - t1, &dv1_tmp, &dv2_tmp, tol, &shooter_deriv);
  // Store the first delta-v.
  dv1 = dv1_tmp;
  // Replace the solved second "delta-v" as the velocity of arrival at R2.
  for (i=0; i<3; i++)
    x2[i+3] = -dv2_tmp[i];
  // Store the derivative from this leg of the transfer in the total derivative vector.
  for (i=0; i<12; i++){
    global_deriv[i] = shooter_deriv[i+12];
    global_deriv[i+12] = shooter_deriv[i+36];
  }

  // Repeat for the second leg of the trajectory.
  vec_type x3 = { R3[0], R3[1], R3[2], 0.0, 0.0, 0.0 };
  ::eom_shooter(&x2, &x3, t3 - t2, &dv1_tmp, &dv2_tmp, tol, &shooter_deriv);
  // Store the second delta-v.
  dv2 = dv1_tmp;
  // Replace the solved second "delta-v" as the velocity of arrival at R3.
  for (i=0; i<3; i++)
    x3[i+3] = -dv2_tmp[i];
  // Store the derivative from this leg of the transfer in the total derivative vector.
  for (i=0; i<48; i++)
    global_deriv[i+24] = shooter_deriv[i];

  // Repeat for the third leg of the transfer.
  vec_type x4 = { R4[0], R4[1], R4[2], 0.0, 0.0, 0.0 };
  ::eom_shooter(&x3, &x4, t4 - t3, &dv1_tmp, &dv2_tmp, tol, &shooter_deriv);
  // Store the third delta-v.
  dv3 = dv1_tmp;
  // Replace the solved second "delta-v" as the velocity of arrival at R4.
  for (i=0; i<3; i++)
    x4[i+3] = -dv2_tmp[i];
  // Store the derivative from this leg of the transfer in the total derivative vector.
  for (i=0; i<48; i++)
    global_deriv[i+72] = shooter_deriv[i];

  // Repeat for the final leg of the transfer.
  double t5 = data->tf;
  vec_type x5 = (*data->xf);  // Create a vector with initial state.
  ::eom_shooter(&x4, &x5, t5 - t4, &dv1_tmp, &dv2_tmp, tol, &shooter_deriv);
  // Store the fourth and fifth delta-v.
  dv4 = dv1_tmp;
  dv5 = dv2_tmp;
  // Store the derivative from this leg of the transfer in the total derivative vector.
  for (i=0; i<12; i++){
    global_deriv[i+120] = shooter_deriv[i];
    global_deriv[i+132] = shooter_deriv[i+24];
  }

  const vector<vec_type> all_dvs = { dv1, dv2, dv3, dv4, dv5 };
  return all_dvs;
}

double nlopt_cost(const vec_type &x, vec_type &dcdx, void *data){
  /*The objective function, which returns the total delta-v for a 5-impulse
   transfer.*/
  uint i, j, k;
  double cost = 0.0;
  vec_type dvdx(144);

  // Type cast transfer data to it's actual type.
  transfer_data *d = (transfer_data *) data;

  // Construct the delta-v's from the state vector.
  const vector<vec_type> dvs = nlopt_build_transfer(&x, d, &dvdx);

  // Sum the cost.
  for(auto dv : dvs)
    cost += dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];

  // Resize the derivative if it is not the correct size.
  if (!dcdx.empty()){
    dcdx.resize(12);
    for (i=0; i<12; i++)
      dcdx[i] = 0.0;

    for (i=0; i<3; i++){
      for (j=0; j<4; j++){
        k=3*j + i; // locator witin a block of dv#/dx partials.
        dcdx[j] += dvs[0][i]*dvdx[k] + dvs[1][i]*(dvdx[24+k] - dvdx[12+k]) - dvs[2][i]*dvdx[48+k];
        dcdx[j+4] += dvs[1][i]*dvdx[k+36] + dvs[2][i]*(dvdx[k+72] - dvdx[k+60]) - dvs[3][i]*dvdx[96+k];
        dcdx[j+8] += dvs[2][i]*dvdx[k+84] + dvs[3][i]*(dvdx[k+120] - dvdx[k+108]) - dvs[4][i]*dvdx[132+k];
      }
    }
    // Multiply all partials by 2.0 since that was skipped in the accumulation step.
    for (i=0; i<12; i++)
      dcdx[i] *= 2.0;
  }
  return cost;
}


void nlopt_constraints(unsigned, double *result, unsigned, const double *x, double *deriv, void* data){
  /*The vector inequality constraint function. Enforces time between
    maneuvers and max delta-v per impulse. */
  uint i, j, k;
  double dv_sq, msi_sq;
  vec_type dvdx(144);  // velocity partials

  // Type cast transfer data to it's actual type.
  transfer_data *d = (transfer_data *) data;
  msi_sq = pow(d->ship.max_single_impulse, 2);

  // Construct the transfer to get the delta-v's.
  const vec_type state(x, x+12);
  const vector<vec_type> dvs = nlopt_build_transfer(&state, d, &dvdx);

  // The first five constraints are that delta-v's do not exceed the max.
  for (uint i=0; i<5; i++){
    dv_sq = pow(dvs[i][0], 2.0) + pow(dvs[i][1], 2.0) + pow(dvs[i][2], 2.0);
    result[i] = dv_sq - msi_sq;
  }

  // Ensure that the transfer times are sequential and spaced by 1 Myrs.
  result[5] = x[0] - x[4] + MYRS_BETWEEN_DELTAV;
  result[6] = x[4] - x[8] + MYRS_BETWEEN_DELTAV;

  // Initialize the constraint jacobian.
  if (deriv){
    vector<vec_type> dcdx(7);  // constraint partials
    for (i=0; i<7; i++)
      dcdx[i] = vec_type(12, 0.0);

    // Compute the jacobian of the delta-v constraints.
    for (i=0; i<3; i++){
      for (j=0; j<4; j++){
        k = 3*j + i;
        dcdx[0][j] += 2.0*dvs[0][i]*dvdx[k];

        dcdx[1][j] += 2.0*dvs[1][i]*(dvdx[k+24] - dvdx[k+12]);
        dcdx[1][j+4] += 2.0*dvs[1][i]*dvdx[k+36];

        dcdx[2][j] -= 2.0*dvs[2][i]*dvdx[k+48];
        dcdx[2][j+4] += 2.0*dvs[2][i]*(dvdx[k+72] = dvdx[k+60]);
        dcdx[2][j+8] += 2.0*dvs[2][i]*dvdx[k+84];

        dcdx[3][j+4] -= 2.0*dvs[3][i]*dvdx[k+96];
        dcdx[3][j+8] += 2.0*dvs[3][i]*(dvdx[k+120] - dvdx[k+108]);

        dcdx[4][j+8] -= 2.0*dvs[4][i]*dvdx[k+132];
      }
    }

    // Compute the jacobian for the time constraints.
    dcdx[5][0] = 1.0;
    dcdx[5][4] = -1.0;
    dcdx[6][4] = 1.0;
    dcdx[6][8] = -1.0;

    // Populate the array with the values computed.
    for (i=0; i<7; i++){
      for (j=0; j<12; j++)
         deriv[i*12 + j] = dcdx[i][j];
    }
  }
  return;
}

// Define the DeltaV class's member functions.
DeltaV::DeltaV(double t, double dvx, double dvy, double dvz) {
  this->t = t;
  this->dvx = dvx;
  this->dvy = dvy;
  this->dvz = dvz;
}

double& DeltaV::operator[](const int index) {
  if (index == 0) { return dvx; }
  if (index == 1) { return dvy; }
  if (index == 2) { return dvz; }
  throw "DeltaV objects only have 3 indices.";
}

void DeltaV::set_from_vector(vec_type dv_vec) {
  if (dv_vec.size() != 3)
    throw "Vector to set delta-v is not 3 elements long.";

  this->dvx = dv_vec[0];
  this->dvy = dv_vec[1];
  this->dvz = dv_vec[2];
}

void DeltaV::set_values(double t, double dvx, double dvy, double dvz) {
  this->t = t;
  this->dvx = dvx;
  this->dvy = dvy;
  this->dvz = dvz;
}

void DeltaV::set_time(double t) { this-> t = t; }
vec_type DeltaV::get_vector() {
  vec_type x = { dvx, dvy, dvz };
  return x;
}
double DeltaV::get_time() { return t; }
double DeltaV::mag() { return sqrt(dvx*dvx + dvy*dvy + dvz*dvz);  }

// Define the Transfer class's member functions.
Transfer::Transfer(const ship_info &ship_type):ship(ship_type){
  this->deltavs = vector<DeltaV>(ship.num_impulses, DeltaV());
}

double Transfer::total_deltav() {
  double cost = 0.0;
  for (auto dv : deltavs)
    cost += dv.mag();
  return cost;
}

bool Transfer::_check_common_constraints() {
  if (!constructed) { return false; } // Check that transfer was constructed.

  // Check that we don't exceed the max number of maneuvers.
  if (deltavs.size() > (uint)ship.num_impulses) { return false; }

  // Check that no delta-v's exceed the max impulse or falls outside the simulation..
  for (auto dv : deltavs) {
    if (dv.mag() > ship.max_single_impulse) { return false; }
    if (dv.get_time() > 90.0) { return false; }
    if (dv.get_time() < 0.0) { return false; }
  }

  // Check that the total delta-v does not exceed the max.
  if (this->total_deltav() > ship.max_deltav) { return false; }

  //Check that maneuvers are spaced by 1 Myrs.
  for (uint i=0; i < deltavs.size(); i++)
    for (uint j=0;  j < deltavs.size(); j++) {
      if (i >= j || deltavs[i].mag() < 1e-15 || deltavs[j].mag() < 1e-15)
        continue;
      if ( abs(deltavs[i].get_time() - deltavs[j].get_time()) < MYRS_BETWEEN_DELTAV)
        return false;
    }
  return true;  // If nothing trips false, maneuvers are valid.
}

void Transfer::reset_transfer() {
  for (auto dv : deltavs)
    dv.set_values();
  constructed = false;
}


// Define the Settler Ship class's member functions.
SettlerShipTransfer::SettlerShipTransfer(Star &star1, Star &star2)
: Transfer(SETTLER_SHIP){
this->star1 = &star1;
this->star2 = &star2;
}

void SettlerShipTransfer::print_transfer() {
  uint i = 1;
  for (DeltaV dv : deltavs){
    printf("  t%d = %f    DV%d: %f, %f, %f kpc/Myr\n", \
           i, dv.get_time(), i, dv[0], dv[1], dv[2]);
    i++;
  }
  printf("Delta-V Total: %f\t(max = %f)\n", total_deltav(), ship.max_deltav);

  if (!valid_transfer())
    printf("INVALID TRANSFER!\n");
}

bool SettlerShipTransfer::valid_transfer() { return this->_check_common_constraints(); }

void SettlerShipTransfer::two_impulse_transfer(double t1, double t2) {
  // Reset any data from currently stored transfers.
  reset_transfer();

  DeltaV& dv1 = deltavs[0];
  DeltaV& dv2 = deltavs.back();

  ::two_impulse_transfer(star1, t1, star2, t2, &dv1, &dv2);

  constructed = true;
}

void SettlerShipTransfer::optimal_two_impulse_transfer(double t1, double max_transfer_time) {
  // Reset any data from currently stored transfers.
  reset_transfer();

  constructed = true;
}

void SettlerShipTransfer::optimal_five_impulse_transfer(double t1, double t2) {
  /*Build an optimal transfer between two stars for a Settler ship while
  meeting all of it's constraints.*/
  vec_type ic(12);
  double optimal_func_value;
  const vec_type xi = (*star1).getState(t1).as_posvel_vector();
  const vec_type xf = (*star2).getState(t2).as_posvel_vector();
  vec_type deriv(192);

  // Reset any data from currently stored transfers.
  reset_transfer();

  // Define the transfer data needed by the optimizer.
  transfer_data data = { t1, t2, &xi, &xf, this->ship };

  // Generate the optimizer.
  // The state vector will be [t2, R2, t3, R3, t4, R4] (12-elements)
  auto optimizer = nlopt::opt(nlopt::LD_SLSQP , 12);
  optimizer.set_min_objective(::nlopt_cost, &data);

  // Define upper and lower bounds for all variables.
  vec_type lb(12, -OUTER_RADIUS_BOUNDARY);
  lb[0] = t1 + MYRS_BETWEEN_DELTAV;
  lb[4] = t1 + MYRS_BETWEEN_DELTAV;
  lb[8] = t1 + MYRS_BETWEEN_DELTAV;
  optimizer.set_lower_bounds(lb);

  vec_type ub(12, OUTER_RADIUS_BOUNDARY);
  ub[0] = t2 - MYRS_BETWEEN_DELTAV;
  ub[4] = t2 - MYRS_BETWEEN_DELTAV;
  ub[8] = t2 - MYRS_BETWEEN_DELTAV;
  optimizer.set_upper_bounds(ub);

  // Define inequality constraint function.
  const vec_type constraint_tol(7, 1e-6);
  optimizer.add_inequality_mconstraint(::nlopt_constraints, (void*)&data, constraint_tol);

  // Define the initial step size
  vec_type init_step(12, 0.1);
  optimizer.set_initial_step(init_step);

  // Set the stopping criteria.
  optimizer.set_stopval(1e-6);
  optimizer.set_ftol_rel(1e-6);
  optimizer.set_ftol_abs(1e-6);
  optimizer.set_xtol_rel(1e-6);
  optimizer.set_xtol_abs(1e-6);
  // optimizer.set_maxeval(10000);
  // optimizer.set_maxtime(30.0); // seconds

  // Create initial conditions. Start with a two-impulse transfer to get the initial delta-v.
  vec_type dv1(3, 0.0), dv2(3, 0.0);
  ::two_impulse_transfer(star1, t1, star2, t2, &dv1, &dv2);

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

  // Call the optimizer.
  nlopt::result status = optimizer.optimize(ic, optimal_func_value);
  cout << "Status: " << status << endl;

  // Store the result in this class's members.
  const vector<vec_type> solved_mnvrs = ::nlopt_build_transfer(&ic, &data, &deriv);
  deltavs[0].set_time(t1);
  deltavs[1].set_time(ic[0]);
  deltavs[2].set_time(ic[4]);
  deltavs[3].set_time(ic[8]);
  deltavs[4].set_time(t2);

  for (uint i=0; i<5; i++)
    deltavs[i].set_from_vector(solved_mnvrs[i]);

  constructed = true;
}


// Define the Fast Ship class's member functions.
FastShipTransfer::FastShipTransfer(Star &star1, Star &star2) :
Transfer(FAST_SHIP) {
  this->star1 = &star1;
  this->star2 = &star2;
}

bool FastShipTransfer::valid_transfer() { return this->_check_common_constraints(); }

void FastShipTransfer::print_transfer() {
  uint i = 1;
  for (DeltaV dv : deltavs){
    printf("  t%d = %f    DV%d: %f, %f, %f kpc/Myr\n", \
           i, dv.get_time(), i, dv[0], dv[1], dv[2]);
    i++;
  }
  printf("Delta-V Total: %f\t(max = %f)\n", total_deltav(), ship.max_deltav);

  if (!valid_transfer())
    printf("INVALID TRANSFER!\n");
}

void FastShipTransfer::two_impulse_transfer(double t1, double t2) {
  // Reset any data from currently stored transfers.
  reset_transfer();

  DeltaV& dv1 = deltavs[0];
  DeltaV& dv2 = deltavs.back();

  dv1.set_time(t1);
  dv2.set_time(t2);
  ::two_impulse_transfer(star1, t1, star2, t2, &dv1, &dv2);

  constructed = true;
}

void FastShipTransfer::optimal_two_impulse_transfer(double t1, double max_transfer_time) {
  // Reset any data from currently stored transfers.
  reset_transfer();

  constructed = true;
}

// Define the Mother Ship class's member functions.
MotherShipTransfer::MotherShipTransfer() : Transfer(MOTHER_SHIP){
}

void MotherShipTransfer::print_transfer() { printf("Implement me...\n"); }

bool MotherShipTransfer::valid_transfer() {
  if (!this->_check_common_constraints()) { return false; }
  // ADD ADDITIONAL VALIDITY CHECKS.
  return true;
}
