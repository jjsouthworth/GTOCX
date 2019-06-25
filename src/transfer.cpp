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

// A struct that contains extra data needed to compute NLOPT objective functions.
typedef struct {
  const double ti, tf;
  const vec_type *xi, *xf;
  ship_info ship;
} transfer_data;


// Define functions that are used to build transfers.
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
  StateVec x2_sv = star2->getState(t1);
  two_impulse_transfer(&x1_sv, t1, &x2_sv, t2, dv1, dv2);
}

void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, DeltaV *dv1, DeltaV *dv2){
    StateVec x1_sv = star1->getState(t1);
    StateVec x2_sv = star2->getState(t1);
    vec_type dv1_vec = dv1->get_vector();
    vec_type dv2_vec = dv2->get_vector();

    // Compute transfer and reassign variables.
    two_impulse_transfer(&x1_sv, t1, &x2_sv, t2, &dv1_vec, &dv2_vec);
    dv1->set_from_vector(dv1_vec);
    dv1->set_time(t1);
    dv2->set_from_vector(dv2_vec);
    dv2->set_time(t2);
}

const vector<vec_type> _nlopt_build_transfer(const vec_type *state_vec, transfer_data *data){
  uint i;
  vec_type dv4(3, 0.0), dv5(3, 0.0);
  const vec_type& x = *state_vec;

  vec_type dv1 = { x[0], x[1], x[2] };
  vec_type dv2 = { x[3], x[4], x[5] };
  vec_type dv3 = { x[6], x[7], x[8] };
  double t2 = x[9];
  double t3 = x[10];
  double t4 = x[11];

  // Propagate the first leg of the trajectory.
  vec_type x_vec = (*data->xi);  // Create a vector with initial state.
  for (i=0; i<3; i++)
    x_vec[i+3] += dv1[i];
  propagate_eom(&x_vec, data->ti, t2);

  // Repeat for the second and third legs.
  for (i=0; i<3; i++)
    x_vec[i+3] += dv2[i];
  propagate_eom(&x_vec, t2, t3);

  for (i=0; i<3; i++)
    x_vec[i+3] += dv3[i];
  propagate_eom(&x_vec, t3, t4);

  // Now build a two-impulse transfer from time of fourth maneuver to final position.
  vec_type xf = (*data->xf);
  two_impulse_transfer(&x_vec, t4, &xf, data->tf, &dv4, &dv5);

  const vector<vec_type> all_dvs = { dv1, dv2, dv3, dv4, dv5 };
  return all_dvs;
}

double _nlopt_cost(const vec_type &x, vec_type &, void *data){
  /*The objective function, which returns the total delta-v for a 5-impulse
   transfer.*/
  double cost = 0.0;

  // Type cast transfer data to it's actual type.
  transfer_data *d = (transfer_data *) data;

  // Construct the delta-v's from the state vector.
  const vector<vec_type> all_dvs = _nlopt_build_transfer(&x, d);

  // Sum the cost.
  for(auto dv : all_dvs)
    cost += dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];

  return cost;
}


void _nlopt_constraints(unsigned, double *result, unsigned, const double *x, double *, void* data){
  /*The vector inequality constraint function. Enforces time between
    maneuvers and max delta-v per impulse. */

  // Type cast transfer data to it's actual type.
  transfer_data *d = (transfer_data *) data;

  // Construct the transfer to solve for the last two maneuvers.
  const vec_type state(x, x+12);
  const vector<vec_type> solved_mnvrs = _nlopt_build_transfer(&state, d);

  // The first five constraints are that delta-v's do not exceed the max.
  double msi_sq = pow(d->ship.max_single_impulse, 2);
  double dv_sq = 0.0;
  for (uint i=0; i<5; i++){
    dv_sq = pow(solved_mnvrs[i][0], 2.0) + \
      pow(solved_mnvrs[i][1], 2.0) + \
      pow(solved_mnvrs[i][2], 2.0);
    result[i] = dv_sq - msi_sq;
  }

  // Ensure that the transfer times are sequential and spaced by 1 Myrs.
  result[5] = d->ti - x[9] + MYRS_BETWEEN_DELTAV;
  result[6] = x[9] - x[10] + MYRS_BETWEEN_DELTAV;
  result[7] = x[10] - x[11] + MYRS_BETWEEN_DELTAV;
  result[8] = x[11] - d->tf + MYRS_BETWEEN_DELTAV;

  //TODO: Add Gradient.
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

  // Reset any data from currently stored transfers.
  reset_transfer();

  // Define the transfer data needed by the optimizer.
  transfer_data data = { t1, t2, &xi, &xf, this->ship };

  // Generate the optimizer.
  // The state vector will be [dV1, dV2, dV3, t2, t3, t4] (12-elements)
  auto optimizer = nlopt::opt(nlopt::LN_COBYLA , 12);
  optimizer.set_min_objective(::_nlopt_cost, &data);

  // Define upper and lower bounds for all variables.
  vec_type lb(12, -ship.max_single_impulse);
  lb[9] = t1 + MYRS_BETWEEN_DELTAV;
  lb[10] = t1 + MYRS_BETWEEN_DELTAV;
  lb[11] = t1 + MYRS_BETWEEN_DELTAV;
  optimizer.set_lower_bounds(lb);

  vec_type ub(12, ship.max_single_impulse);
  ub[9] = t2 - MYRS_BETWEEN_DELTAV;
  ub[10] = t2 - MYRS_BETWEEN_DELTAV;
  ub[11] = t2 - MYRS_BETWEEN_DELTAV;
  optimizer.set_upper_bounds(ub);

  // Define inequality constraint function.
  vec_type constraint_tol(2, 1e-13);
  optimizer.add_inequality_mconstraint(::_nlopt_constraints, (void*)&data, constraint_tol);

  // Set the stopping criteria.
  optimizer.set_stopval(0.0);
  optimizer.set_ftol_rel(1e-10);
  optimizer.set_ftol_abs(1e-10);
  optimizer.set_xtol_rel(1e-10);
  optimizer.set_xtol_abs(1e-10);
  optimizer.set_maxeval(1000);
  optimizer.set_maxtime(60.0); // seconds

  // Create initial conditions. Start with a two-impulse transfer to set the initial delta-v.
  vec_type dv1(3), dv2(3);
  ::two_impulse_transfer(star1, t1, star2, t2, &dv1, &dv2);
  // Scale the result so it is within the bounds.
  double dv1_mag = sqrt(pow(dv1[0], 2.0) + pow(dv1[1], 2.0) + pow(dv1[2], 2.0));
  double scale = 1.0;
  if (dv1_mag > ship.max_single_impulse)
    scale = ship.max_single_impulse / dv1_mag;

  ic[0] = dv1[0] * scale;
  ic[1] = dv1[1] * scale;
  ic[2] = dv1[2] * scale;

  // Set times for the intermediate maneuvers at evenly spaced intervals.
  double delta_t = t2 - t1;
  ic[9] = t1 + delta_t*0.25;
  ic[10] = t1 + delta_t*0.5;
  ic[11] = t1 + delta_t*0.75;

  // Call the optimizer.
  nlopt::result status = optimizer.optimize(ic, optimal_func_value);

  // Store the result in this class's members.
  const vector<vec_type> solved_mnvrs = ::_nlopt_build_transfer(&ic, &data);
  deltavs[0].set_time(t1);
  deltavs[1].set_time(ic[9]);
  deltavs[2].set_time(ic[10]);
  deltavs[3].set_time(ic[11]);
  deltavs[4].set_time(t2);

  for (uint i; i<5; i++)
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
