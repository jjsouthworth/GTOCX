#pragma once
#include <vector>
#include "common.h"
#include "star.h"

typedef std::vector<double> vec_type;

class DeltaV {
public:

  DeltaV(double t=0.0, double dvx=0.0, double dvy=0.0, double dvz=0.0);

  double& operator[] (const int index);

  void set_values(double t=0.0, double dvx=0.0, double dvy=0.0, double dvz=0.0);
  void set_from_vector(vec_type dv_vec);
  void set_time(double t);
  vec_type get_vector();
  double get_time();
  double mag();
protected:
  double t, dvx, dvy, dvz;
};


class Transfer {
public:
  Transfer(const ship_info &ship_type);
  void reset_transfer();
  double total_deltav();
  virtual void print_transfer() { };
  virtual bool valid_transfer() { return false; };
  std::vector<DeltaV> deltavs;
protected:
  bool constructed;
  ship_info ship;
  bool _check_common_constraints();
};


class SettlerShipTransfer: public Transfer {
public:
  SettlerShipTransfer(Star &star1, Star &star2);
  void print_transfer() override;
  bool valid_transfer() override;
  void two_impulse_transfer(double t1, double t2);
  void optimal_two_impulse_transfer(double t1, double max_transfer_time);
  void optimal_five_impulse_transfer(double t1, double t2);
protected:
  Star *star1, *star2;
};

class FastShipTransfer: public Transfer {
public:
  FastShipTransfer(Star &star1, Star &star2);
  void print_transfer() override;
  bool valid_transfer() override;
  void two_impulse_transfer(double t1, double t2);
  void optimal_two_impulse_transfer(double t1, double max_transfer_time);
protected:
  Star *star1, *star2;
};


class MotherShipTransfer: public Transfer {
    MotherShipTransfer();
    void print_transfer() override;
    bool valid_transfer() override;
};


int eom_shooter(double *x1, double *x2, double dt, vec_type *dv1, vec_type *dv2, double tol = 1e-13);

void build_settler_ship_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2);

void two_impulse_transfer(vec_type *x1, double t1, vec_type *x2, double t2, vec_type *dv1, vec_type *dv2);
void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, DeltaV *dv1, DeltaV *dv2);
void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, vec_type *dv1, vec_type *dv2);
void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, DeltaV *dv1, DeltaV *dv2);
void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2);

// A struct that contains extra data needed to compute NLOPT objective functions.
typedef struct {
  const double ti, tf;
  const vec_type *xi, *xf;
  ship_info ship;
} transfer_data;


int eom_shooter(vec_type *x1, vec_type *x2, double dt, vec_type *dv1, vec_type *dv2, double tol=1e-13, vec_type *deriv=NULL);
const std::vector<vec_type> nlopt_build_transfer(const vec_type *state_vec, transfer_data *data, vec_type *deriv);
double nlopt_cost(const vec_type &x, vec_type &dcdx, void *data);
void nlopt_constraints(unsigned, double *result, unsigned, const double *x, double *deriv, void* data);
