#pragma once
#include <vector>
#include "star.h"

typedef std::vector<double> vec_type;

int eom_shooter(double *x1, double *x2, double dt, vec_type *dv1, vec_type *dv2, double tol = 1e-13);

void build_settler_ship_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2);


void two_impulse_transfer(StateVec *x1, double t1, StateVec *x2, double t2, vec_type *dv1, vec_type *dv2);
void two_impulse_transfer(vec_type *x1, double t1, vec_type *x2, double t2, vec_type *dv1, vec_type *dv2);
void two_impulse_transfer(Star *star1, double t1, Star *star2, double t2, vec_type *dv1, vec_type *dv2);
