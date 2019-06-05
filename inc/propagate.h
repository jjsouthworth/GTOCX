#pragma once
#include <vector>

double _kr(double r);

void eqns_of_motion(const std::vector<double> &x, std::vector<double> &dxdt, double t);
void eqns_of_motion_stm(const std::vector<double> &x, std::vector<double> &dxdt, double t);

size_t propagate_eom(std::vector<double> *x, double t0, double t1);
size_t propagate_eom_stm(std::vector<double> *x, double t0, double t1);
