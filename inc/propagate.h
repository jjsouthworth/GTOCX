#pragma once
#include <vector>

void eqns_of_motion(const std::vector<double> &x, std::vector<double> &dxdt, double t);
void eqns_of_motion_stm(const std::vector<double> &x, std::vector<double> &dxdt, double t);
