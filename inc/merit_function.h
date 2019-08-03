#pragma once

class Solution {
  public:

  	// Constructor
  	Solution();  //TODO: populate this

    // Member functions
    void export_solution_file();
    void parse_solution_file();
    void compute_merit_function();

};

void merit_function(Solution inputs, double merit);
//void f_and_g(double r, double theta, double f_r, double f_theta, double g_r, double g_theta);
//void f(double x, double mu, double s, double f);
