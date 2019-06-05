#include "transfer.h"
#include "star.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

double norm(vector<double> x){
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

int main(){
  // Initialize variables.
  string star_file = "../data/stars.txt";
  double t0 = 0.0;
  double t1 = 10.0;
  double dv_total = 0.0;
  vector<double> dv1(3, 0.0), dv2(3, 0.0);

  // Load all stars.
  Galaxy galaxy;
  Star sol, star;
  galaxy.loadStars(star_file);
  cout << "Loaded: " << galaxy.starVec.size() << " stars." << endl;

  // Build transfers from Sol to all other stars.
  cout <<"\nBuilding transfers from t=" << t0 << " Myrs to t=" << t1 << " MYrs." << endl;
  sol = galaxy.starVec[0];
  for (uint i=1; i < galaxy.starVec.size(); ++i){
    star = galaxy.starVec[i];
    two_impulse_transfer(&sol, 0.0, &star, 10.0, &dv1, &dv2);
    dv_total = norm(dv1) + norm(dv2);
    printf("  Sol -> Star %d:\n dV = %f kpc/Myr\n dV1 = %.6f, %.6f, %.6f kpc/Myr \
              \n dV2 = %.6f, %.6f, %.6f kpc/Myr\n\n",
           i, dv_total, dv1[0], dv1[1], dv1[2], dv2[0], dv2[1], dv2[2]);
  }

  return 0;
}
