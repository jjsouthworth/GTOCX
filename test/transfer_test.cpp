#include "transfer.h"
#include "star.h"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

int main(){
  // Initialize variables.
  string star_file = "../data/stars.txt";
  double t0 = 0.0;
  double t1 = 10.0;
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
    printf("  Star %d: dV1 = %.6f, %.6f, %.6f, dV2 = %.6f, %.6f, %.6f",
           i, dv1[0], dv1[1], dv1[2], dv2[0], dv2[1], dv2[2]);
  }

  return 0;
}
