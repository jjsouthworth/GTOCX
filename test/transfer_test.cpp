#include "transfer.h"
#include "star.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int main(){
  // Initialize variables.
  string star_file = "../data/stars.txt";
  double t0 = 0.0;
  double t1 = 10.0;
  int i1, i2;

  // Load all stars.
  Galaxy galaxy;
  Star star1, star2;
  galaxy.loadStars(star_file);
  cout << "Loaded: " << galaxy.starVec.size() << " stars." << endl;

  // Build transfers from Sol to all other stars.
  cout <<"\nBuilding transfers from t=" << t0 << " Myrs to t=" << t1 << " MYrs." << endl;
  for (uint i=1; i < 10; ++i){
    // Create a random transfer.
    i1 = rand() % (NUM_STARS - 1) + 1;
    i2 = rand() % (NUM_STARS - 1) + 1;
    star1 = galaxy[i1];
    star2 = galaxy[i2];
    printf("====================================\n");
    printf("Star %d -> Star %d\n", i1, i2);
    printf("====================================\n");

    printf("FAST SHIP TRANSFER\n");
    auto fs_transfer = FastShipTransfer(star1, star2);

    printf("\nEmpty Transfer:\n");
    fs_transfer.print_transfer();

    printf("\nTwo-Impulse Transfer\n");
    cout << "Enter two_impulse_transfer:" << endl;
    fs_transfer.two_impulse_transfer(t0, t1);
    fs_transfer.print_transfer();

    printf("-----------------------------------\n");
    printf("SETTLER SHIP TRANSFERS\n");
    auto ss_transfer = SettlerShipTransfer(star1, star2);

    printf("Empty Transfer:\n");
    ss_transfer.print_transfer();

    printf("\nTwo-Impulse Transfer\n");
    ss_transfer.two_impulse_transfer(t0, t1);
    ss_transfer.print_transfer();

    printf("\nFive-Impulse-Transfer\n");
    ss_transfer.optimal_five_impulse_transfer(t0, t1);
    ss_transfer.print_transfer();
  }

  return 0;
}
