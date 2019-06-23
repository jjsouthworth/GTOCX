#include <iostream>
#include <vector>
#include <numeric>
#include <thread>
#include <fstream>
#include <chrono>

#include "common.h"
#include "transfer.h"
#include "star.h"

int main(int argc, char** argv) {

  std::string star_file = "../data/stars.txt";
  if (argc > 1) {
    star_file = argv[1];
  }

  // Load and initialize star data
  Galaxy galaxy;
  galaxy.loadStars(star_file);
  const int num_stars = galaxy.starVec.size();
  std::cout << "Loaded: " << num_stars << " stars.\n";

  // Set the random number generator seed to the current time.
  Star base_star = galaxy[0];

  // For every other star in the galaxy, check to see if a Settler
  // ship is able to reach it in 90 Myr
  std::vector<int> stars_possible(num_stars, 0);
  stars_possible[0] = 1; // Can get to our own star in zero time
  unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
  std::cout << "Num threads = " << concurentThreadsSupported << std::endl;
  auto start = std::chrono::system_clock::now();
#pragma omp parallel for num_threads(concurentThreadsSupported)
  // Start loop at 1, 0 is sol, the star to test
  for (int ii = 1; ii < 1000/*num_stars*/; ++ii) {
    Star& other_star = galaxy.starVec[ii];
    SettlerShipTransfer sf_transfer(base_star, other_star);
    sf_transfer.two_impulse_transfer(0, 90);
    if (sf_transfer.valid_transfer()) {
      stars_possible[ii] = 1;
      // std::cout << "VALID TRANSFER\n";
    }
    if (ii % 50 == 0) {
      std::cout << "Calculating star " << ii << std::endl;
    }
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Total time = " << elapsed.count() << '\n';

  int num_possible = std::accumulate(stars_possible.begin(), stars_possible.end(), 0);
  std::cout << "Number stars calculated: " << num_stars
            << ", and " << num_possible << " have valid transfers\n";
  std::ofstream ofile;
  ofile.open("possible_stars.txt", std::ios::out);
  for (int possible : stars_possible) {
    ofile << possible << "\n";
  }
  ofile.close();
  return 0;
}
