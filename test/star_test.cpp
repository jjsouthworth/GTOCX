#include <iostream>
#include <string>

#include <boost/program_options.hpp>
//#include <Eigen/Dense>

namespace po = boost::program_options;

#include "star.h"

int main(int argc, char** argv) {
  std::string star_file="../data/stars.txt";
  /*
  po::options_description desc("Options in test");
  desc.add_options()
      ("help,h", "Write help message")
      ("starfile,f", po::value(&star_file)->default_value("../data/stars.txt"))
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }
  */

  Galaxy allStars;
  std::vector<Star*> sol_hood;
  allStars.loadStars(star_file);
  std::cout << "Loaded: " << allStars.starVec.size() << " stars." << std::endl;
  StateVec starState;
  for (uint ii=0; ii<allStars.starVec.size(); ++ii) {
    if (ii % 10000 == 0) {
      starState = allStars.starVec[ii].getState(90.0);
      std::cout << "x: " << starState.x() << " y: " << starState.y() << \
          " z: " << starState.z() << std::endl;
      std::cout << "vx: " << starState.vx() << " vy: " << starState.vy() << \
          " vz: " << starState.vz() << std::endl;
      //printStar(allStars.starVec[ii]);
    }
  }
  sol_hood = allStars.nearestKStars(0, 0.0, 5);
  for (uint ii=0; ii<sol_hood.size(); ++ii) {
    printStar(*sol_hood[ii]);
  }
  
  return 0;
}
