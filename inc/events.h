#pragma once

#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>
#include <queue>
#include <utility>

#include "common.h"
#include "transfer.h"
#include "star.h"

class ShipLog {
    public:
      // Constructor
      ShipLog(double start_time, Star& starting_star, Transfer& transfer);

      // Member Functions
      void print_log();
      void add_colonized_star(double time, Star& star);
      std::vector<Star> get_colonized_stars();

    protected:
      std::vector<std::pair<double, Star>> colonized_stars;
      Star starting_star;
      double start_time;
      Transfer transfer;
};
