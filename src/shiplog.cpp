#include <utility>
#include <vector>

#include "../inc/shiplog.h"
#include "star.h"
#include "transfer.h"

using namespace std;

ShipLog::ShipLog(double start_time, const ship_info& dummy):
  start_time_(start_time), transfer_(dummy){}

ShipLog::ShipLog(double start_time, Star& starting_star, Transfer& transfer):
  start_time_(start_time),starting_star_(starting_star),transfer_(transfer){}

void ShipLog::print_log(){
  cout << "Implment me..." << endl;
}

void ShipLog::add_colonized_star(double time, Star& star){
  pair<double, Star> colonize_event;
  colonize_event = make_pair(time, star);
  colonized_stars_.push_back(colonize_event);
}

vector<Star> ShipLog::get_colonized_stars(void){
  vector<Star> stars(colonized_stars_.size());
  for (auto colonize_event : colonized_stars_)
    stars.push_back(colonize_event.second);

  return stars;
}

double ShipLog::start_time() const {
    return start_time_;
}
