#include <utility>
#include <vector>

#include "events.h"
#include "star.h"
#include "transfer.h"

using namespace std;

ShipLog::ShipLog(double start_time, Star& starting_star, Transfer& transfer):
  start_time(start_time),starting_star(starting_star),transfer(transfer){}

void ShipLog::print_log(){
  cout << "Implment me..." << endl;
}

void ShipLog::add_colonized_star(double time, Star& star){
  pair<double, Star> colonize_event;
  colonize_event = make_pair(time, star);
  colonized_stars.push_back(colonize_event);
}

vector<Star> ShipLog::get_colonized_stars(void){
  vector<Star> stars(colonized_stars.size());
  for (auto colonize_event : colonized_stars)
    stars.push_back(colonize_event.second);

  return stars;
}
