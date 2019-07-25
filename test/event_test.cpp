#include <vector>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real.hpp>

#include "event_processor.h"
#include "shiplog.h"
#include "common.h"

int main(void){

    boost::random::mt19937 gen;

    // generate some initial events within the first 5 million years (time chosen arbitarily)
    boost::random::uniform_int_distribution<> dist1(1, 6); // up to 6 initial events, at least 1
    size_t n_events = dist1(gen);

    Event_Processor processor;
    boost::random::uniform_real_distribution<> dist2(0.0, 5.0E6); // randomly choose time within first 5M years
    ship_info dummy;
    for(size_t ii = 0; ii < n_events; ++ii) {
        ShipLog log(dist2(gen), dummy);
        processor.push(log);
    }

    // now, pull these events and randomly generate new events based on processed events
    const double dt = 1.0E6;
    std::vector<ShipLog> processed_events;
    boost::random::uniform_real_distribution<> dist3(0, 1.0);
    while(processor.size()) {
        processor.evolve(dt);
        std::vector<ShipLog> events = processor.pop();

        // randomly generate subsequent events
        for(size_t ii = 0; ii < events.size(); ++ii) {
            double t = events[ii].start_time() + 3.0E6; // 3M years from current event start_time
            if(t > 90.0E6) continue; // skip events that occur after 90M years

            if(dist3(gen) < 0.75) { // 75% of processed events will generate new events
                ShipLog log(t, dummy);
                // add 3 new events (settle)
                processor.push(log);
                processor.push(log);
                processor.push(log);
            }
        }

        processed_events.insert(processed_events.end(), events.begin(), events.end());
    }

    std::cout << "n_events: " << n_events << std::endl;
    std::cout << "processed_events: " << processed_events.size() << std::endl;
    for(size_t ii = 0; ii < processed_events.size(); ++ii) {
        std::cout << processed_events[ii].start_time() << std::endl;
    }

    return 0;
}
