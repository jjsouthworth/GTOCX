#include <memory>

#include "event_processor.h"
#include "vessels.h"
#include "star.h"

std::vector<Event> process_settle(Event& e) {
	std::vector<Event> result;

	// placeholder...we need a way to get vessel and star directly from event being proccesed
	std::shared_ptr<Vessel> vessel;
	std::shared_ptr<Star> star;

	std::vector<Event> events = vessel->settle(star);
	result.insert(result.end(), events.begin(), events.end());

	return result;
}

std::vector<Event> Event_Processor::process(Event& e) {
	switch(e.type) {
		case LAUNCH:
			break;
		case LAND:
			break;
		case MANUVER:
			break;
		case RENDEVOUS:
			break;
		case SETTLE:
			return process_settle(e);
		case ORIGIN:
			break;
		default:
			break;
	}
}

