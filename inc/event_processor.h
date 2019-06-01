#pragma once

#include "events.h"

class Event_Processor {
public:
	Event_Processor() = default;

	std::vector<Event> process(Event& e);
};
