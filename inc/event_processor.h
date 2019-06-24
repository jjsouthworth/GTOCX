#pragma once

#include "events.h"

class Event_Processor {
public:
	Event_Processor() = default;

	// push event onto heap
	void push(const Event& e);

	// evolve system by moving all events within dt of next event into buffer
	void evolve(double dt);

	// get access to buffer (FIFO)
	const std::deque<Event>& buffer() const;

	// pop all first event from buffer (and all concurrent events)
	std::vector<Event> pop();

private:
	struct comparator {
		bool operator()(const Event& lhs, const Event& rhs){return lhs.time < rhs.time;};
		bool operator()(const Event& lhs, double rhs){return lhs.time < rhs;};
		bool operator()(double lhs, const Event& rhs){return lhs < rhs.time;};
	};
	std::priority_queue<Event, std::vector<Event>, comparator> heap_;
	std::deque<Event> buffer_;
};
