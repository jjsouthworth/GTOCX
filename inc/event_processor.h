#pragma once

#include "shiplog.h"

class Event_Processor {
public:
	Event_Processor() = default;

	// push ShipLog/event onto heap
	void push(const ShipLog& e);

	// evolve system by moving all ShipLogs/events within dt of next ShipLog/event into buffer
	void evolve(double dt);

	// get access to buffer (FIFO)
	const std::deque<ShipLog>& buffer() const;

	// pop all first ShipLog/event from buffer (and all concurrent ShipLogs/events)
	std::vector<ShipLog> pop();

	// number of ShipLogs/events
	size_t size() const;

private:
	struct comparator {
		bool operator()(const ShipLog& lhs, const ShipLog& rhs){return lhs.start_time() < rhs.start_time();};
		bool operator()(const ShipLog& lhs, double rhs){return lhs.start_time() < rhs;};
		bool operator()(double lhs, const ShipLog& rhs){return lhs < rhs.start_time();};
	};
	std::priority_queue<ShipLog, std::vector<ShipLog>, comparator> heap_;
	std::deque<ShipLog> buffer_;
};
