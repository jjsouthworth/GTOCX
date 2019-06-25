#include <memory>
#include <algorithm>

#include "event_processor.h"
#include "star.h"

void Event_Processor::push(const Event& e) {
	heap_.push(e);
}

void Event_Processor::evolve(double dt) {
	if(heap_.empty()) return;

	// if buffer is empty, grab next element from heap
	if(buffer_.empty()) {
		buffer_.push_back(heap_.top());
		heap_.pop();
	}

	// grab event time at front of buffer
	double t0 = buffer_.front().time;

	// find all events in heap that occur within dt of t0
	while(!heap_.empty()) {
		if(heap_.top().time - t0 < dt) {
			buffer_.push_back(heap_.top());
			heap_.pop();
		} else {
			break;
		}
	}
}

const std::deque<Event>& Event_Processor::buffer() const {
	return buffer_;
}

std::vector<Event> Event_Processor::pop() {
	if(buffer_.empty()) return std::vector<Event>();

	auto bounds = std::equal_range(buffer_.begin(), buffer_.end(),
			               buffer_.front().time, comparator());

	std::vector<Event> result = std::vector<Event>(bounds.first, bounds.second);
	buffer_.erase(bounds.first, bounds.second);

	return result;
}
