#pragma once

#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>
#include <queue>

#include "common.h"

enum event_type{
    LAUNCH,
    LAND,
    MANUVER,
    RENDEVOUS,
    SETTLE,
    ORIGIN,
    NONE,
};

class Event {
    public:

        // Variables
        int id;
        double time;
        StateVec state;
        event_type type;

        // Member Functions
        Event(int,event_type,double,StateVec); //No custom desctructor required(R3)
        void print_event();
};

class Event_Queue {
public:
	Event_Queue() = default;

	void push(const Event& e);
	std::vector<Event> pop(double dt);

private:
	struct comparator {
		bool operator()(const Event& lhs, const Event& rhs){
			return lhs.time < rhs.time;
		};
	};
	std::priority_queue<Event, std::vector<Event>, comparator> heap_;
};

class Event_List{
    public:
        Event_List();
        std::vector<Event> events;
};
