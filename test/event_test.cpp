#include "common.h"
#include "events.h"

int main(void){
    StateVec state {0,0,0,0,0,0,0,0,0};
    Event_List events;
    events.events[0].print_event();
    return 0;
}
