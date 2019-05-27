#include "events.h"
#include "star.h"

Event::Event(int ID, event_type Type, double Time, StateVec State){
    id = ID;
    type = Type;     
    state = State;
    time = Time;
}

void Event::print_event(){
    std::string stype;
    switch(type){
        case LAUNCH:   
            stype = "Launch"; 
            break;
        case LAND:     
            stype = "Land";     
            break;
        case MANUVER:  
            stype = "Manuver";   
            break;
        case RENDEVOUS:
            stype = "Rendevous"; 
            break;
        case SETTLE:   
            stype = "Settle";    
            break;
        case ORIGIN:   
            stype = "Origin";    
            break;
        default:       
            stype = "None";      
            break;
    }
    printf("Event id:\t%d\t%s\n",id,stype.c_str());
    printf("Time:\t\t%f\n",time);
    printf("Location:\t%f\t%f\t%f\n",state.x, state.y, state.z);
    printf("Motion:\t\t%f\t%f\t%f\n",state.vx,state.vy,state.vz);
    printf("Acceleration:\t%f\t%f\t%f\n",state.ax,state.ay,state.az);
}

Event_List::Event_List(){
    // TODO, initialize this to Sol coordinates at this time
    Star sol = Star(0,8.34,180.0,0.0,0.0,-162.47249203851774016,1);
    StateVec state = sol.getState(0.0);
    Event Sol (0,ORIGIN,0,state);
    events.push_back(Sol);
}
