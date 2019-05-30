#include <stdexcept>
#include <memory>

#include "vessels.h"

// Vessel
Vessel::Vessel(const StateVec& sv, int imps, double msi, double mtdV) :
    sv_(sv), dV_(0.0), impulses_(imps), stars_(),
    max_single_impulse_(msi), max_total_dV_(mtdV) {}

// MotherShip
MotherShip::MotherShip(const StateVec& sv) :
    Vessel(sv, 3, 200*KM_PER_SEC_TO_KPC_PER_MYR,
                  500*KM_PER_SEC_TO_KPC_PER_MYR),
    pods_(10) {}
                      
bool MotherShip::has_pods() const {
    return pods_ > 0;
}

std::vector<Event> MotherShip::settle(std::shared_ptr<Star> star) {
	std::vector<Event> result;

	// generate event to create settlement pod
	Event e1 = drop_pod();
	result.push_back(e1);

	// generate a rendevous event for settlement pod
	Event e2(0,NONE,0.0,sv_); // need an interface to do this...
	result.push_back(e2);

	// generate a settle event for settlement pod
	Event e3(0,NONE,0.0,sv_); // need an interface to do this...
	result.push_back(e3);

    return result;
}

Event MotherShip::drop_pod() {
    if(!has_pods()) {
        throw std::runtime_error("No remaining SettlemetPods remaining in Mothership");
    }

    --pods_;

    return Event(0,NONE,0.0,sv_);
}

// SettlementPod
SettlementPod::SettlementPod(const StateVec& sv) :
        Vessel(sv, 1, 300*KM_PER_SEC_TO_KPC_PER_MYR,
                      300*KM_PER_SEC_TO_KPC_PER_MYR) {}

// FastShip
FastShip::FastShip(const StateVec& sv) :
        Vessel(sv, 2, 1500*KM_PER_SEC_TO_KPC_PER_MYR,
                      1500*KM_PER_SEC_TO_KPC_PER_MYR) {}

// SettlerShip     
SettlerShip::SettlerShip(const StateVec& sv) :
        Vessel(sv, 5, 175*KM_PER_SEC_TO_KPC_PER_MYR,
                      400*KM_PER_SEC_TO_KPC_PER_MYR) {}
