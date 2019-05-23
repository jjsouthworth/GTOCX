#include <stdexcept>

#include "vessels.h"

// Vessel
Vessel::Vessel(const StateVec& sv, int imps, double msi, double mtdV) :
    sv_(sv), dV_(0.0), impulses_(imps), stars_(),
    max_single_impulse_(msi), max_total_dV_(mtdV) {}

Vessel::~Vessel() {}

// MotherShip
MotherShip::MotherShip(const StateVec& sv) :
    Vessel(sv, 3, 200*KM_PER_SEC_TO_KPC_PER_MYR,
                  500*KM_PER_SEC_TO_KPC_PER_MYR),
    pods_(10) {}
                      
bool MotherShip::has_pods() const {
    return pods_ > 0;
}

SettlementPod MotherShip::drop_pod() {
    if(!has_pods()) throw std::runtime_error("No remaining SettlemetPods remaining in Mothership");
    --pods_;
    return SettlementPod(this->sv_);
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
