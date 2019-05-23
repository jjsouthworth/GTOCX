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

std::shared_ptr<Vessel> MotherShip::settle(std::shared_ptr<Star> star) {
    std::shared_ptr<SettlementPod> pod = drop_pod();

    if(!pod->rendezvous(star)) {
        throw std::runtime_error("SettlemetPod failed to rendezvous with star");
    }
    pod->settle(star);

    return std::static_pointer_cast<Vessel>(pod);
}

std::shared_ptr<SettlementPod> MotherShip::drop_pod() {
    if(!has_pods()) {
        throw std::runtime_error("No remaining SettlemetPods remaining in Mothership");
    }

    --pods_;

    return std::make_shared<SettlementPod>(this->sv_);
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
