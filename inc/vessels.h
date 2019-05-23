#pragma once

#include <vector>

#include "common.h"
#include "star.h"

class Vessel {
public:
    Vessel() = delete;
    Vessel(const StateVec& sv, int imps, double msi, double mtdV);

    virtual ~Vessel() = 0; // force Vessel to be abstract

    virtual bool rendezvous(Star* star);

    protected:
      StateVec sv_;
      double dV_;
      int impulses_;

      std::vector<Star*> stars_;

      const double max_single_impulse_;
      const double max_total_dV_;
};

// forward declare all vessel types
class MotherShip;
class SettlementPod;
class FastShip;
class SettlerShip;

// declarations
class MotherShip : public Vessel {
public:
    MotherShip() = delete;
    MotherShip(const StateVec& sv);

    bool has_pods() const;
    SettlementPod drop_pod();

private:
    int pods_;
};

class SettlementPod : public Vessel {
public:
    SettlementPod() = delete;
    SettlementPod(const StateVec& sv);
};

class FastShip : public Vessel {
public:
    FastShip() = delete;
    FastShip(const StateVec& sv);
};

class SettlerShip : public Vessel {
public:
    SettlerShip() = delete;
    SettlerShip(const StateVec& sv);
};
