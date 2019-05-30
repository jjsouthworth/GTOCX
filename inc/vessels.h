#pragma once

#include <vector>
#include <memory>

#include "common.h"
#include "events.h"
#include "star.h"

class Vessel {
public:
    Vessel() = delete;
    Vessel(const StateVec& sv, int imps, double msi, double mtdV);

    virtual ~Vessel() = default;

    virtual std::vector<Event> settle(std::shared_ptr<Star> star) = 0;

    protected:
      StateVec sv_;
      double dV_;
      int impulses_;

      std::vector< std::shared_ptr<Star> > stars_;

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
    std::vector<Event> settle(std::shared_ptr<Star> star);

private:
    int pods_;

    Event drop_pod();
};

class SettlementPod : public Vessel {
public:
    SettlementPod() = delete;
    SettlementPod(const StateVec& sv);

    std::vector<Event> settle(std::shared_ptr<Star> star);
};

class FastShip : public Vessel {
public:
    FastShip() = delete;
    FastShip(const StateVec& sv);

    std::vector<Event> settle(std::shared_ptr<Star> star);
};

class SettlerShip : public Vessel {
public:
    SettlerShip() = delete;
    SettlerShip(const StateVec& sv);

    std::vector<Event> settle(std::shared_ptr<Star> star);
};
