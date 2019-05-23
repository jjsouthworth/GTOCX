#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include "common.h"

class Star {
public:
    // Constructors & Destructors
    Star();
    Star(int id, double R, double i, double Omega, double phi,\
        double theta_f, bool isSettled);
    ~Star();

    // variables
    int id;
    double R; // kpc
    double i; // all angles in rad
    double Omega;
    double phi;
    double theta_f;
    bool isSettled;

    // member functions
    StateVec getState(double t);
};

class Galaxy {
public:
    // Constructors & Destructors
    Galaxy();
    ~Galaxy();

    // variables
    std::vector<Star> starVec;
    std::vector<bool> settled;

    // member functions
    void loadStars(std::string filename);
};

// non-member functions
void printStar(Star star);

