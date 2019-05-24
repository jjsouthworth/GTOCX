#pragma once

struct StateVec {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
};

// Physical Constants
static const double PI = 3.141592653589793;

// Scenario Constants
static const int NUM_STARS = 100001; // includes Sol
static const double K8 = -1.94316e-12; // (km/s)^-1/kpc^8
static const double K7 = 3.7516e-10; // (km/s)^-1/kpc^7
static const double K6 = -2.70559e-8; // (km/s)^-1/kpc^6
static const double K5 = 9.70521e-7; // (km/s)^-1/kpc^5
static const double K4 = -1.88428e-5; // (km/s)^-1/kpc^4
static const double K3 = 0.000198502; // (km/s)^-1/kpc^3
static const double K2 = -0.0010625; // (km/s)^-1/kpc^2
static const double K1 = 0.0023821; // (km/s)^-1/kpc
static const double K0 = 0.00287729; // (km/s)^-1
static const double K[9] = { K0, K1, K2, K3, K4, K5, K6, K7, K8};
 
static const long long KM_PER_KPC = 30856775814671900; // km
static const double YR_PER_MYR = 1e6; // yr
static const int SEC_PER_YR = 31557600; // sec
static const double KM_PER_SEC_TO_KPC_PER_MYR = 1E6*SEC_PER_YR / KM_PER_KPC;
