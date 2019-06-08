#pragma once

#include <array>
#include <vector>
#include <iostream>

class StateVec {
 public:
  StateVec () :
      sv_data()
  {}

  StateVec (const std::array<double, 6>& data)
  {
    sv_data[0] = data[0];
    sv_data[1] = data[1];
    sv_data[2] = data[2];
    sv_data[3] = data[3];
    sv_data[4] = data[4];
    sv_data[5] = data[5];
    sv_data[6] = sv_data[7] = sv_data[8] = 0;
  }

  StateVec (double x, double y, double z,
            double vx, double vy, double vz,
            double ax, double ay, double az)
  {
    sv_data[0] = x;
    sv_data[1] = y;
    sv_data[2] = z;
    sv_data[3] = vx;
    sv_data[4] = vy;
    sv_data[5] = vz;
    sv_data[6] = ax;
    sv_data[7] = ay;
    sv_data[8] = az;
  }

  // Getter methods

  inline double x() const {
	return sv_data[0];
  }

  inline double y() const {
	return sv_data[1];
  }

  inline double z() const {
	return sv_data[2];
  }

  inline double vx() const {
	return sv_data[3];
  }

  inline double vy() const {
	return sv_data[4];
  }

  inline double vz() const {
	return sv_data[5];
  }

  inline double ax() const {
	return sv_data[6];
  }

  inline double ay() const {
	return sv_data[7];
  }

  inline double az() const {
	return sv_data[8];
  }

  // Setter methods

  inline void set_x(double x) {
    sv_data[0] = x;
  }

  inline void set_y(double y) {
    sv_data[1] = y;
  }

  inline void set_z(double z) {
    sv_data[2] = z;
  }

  inline void set_vx(double vx) {
    sv_data[3] = vx;
  }

  inline void set_vy(double vy) {
    sv_data[4] = vy;
  }

  inline void set_vz(double vz) {
    sv_data[5] = vz;
  }

  inline void set_ax(double ax) {
    sv_data[6] = ax;
  }

  inline void set_ay(double ay) {
    sv_data[7] = ay;
  }

  inline void set_az(double az) {
    sv_data[8] = az;
  }

  inline std::array<double, 6> as_posvel_array() {
    std::array<double, 6> tmp = {sv_data[0], sv_data[1], sv_data[2],
                                 sv_data[3], sv_data[4], sv_data[5]};
    return tmp;
  }

  inline std::vector<double> as_posvel_vector() {
    std::vector<double> tmp = {sv_data[0], sv_data[1], sv_data[2],
                               sv_data[3], sv_data[4], sv_data[5]};
    return tmp;
  }

  void print_state(const std::ostream&) const {
    printf("%.3f, %.3f, %.3f, %.6f, %.6f, %.6f\n",sv_data[0], sv_data[1], sv_data[2], sv_data[3], sv_data[4], sv_data[5]);
  }

 private:
  std::array<double, 9> sv_data;
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

static const double MYRS_BETWEEN_DELTAV = 1.0;

// Define structures that hold information about each ship type.
typedef struct {
  int num_impulses;
  double max_single_impulse;
  double max_deltav;
  int num_pods;
} ship_info;

const ship_info MOTHER_SHIP = {3, 200.0*KM_PER_SEC_TO_KPC_PER_MYR,500.0*KM_PER_SEC_TO_KPC_PER_MYR, 10};
const ship_info SETTLEMENT_POD = {1, 300.0*KM_PER_SEC_TO_KPC_PER_MYR, 300.0*KM_PER_SEC_TO_KPC_PER_MYR, 0};
const ship_info FAST_SHIP = {2, 1500.0*KM_PER_SEC_TO_KPC_PER_MYR, 1500.0*KM_PER_SEC_TO_KPC_PER_MYR, 0};
const ship_info SETTLER_SHIP = {5, 175.0*KM_PER_SEC_TO_KPC_PER_MYR, 400.0*KM_PER_SEC_TO_KPC_PER_MYR, 0};
