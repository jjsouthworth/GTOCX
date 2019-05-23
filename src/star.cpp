#include "star.h"

// Star Class methods
Star::Star() {
    id = -1;
    R = -1.0;
    i = -1.0;
    Omega = -1.0;
    phi = -1.0;
    theta_f = -1.0;
    isSettled = false;
}

Star::Star(int id, double R, double i, double Omega, double phi,\
        double theta_f, bool isSettled=false) {
    this->id = id;
    this->R = R;
    this->i = i;
    this->Omega = Omega;
    this->phi = phi;
    this->theta_f = theta_f;
    this->isSettled = isSettled;
}

Star::~Star() {}

StateVec Star::getState(double t) {
    StateVec starState;
    // intermediate & repeated values
    double vc = 1 / (K8*pow(R,8) + K7*pow(R,7) + K6*pow(R,6) + K5*pow(R,5)\
        + K4*pow(R,4) + K3*pow(R,3) + K2*pow(R,2) + K1*R + K0); // (km/s)
    vc *= ((SEC_PER_YR * YR_PER_MYR) / KM_PER_KPC); // kpc/Myr
    double n = vc / R; // Myr^-1
    double nt_phi = n*t + phi;
    double c_nt_phi = cos(nt_phi);
    double s_nt_phi = sin(nt_phi);

    // populate state
    starState.x = R * (c_nt_phi*cos(Omega) - \
        s_nt_phi*cos(i)*sin(Omega)); // kpc
    starState.y = R * (c_nt_phi*sin(Omega) + \
        s_nt_phi*cos(i)*cos(Omega)); // kpc
    starState.z = R * (s_nt_phi*sin(i)); //kpc
    starState.vx = vc * (-s_nt_phi*cos(Omega) - \
        c_nt_phi*cos(i)*sin(Omega)); // kpc/Myr
    starState.vy = vc * (-s_nt_phi*sin(Omega) + \
        c_nt_phi*cos(i)*cos(Omega)); // kpc/Myr
    starState.vz = vc * (c_nt_phi*sin(i)); //kpc/Myr

    return starState;
}

// Galaxy Class methods
Galaxy::Galaxy() {
    settled = std::vector<bool> (NUM_STARS, false);
    settled[0] = true;
}
Galaxy::~Galaxy() {}

void Galaxy::loadStars(std::string filename) {
    // vars
    std::ifstream infile;
    std::string line;
    Star tempStar;
    int id;
    double R;
    double i;
    double Omega;
    double phi;
    double theta_f;
    bool isSettled = false;

    // read files
    infile.open(filename);
    if (!infile.is_open()) {
        std::cout << "Could not open file" << std::endl;
        return;
    }
    std::getline(infile, line); // skip header
    while (std::getline(infile, line)) {
        std::istringstream values(line);
        std::string value;
        std::getline(values, value, ','); // get star id
        id = atoi(value.c_str());
        std::getline(values, value, ',');
        R = std::stod(value);
        std::getline(values, value, ',');
        i = std::stod(value) * (PI/180.0);
        std::getline(values, value, ',');
        Omega = std::stod(value) * (PI/180.0);
        std::getline(values, value, ',');
        phi = std::stod(value) * (PI/180.0);
        std::getline(values, value, ',');
        theta_f = std::stod(value) * (PI/180.0);
	isSettled = (id == 0) ? true : false;
	tempStar = Star(id, R, i, Omega, phi, theta_f, isSettled);
        starVec.push_back(tempStar);
    }

    infile.close();
}

// non-member functions
void printStar(Star star) {
    std::cout << "Star ID: " << star.id << std::endl;
    std::cout << "Radius: " << star.R << std::endl;
    std::cout << "Inclination: " << star.i << std::endl;
    std::cout << "Omega: " << star.Omega << std::endl;
    std::cout << "Phi: " << star.phi << std::endl;
    std::cout << "Final Theta: " << star.theta_f << std::endl;
    std::cout << "Settled: " << star.isSettled << std::endl;
    std::cout << std::endl << std::endl;
}
