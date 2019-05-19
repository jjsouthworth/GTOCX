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

// StarCollection Class methods
StarCollection::StarCollection() {}
StarCollection::~StarCollection() {}

void StarCollection::loadStars(std::string filename) {
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
        i = std::stod(value);
        std::getline(values, value, ',');
        Omega = std::stod(value);
        std::getline(values, value, ',');
        phi = std::stod(value);
        std::getline(values, value, ',');
        theta_f = std::stod(value);
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
