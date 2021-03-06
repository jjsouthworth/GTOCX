#include "star.h"
#include "propagate.h"


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
    double vc = KM_PER_SEC_TO_KPC_PER_MYR / _kr(R);  // kpc/Myr
    double n = vc / R; // Myr^-1
    double nt_phi = n*t + phi;
    double c_nt_phi = cos(nt_phi);
    double s_nt_phi = sin(nt_phi);
    double c_omega = cos(Omega);
    double s_omega = sin(Omega);
    double cos_i = cos(i);
    double sin_i = sin(i);

    // populate state
    starState.set_x(R * (c_nt_phi*c_omega - s_nt_phi*cos_i*s_omega)); // kpc
    starState.set_y(R * (c_nt_phi*s_omega +  s_nt_phi*cos_i*c_omega)); // kpc
    starState.set_z(R * (s_nt_phi*sin_i)); //kpc
    starState.set_vx(vc * (-s_nt_phi*c_omega - c_nt_phi*cos_i*s_omega)); // kpc/Myr
    starState.set_vy(vc * (-s_nt_phi*s_omega + c_nt_phi*cos_i*c_omega)); // kpc/Myr
    starState.set_vz(vc * (c_nt_phi*sin_i));  //kpc/Myr

    return starState;
}

// Galaxy Class methods
Galaxy::Galaxy() {
    settled = std::vector<bool> (NUM_STARS, false);
    settled[0] = true;
}
Galaxy::~Galaxy() {}

Star &Galaxy::operator[](int index) {
  if (index >= (int)starVec.size()){
    std::cout << "Galaxy index out of bounds. Exiting.";
    exit(1);
  }
  return starVec[index];
}

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

std::vector<Star*> Galaxy::nearestKStars(int starInd, double epoch, int k=1) {
// process stuff, then call other overload of function
    StateVec state = starVec[starInd].getState(epoch);
    std::vector<Star*> neighbors = nearestKStars(state, epoch, k);
    return neighbors;
}

std::vector<Star*> Galaxy::nearestKStars(StateVec state, double epoch, int k=1) {
    std::vector<Star*> neighbors;
    std::map<double, int> distances; //keys=distances, values=indices
    for (uint ii=0; ii<starVec.size(); ++ii) {
        StateVec tmpState = starVec[ii].getState(epoch);
        double mag = state.norm(tmpState);

	// track distances
	distances[mag] = ii;
    }

    // add stars into neighbors vector if norm is one of k-smallest
    // skip any vector that is essentially co-located or settled
    //std::sort(distances.begin(), distances.end());
    std::map<double,int>::iterator index = distances.begin();
    while (neighbors.size() < k && index != distances.end()) {
        if (index->first < 1E-6) {index++; continue;}
        if (starVec[index->second].isSettled) {index++; continue;}
        neighbors.push_back(&starVec[index->second]);
        index++;
	std::cout << "starID: " << index->second << " distance from sol: " << index->first << std::endl;
    }
    return neighbors;
}

// end Galaxy class
//////////////////////////////////////////////////////////////////

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
