#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

class Star {
public:
    // Constructors & Destructors
    Star();
    Star(int id, double R, double i, double Omega, double phi,\
        double theta_f, bool isSettled);
    ~Star();

    // variables
    int id;
    double R;
    double i;
    double Omega;
    double phi;
    double theta_f;
    bool isSettled;

    // member functions
    
};

class StarCollection {
public:
    // Constructors & Destructors
    StarCollection();
    ~StarCollection();

    // variables
    std::vector<Star> starVec;

    // member functions
    void loadStars(std::string filename);
};

// non-member functions
void printStar(Star star);

