#include <iostream>
#include "star.h"

int main() {
    Galaxy allStars;
    allStars.loadStars("/Users/joshandkris/projects/gtocx/data/stars.txt");
    std::cout << "Loaded: " << allStars.starVec.size() << " stars." << std::endl;
    StateVec starState;
    for (uint ii=0; ii<allStars.starVec.size(); ++ii) {
        if (ii % 10000 == 0) {
            starState = allStars.starVec[ii].getState(90.0);
            std::cout << "x: " << starState.x << " y: " << starState.y << \
                " z: " << starState.z << std::endl;
            std::cout << "vx: " << starState.vx << " vy: " << starState.vy << \
                " vz: " << starState.vz << std::endl;
            //printStar(allStars.starVec[ii]);
        }
    }

    return 0;
}
