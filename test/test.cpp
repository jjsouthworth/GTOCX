#include <iostream>
#include "star.h"

int main() {
    StarCollection allStars;
    allStars.loadStars("/Users/joshandkris/projects/gtocx/data/stars.txt");
    std::cout << "Loaded: " << allStars.starVec.size() << " stars." << std::endl;
    for (uint ii=0; ii<allStars.starVec.size(); ++ii) {
        if (ii % 10000 == 0) {
            printStar(allStars.starVec[ii]);
        }
    }

    return 0;
}
