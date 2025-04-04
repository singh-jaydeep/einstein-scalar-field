#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "init.hpp"
#include "constants.hpp"

using namespace Const;

using namespace std;
using dvector = vector<double>;

void initialize_xi(dvector &field_xi,const double &amplitude){
    for (int i =0; i <= nR-1; i++){
        field_xi[i] = amplitude*gausAmpl*exp(-pow((i*dR - gausMean),2)/ gausWidth) - gausAmpl*exp(-pow(gausMean,2)/ gausWidth);
    };
};

void initialize_pi(dvector &field_xi, dvector &field_pi){
    for (int i =0; i <= nR-1; i++){
        field_pi[i] = field_xi[nR-1];
    };

};