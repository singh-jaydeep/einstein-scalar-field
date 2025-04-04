#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
using dvector = vector<double>;

void initialize_xi(dvector &field_xi, const double &amplitude);
void initialize_pi(dvector &field_xi, dvector &field_pi);