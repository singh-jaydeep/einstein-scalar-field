#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
using dvector = vector<double>;

void initializeCSV(const dvector &field_xi, const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi, const dvector &expansion,const dvector &mass);
void updateCSV(const dvector &field_xi, const dvector &field_pi,const dvector &alpha, const dvector &beta, const dvector &psi, const dvector &expansion,const dvector &mass);
void record_parameters_CSV(const double &amplitude, const double &supercrit, const double &max_curv_scalar);
