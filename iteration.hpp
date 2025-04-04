#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
using dvector = vector<double>;

void eomstep_xi(dvector &field_xi_next, const dvector &field_xi,const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi);
void eomstep_pi(dvector &field_pi_next, const dvector &field_xi,const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi);
void iteration_loop(dvector &field_xi_next, dvector &field_pi_next, dvector &field_xi, dvector &field_pi, dvector &alpha, dvector &beta, dvector &psi, dvector &mass, dvector &expansion, double &max_curv_scalar, bool &supercrit);