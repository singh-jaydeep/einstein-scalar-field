#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
using dvector = vector<double>;

void solve_elliptics(const dvector &field_xi, const dvector &field_pi, dvector &alpha, dvector &beta, dvector &psi);
void solve_tridiag(dvector &output, const dvector &diag1,const dvector &diag2,const dvector &diag3, const dvector &rhs);
void eqn_alpha(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha,
                                                            const dvector &beta, const dvector &psi);
void eqn_psi(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha,
                                                            const dvector &beta, const dvector &psi);
void eqn_beta(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha,
                                                            const dvector &beta, const dvector &psi);
void compute_mass(dvector &output, const dvector &alpha, const dvector &beta, const dvector &psi);
void compute_expansion(dvector &output,const dvector &alpha, const dvector &beta, const dvector &psi);
void compute_curvscalar(dvector &output, const dvector &field_xi, const dvector &field_pi, const dvector &psi);

