#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "elliptics.hpp"
#include "grid_functions.hpp"
#include "constants.hpp"
#include "csv_management.hpp"
#include "init.hpp"
#include "iteration.hpp"

using namespace Const;
using namespace std;
using dvector = vector<double>;




int main(){

  // declare the field vector, and initialize the geometry. 
  dvector field_xi(nR,0); dvector field_pi(nR,0); dvector field_xi_next(nR,0); dvector field_pi_next(nR,0);
  dvector alpha(nR,1); dvector beta(nR,0.1); dvector psi(nR,1);
  dvector mass(nR,0); dvector mass_aspect(nR,0);
  dvector expansion(nR,0); double max_curv_scalar=0;

  bool supercrit = false;

  // initialize the data and geometry

  initialize_xi(field_xi,amplitude); initialize_pi(field_xi, field_pi);
  solve_elliptics(field_xi,field_pi,alpha,beta,psi);
  compute_mass(mass, alpha,beta,psi);
  compute_expansion(expansion,alpha,beta,psi);

  //initialize storage
  initializeCSV(field_xi,field_pi,alpha,beta,psi, expansion, mass);

  //main iteration loop
  iteration_loop(field_xi_next, field_pi_next, field_xi, field_pi, alpha, beta, psi, mass, expansion, max_curv_scalar, supercrit);
  record_parameters_CSV(amplitude, supercrit, max_curv_scalar);

  cout << "well you're all done now! what did you expect?" << "\n";
  if(supercrit){
    cout << "this run was supercritical, and there was a tiny implosion" << "\n";
  } else {
    cout << "this run was subcritical, and achieved a maximum curvature scalar value of " << max_curv_scalar << "\n";
  };
  return 0;
}











