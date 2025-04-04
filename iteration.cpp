#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "grid_functions.hpp"
#include "constants.hpp"
#include "elliptics.hpp"
#include "csv_management.hpp"

using namespace std;
using dvector = vector<double>;

using namespace Const;

void iteration_loop(dvector &field_xi_next, dvector &field_pi_next, dvector &field_xi, dvector &field_pi, dvector &alpha, dvector &beta, dvector &psi, dvector &mass, dvector &expansion,double &max_curv_scalar, bool &supercrit){
	dvector curv_scalar(nR,0); double temp_curv_scalar; double apparent_horizon_index = 0;

	for(int i =1; i <= nT; i++){
		  if (LOUD){
			  cout << "iteration number " << i << "\n";
		  };

	      eomstep_xi(field_xi_next,field_xi,field_pi,alpha,beta,psi);
	      eomstep_pi(field_pi_next, field_xi,field_pi,alpha,beta,psi);
	      field_pi_next[nR-1] = -1*field_xi_next[nR-1]; // enforce the outgoing radiation condition

	      field_pi = field_pi_next;
	      field_xi = field_xi_next;


	      solve_elliptics(field_xi,field_pi,alpha,beta,psi);
	      compute_mass(mass,alpha,beta,psi);
	      compute_curvscalar(curv_scalar, field_xi,field_pi, psi);
	      temp_curv_scalar = max_value(curv_scalar);
	      if(temp_curv_scalar > max_curv_scalar){
	    	 max_curv_scalar = temp_curv_scalar;
	      };

	      compute_expansion(expansion, alpha,beta,psi);

	      if(contains_negative(expansion, apparent_horizon_index)){
	    	  	if (!supercrit){
	    	  		cout << "a trapped surface has formed at radius " << apparent_horizon_index*dR << "\n";
	    	  	};

	    	  	supercrit = true;
	      } else {
	    	  	supercrit = false;
	      };

	      if(supercrit && SINGULARITY_AVOIDING){
	    	  break;
	      };

	      // update CSV if condition is met
	      if (i % recordInterval == 0){
	    	  updateCSV(field_xi, field_pi,alpha,beta,psi, expansion,mass);
	    	  cout << "new max curvature scalar is " << max_curv_scalar << "\n";
	      };


	  };
};



void eomstep_xi(dvector &field_xi_next, const dvector &field_xi,const dvector &field_pi, const dvector &alpha,const dvector &beta, const dvector &psi){
    dvector rhs_preder(nR);
    dvector rhs_der(nR);
    for (int i=0; i <= nR-1; i++){
        rhs_preder[i] = alpha[i]*field_pi[i] / (pow(psi[i],2)) + beta[i]*field_xi[i];
    };
    compute_dr(rhs_der,rhs_preder);

    field_xi_next[0] = 0;
    for(int i =1; i <= nR-2; i++){
        field_xi_next[i] = .5*(field_xi[i-1]+field_xi[i+1]) + dT * rhs_der[i];
    };
    field_xi_next[nR-1] = 2*field_xi_next[nR-2] -field_xi_next[nR-3];
};

void eomstep_pi(dvector &field_pi_next, const dvector &field_xi,const dvector &field_pi, const dvector &alpha,const dvector &beta, const dvector &psi){
    dvector dbeta(nR), dpsi(nR), dxi(nR);
    compute_dr(dbeta,beta); compute_dr(dpsi,psi); compute_dr(dxi,field_xi);
    dvector temp1(nR); dvector dtemp1(nR);

    for (int i =0; i<= nR-1; i++){
        temp1[i] = pow(psi[i],4)*(beta[i]*field_pi[i]+ (1/pow(psi[i],2))*alpha[i]*field_xi[i]);
    };
    compute_dr(dtemp1,temp1);

    dvector rhs_der(nR);
    for (int i =1; i <= axisSplitnR-1; i++){
      rhs_der[i] = (1/pow(psi[i],4))*dtemp1[i] + 2*field_pi[i]*dbeta[i] + (1/pow(psi[i],2))*2*alpha[i]*dxi[i] - 2/3*field_pi[i]*dbeta[i]
                    - 4/3*field_pi[i]*dbeta[i] - 4*field_pi[i]*beta[i]*dpsi[i]/psi[i];
      field_pi_next[i] = .5*(field_pi[i-1]+field_pi[i+1]) + dT * rhs_der[i];
    };

    for (int i=axisSplitnR; i <= nR-2; i++){
      rhs_der[i] = (1/pow(psi[i],4))*dtemp1[i] + 2/(dR*i)*(beta[i]*field_pi[i] + alpha[i]*field_xi[i]*(1/pow(psi[i],2))) - 2/3*field_pi[i]*dbeta[i]
                  -4/3*1/(dR*i)*field_pi[i]*beta[i]  - 4*field_pi[i]*beta[i]*dpsi[i]/psi[i];
      field_pi_next[i] = .5*(field_pi[i-1]+field_pi[i+1])  + dT * rhs_der[i];
    };
    field_pi_next[0] = field_pi_next[2]; field_pi_next[1] = field_pi_next[2];
};