
#include "elliptics.hpp"
#include "constants.hpp"
#include "grid_functions.hpp"


using namespace Const;



void solve_elliptics(const dvector &field_xi, const dvector &field_pi, dvector &alpha, dvector &beta, dvector &psi){
	dvector perturb_alpha(nR,0); dvector perturb_psi(nR,0); dvector perturb_beta(nR,0);

	//initialization for iteration loop
	  dvector rhs_psi(nR,0); dvector diag1_psi(nR,0); dvector diag2_psi(nR,0); dvector diag3_psi(nR,0);
	  dvector rhs_beta(nR,0); dvector diag1_beta(nR,0); dvector diag2_beta(nR,0); dvector diag3_beta(nR,0);
	  dvector rhs_alpha(nR,0); dvector diag1_alpha(nR,0); dvector diag2_alpha(nR,0); dvector diag3_alpha(nR,0);


	  eqn_alpha(diag1_alpha,diag2_alpha,diag3_alpha,rhs_alpha,field_xi,field_pi,alpha,beta,psi);
	  eqn_psi(diag1_psi, diag2_psi, diag3_psi, rhs_psi, field_xi, field_pi, alpha, beta, psi);
	  eqn_beta(diag1_beta, diag2_beta,diag3_beta, rhs_beta, field_xi, field_pi, alpha, beta, psi);
	  double resid = magnitude(rhs_alpha) + magnitude(rhs_psi) + magnitude(rhs_beta);



	  for (int i =0; i <= maxiters; i++){
	      if (resid < tolerance){
	        break;
	      };
	      solve_tridiag(perturb_alpha, diag1_alpha,diag2_alpha,diag3_alpha,rhs_alpha);
	      solve_tridiag(perturb_psi, diag1_psi,diag2_psi,diag3_psi,rhs_psi);
	      solve_tridiag(perturb_beta, diag1_beta, diag2_beta, diag3_beta, rhs_beta);



	      add_vect(alpha,perturb_alpha); add_vect(psi,perturb_psi); add_vect(beta, perturb_beta);


	      eqn_alpha(diag1_alpha,diag2_alpha,diag3_alpha,rhs_alpha,field_xi,field_pi,alpha,beta,psi);
	      eqn_psi(diag1_psi, diag2_psi, diag3_psi, rhs_psi, field_xi, field_pi, alpha, beta, psi);
	      eqn_beta(diag1_beta, diag2_beta,diag3_beta, rhs_beta, field_xi, field_pi, alpha, beta, psi);
	      resid = magnitude(rhs_alpha)+ magnitude(rhs_psi) + magnitude(rhs_beta);

	  };
	  if (resid < tolerance && LOUD){
	    cout << "converged with residual " << resid  <<"\n";
	  } else if (resid >= tolerance && LOUD){
	    cout << "failed to converge in allotted steps" << "\n";
	  };

};




void solve_tridiag(dvector &output, const dvector &diag1,const dvector &diag2,const dvector &diag3, const dvector &rhs){
    int size = rhs.size();

    dvector temp1(size); dvector temp2(size);

    temp1[0] = diag3[0] / diag2[0];
    temp2[0] = rhs[0] / diag2[0];

    for (int i=1; i< size; i++) {
    double m = 1.0 / (diag2[i] - diag1[i] * temp1[i-1]);

    if (abs(diag2[i] - diag1[i] * temp1[i-1]) < .00005 ){
      cout << "problem" << "\n";
    };

    temp1[i] = diag3[i] * m;
    temp2[i] = (rhs[i] - diag1[i] * temp2[i-1]) * m;
  };

  output[size-1] = temp2[size-1];
  for (int i= size-1; i-- > 0; ) {
    output[i] = temp2[i] - temp1[i] * output[i+1];
  };
};

void eqn_alpha(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi){
    dvector dalpha(nR,0); dvector ddalpha(nR,0); dvector dbeta(nR,0); dvector dpsi(nR,0);
    compute_dr(dalpha, alpha); compute_ddr(ddalpha,alpha); compute_dr(dbeta,beta); compute_dr(dpsi, psi);

    double m = 0; double a = 0; double b =0; double c = 0;
    double offset = 0; // for testing

    // compute the right hand side for the iteration
    for (int i = 1; i < nR-1; i++){
       m = dbeta[i] - beta[i]/(offset + i*dR);
       out_rhs[i] = -1*( ddalpha[i] + 2*dalpha[i]*(1/(offset + i*dR)  + dpsi[i]/psi[i] ) - 2/(3*alpha[i])*pow(psi[i],4)*  pow(m,2)  -8*3.1415*pow(field_pi[i],2)*alpha[i]);
    };
    out_rhs[0] = alpha[0]-alpha[1];
    out_rhs[nR-1] = 1- alpha[nR-1];

    // compute the diagonals of the jacobian

    for (int i =1; i < nR-1; i++){
        a = 2/(offset + i*dR) + 2*dpsi[i]/psi[i]; b = 2/3*pow(psi[i],4) * pow((dbeta[i] - beta[i]/(offset + i*dR)),2);
        c = 8*3.1415*pow(field_pi[i],2);
        diag1[i] = 1/pow(dR,2) - 1/(2*dR)*a;
        diag2[i] = -2/pow(dR,2) + 1/pow(alpha[i],2)*b - c;
        diag3[i] = 1/pow(dR,2) + 1/(2*dR)*a;
    };
    diag1[0] = 0; diag1[nR-1] = 0;
    diag2[0] = -1; diag2[nR-1] = 1;
    diag3[0] = 1; diag3[nR-1] = 0;
};

void eqn_psi(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi){
    dvector dpsi(nR,0); dvector ddpsi(nR,0); dvector dbeta(nR,0);
    compute_dr(dpsi,psi); compute_ddr(ddpsi, psi); compute_dr(dbeta,beta);

    double m = 0;

    for (int i = 1; i < nR-1; i++){
        m = 1/alpha[i] * (dbeta[i] - beta[i]/(i*dR));
        out_rhs[i] = -1* (ddpsi[i] + dpsi[i]*2/(i*dR) + 1/12* pow(psi[i],5) *pow(m,2)+ 3.1415*psi[i]*(pow(field_xi[i],2) +pow(field_pi[i],2) ));

        diag1[i] = 1/pow(dR,2) - 1/(i*dR*dR);
        diag2[i] = -2/pow(dR,2) + 5/12* pow(psi[i],4)*pow(m,2) + 3.1415*(pow(field_xi[i],2) +pow(field_pi[i],2));
        diag3[i] = 1/pow(dR,2) + 1/(i*dR*dR);
    };
    out_rhs[0] = psi[0]-psi[1];
    out_rhs[nR-1] = 1- psi[nR-1];

    diag1[0] = 0; diag1[nR-1] = 0;
    diag2[0] = -1; diag2[nR-1] = 1;
    diag3[0] = 1; diag3[nR-1] = 0;
};

void eqn_beta(dvector &diag1, dvector &diag2, dvector &diag3, dvector &out_rhs, const dvector &field_xi, const dvector &field_pi, const dvector &alpha, const dvector &beta, const dvector &psi){
    dvector dpsi(nR,0); dvector dalpha(nR,0); dvector dbeta(nR,0); dvector ddbeta(nR,0);
    compute_dr(dpsi,psi); compute_dr(dalpha,alpha); compute_dr(dbeta,beta); compute_ddr(ddbeta,beta);

    double offset = 0.0; // for testing

    double m1 =0; double m2 = 0;

    for (int i = 1; i < nR-1; i++){
        m1 = 2/(offset + i*dR) + 6*dpsi[i]/psi[i] - dalpha[i]/alpha[i];
        m2 = 12*3.1415 * alpha[i] * field_xi[i] * field_pi[i] / pow(psi[i],2);

        out_rhs[i] = -1 *  (ddbeta[i] + (dbeta[i] - beta[i]/(offset + i*dR))*m1 + m2 );
        diag1[i] = 1/(pow(dR,2)) - 1/(2*dR)*m1;
        diag2[i] = -2/(pow(dR,2)) - 1/(offset+i*dR)*m1;
        diag3[i] = 1/(pow(dR,2)) + 1/(2*dR)*m1;
    };
    out_rhs[0] = -1*beta[0]; out_rhs[nR-1] = -1*beta[nR-1];
    diag1[0] = 0; diag1[nR-1] = 0;
    diag2[0] = 1; diag2[nR-1] = 1;
    diag3[0] = 0; diag3[nR-1] = 0;

};


void compute_mass(dvector &output, const dvector &alpha, const dvector &beta, const dvector &psi){
	dvector dbeta(nR,0); dvector dpsi(nR,0); dvector rpsi(nR,0); dvector drpsi(nR,0);
	double m = 0;

	for (int i =0; i <= nR-1; i++){
		rpsi[i] = (i*dR)*psi[i];
	};
	compute_dr(dbeta,beta);
	compute_dr(dpsi,psi);
	compute_dr(drpsi,rpsi);

	for (int i =0; i <= nR-1; i++){
		m = i*dR*dbeta[i] - beta[i];
		output[i] = i*dR*pow(psi[i],6)/(18*pow(alpha[i],2))*pow(m,2)-2*pow(i*dR,2)*dpsi[i]*drpsi[i];
	};
};

void compute_expansion(dvector &output,const dvector &alpha, const dvector &beta, const dvector &psi){
	dvector radius(nR,0); dvector temp1(nR,0); dvector temp2(nR,0); dvector dtemp1(nR,0); dvector dtemp2(nR,0);
	for (int i=1; i <= nR-1; i++){
		radius[i] = i*dR;
		temp1[i] = beta[i]/(i*dR);
		temp2[i] = i*dR*pow(psi[i],2);
	}
	compute_dr(dtemp1,temp1); compute_dr(dtemp2,temp2);
	for (int i =0; i <= nR-1;i++){
		output[i] = 2*i*dR/(3*alpha[i]) * dtemp1[i] +2/(i*dR*pow(psi[i],4)) * dtemp2[i];
	};
	output[0] = 2; output[1] = 2;
};


void compute_curvscalar(dvector &output, const dvector &field_xi, const dvector &field_pi, const dvector &psi){
	for (int i =1; i <= nR-1; i++){
		output[i] = 8*3.1415* (pow(field_xi[i],2) - pow(field_pi[i],2)) / (pow(psi[i],4));
	}
};




