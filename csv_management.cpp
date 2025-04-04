#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <filesystem>

using namespace std;
using dvector = vector<double>;

#include "csv_management.hpp"
#include "constants.hpp"


using namespace Const;

// Verify the directory is correct in constants.hpp

void initializeCSV(const dvector &field_xi, const dvector &field_pi,const dvector &alpha, const dvector &beta, const dvector &psi, const dvector &expansion, const dvector &mass){


      ofstream file(direc + "/EKG_simulation_xi.csv");
      for(int i = 0; i <= nR-1; i++){
          file<< field_xi[i] <<",";
      };
      file << "\n";
      file.close();

      ofstream file1(direc +"/EKG_simulation_pi.csv");
      for(int i = 0; i <= nR-1; i++){
          file1<< field_pi[i] <<",";
      };
      file1 << "\n";
      file1.close();

      ofstream file2(direc +"/EKG_simulation_alpha.csv");
      for(int i = 0; i <= nR-1; i++){
           file2<< alpha[i] <<",";
      };
      file2 << "\n";
      file2.close();

      ofstream file3(direc +"/EKG_simulation_beta.csv");
      for(int i = 0; i <= nR-1; i++){
           file3<< beta[i] <<",";
      };
      file3 << "\n";
      file3.close();

      ofstream file4(direc +"/EKG_simulation_psi.csv");
      for(int i = 0; i <= nR-1; i++){
          file4<< psi[i] <<",";
      };
      file4 << "\n";
      file4.close();

      ofstream file5(direc +"/EKG_simulation_expansion.csv");
      for(int i = 0; i <= nR-1; i++){
          file5<< expansion[i] <<",";
      };
      file5 << "\n";
      file5.close();

      ofstream file6(direc +"/EKG_simulation_mass.csv");
            for(int i = 0; i <= nR-1; i++){
                file6<< mass[i] <<",";
            };
      file6 << "\n";
      file6.close();
};


void updateCSV(const dvector &field_xi, const dvector &field_pi,const dvector &alpha, const dvector &beta, const dvector &psi, const dvector &expansion,const dvector &mass){

	ofstream file(direc +"/EKG_simulation_xi.csv", ios_base::app);
	      for(int i = 0; i <= nR-1; i++){
	          file<< field_xi[i] <<",";
	      };
	      file << "\n";
	      file.close();

	      ofstream file1(direc +"/EKG_simulation_pi.csv", ios_base::app);
	            for(int i = 0; i <= nR-1; i++){
	                file1<< field_pi[i] <<",";
	            };
	            file1 << "\n";
	            file1.close();

	      ofstream file2(direc +"/EKG_simulation_alpha.csv", ios_base::app);
	            for(int i = 0; i <= nR-1; i++){
	                 file2<< alpha[i] <<",";
	            };
	            file2 << "\n";
	            file2.close();

	      ofstream file3(direc +"/EKG_simulation_beta.csv", ios_base::app);
	            for(int i = 0; i <= nR-1; i++){
	                 file3<< beta[i] <<",";
	            };
	            file3 << "\n";
	            file3.close();

	      ofstream file4(direc +"/EKG_simulation_psi.csv", ios_base::app);
	            for(int i = 0; i <= nR-1; i++){
	                file4<< psi[i] <<",";
	            };
	            file4 << "\n";
	            file4.close();
	     ofstream file5(direc +"/EKG_simulation_expansion.csv", ios_base::app);
	         	for(int i = 0; i <= nR-1; i++){
	            	file5<< expansion[i] <<",";
	         	};
	         	file5 << "\n";
	            file5.close();
	     ofstream file6(direc +"/EKG_simulation_mass.csv", ios_base::app);
	            for(int i = 0; i <= nR-1; i++){
	            	 file6<< mass[i] <<",";
	            };
	            file6 << "\n";
	            file6.close();


};

void record_parameters_CSV(const double &amplitude, const double &supercrit, const double &max_curv_scalar){

	 ofstream file(direc +"/EKG_simulation_parameters.csv");
	      file << "R range " << "T range " << "R step size " << "T step size " << "Initial data amplitude " << "Supercriticality " << "Max curvature scalar " << "recordInterval " <<"\n";
	      file << rangeR << ", " << rangeT << ", " << dR << ", " << dT << ", " << amplitude << ", " << supercrit << ", " << max_curv_scalar << recordInterval << "\n";
	      file << "\n";
	      file.close();
};




