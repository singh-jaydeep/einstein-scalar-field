#include "grid_functions.hpp"
#include "constants.hpp"

using namespace Const;


void compute_dr(dvector &output, const dvector &field){
    int size = field.size();
    for (int i = 1; i < size-1; i++){
        output[i] = 1/(2*dR) *(field[i+1]-field[i-1]);
    };
    output[0] = 1/(2*dR)*(-3*field[0] + 4*field[1] - field[2]); // second order accurate, one-sided finite differences
    output[size-1] = 1/(2*dR)*(3*field[size-1] - 4*field[size-2] + field[size-3]);
};

void compute_ddr(dvector &output, const dvector &field){
    int size = field.size();
    for (int i = 1; i < size-1; i++){
        output[i] = 1/(pow(dR,2)) *(field[i+1]-2*field[i]+field[i-1]);
    };
    output[0] = 1/(pow(dR,2))*(2*field[0]-5*field[1]+4*field[2]-field[3]);
    output[size-1] = 1/(pow(dR,2))*(2*field[size-1]-5*field[size-2]+4*field[size-3]-field[size-4]);
};

double magnitude(dvector &input){
  int len = input.size();
  double store = 0;
  for (int i =0; i<= len-1; i++){
      store += pow(input[i],2);
  }
  return sqrt(store);
};

void add_vect(dvector &vec1, dvector &vec2){ // important: overrides vec1
    int len = vec1.size();
    assert(len = vec2.size());
    for (int i =0; i <= len-1; i++){
      vec1[i] = vec1[i] + vec2[i];
    };
};

void print_vec(dvector &vec){
   int len = vec.size();
   for (int i =0; i<= len-1; i++){
     cout << vec[i] << "\n";
   };
};

double max_value(const dvector &vec){
	double max = 0;
	int len = vec.size();
	for (int i =0; i <= len-1; i++){
		if (abs(vec[i]) > max){
			max = abs(vec[i]);
		};
	};
	return max;
};

bool contains_negative(const dvector &vec, double &index){
	int len = vec.size();
	for (int i =0; i <= len-1; i++){
		if( vec[i] < 0  ){
			index = i;
			return true;
		};
	};
	return false;
};

void reset_vec(dvector &vec, const double val){
	int len = vec.size();
	for (int i =0; i<= len-1; i++){
		vec[i] = val;
	};
};
