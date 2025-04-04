#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
using dvector = vector<double>;


void compute_dr(dvector &output, const dvector &field);
void compute_ddr(dvector &output, const dvector &field);
double magnitude(dvector &input);
void add_vect(dvector &vec1, dvector &vec2);
void print_vec(dvector &vec);
double max_value(const dvector &vec);
bool contains_negative(const dvector &vec, double &index);
void reset_vec(dvector &vec, const double val);
