

namespace Const{
  const int rangeR = 40;
  const int rangeT = 200;
  const double dR = .03;
  const double dT = .01;
  const int nR = rangeR / dR;
  const int nT = rangeT / dT;

  const double gausAmpl = 0.035;
  const double gausWidth = 5;
  const double gausMean = 20.0;

  double amplitude = .5771; //rescaling to probe critical collapse

  const int recordInterval = 20;

  const int axisSplitR = 0;
  const int axisSplitnR = axisSplitR / dR;

  const double tolerance = .000001;
  const int maxiters = 150;


  const bool LOUD = true;
  const bool SINGULARITY_AVOIDING = true; // if true, simulation will terminate immediately upon locating an apparent horizon


  string direc = "/Users/Singh/Desktop/ekg_c++_updated"; // Update with folder you want to write to
};
