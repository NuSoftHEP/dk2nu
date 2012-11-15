#include <iostream>
#include <cassert>

namespace bsim {
  class Decay;
  class Dk2Nu;
  class DkMeta;
}

#include "TMath.h"
#include "TVector3.h"

/// bsim namespace for beam simulation classes and functions
namespace bsim { 

  /// workhorse routine
  int calcEnuWgt(const bsim::Decay& decay, const TVector3& xyz,
                 double& enu, double& wgt_xy);

  /// alternate interface
  int calcEnuWgt(const bsim::Dk2Nu* dk2nu, const TVector3& xyz,
                 double& enu, double& wgt_xy);

  /// user interface
  void calcLocationWeights(const bsim::DkMeta* dkmeta, bsim::Dk2Nu* dk2nu);

} // end-of-namespace "bsim"
