#include <iostream>
#include <cassert>

#include "dkmeta.h"
#include "dk2nu.h"

#include "TMath.h"
#include "TVector3.h"

// forward declarations
int CalcEnuWgt(const dk2nu* dk2nuObj, const TVector3& xyz,
               double& enu, double& wgt_xy);

// user interface
void calcLocationWeights(dkmeta* dkmetaObj, dk2nu* dk2nuObj);
