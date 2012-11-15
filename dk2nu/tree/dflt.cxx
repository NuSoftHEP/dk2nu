#include "dflt.h"
#include "TMath.h"

Bool_t bsim::IsDefault(Float_t fval)
{
   return (TMath::Abs(fval - kDfltFloat) < 1e-4);
}

Bool_t bsim::IsDefault(Int_t ival) 
{
   return (ival == kDfltInt);
}

Bool_t bsim::IsDefault(UInt_t uival)
{
   return (uival == kDfltUInt);
}

Bool_t bsim::IsDefault(Double_t dval) 
{
   return (TMath::Abs(dval - kDfltDouble) < 1e-4);
}

Bool_t bsim::IsDefault(Bool_t bval) 
{
   return (bval == kDfltBool);
}
