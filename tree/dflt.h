#include "TObject.h"   // to get root types

/// bsim namespace for beam simulation classes and functions
namespace bsim
{
    Bool_t IsDefault(Float_t  intval);
    Bool_t IsDefault(Int_t    intval);
    Bool_t IsDefault(UInt_t   intval);
    Bool_t IsDefault(Double_t intval);
    Bool_t IsDefault(Bool_t   intval);

    static const Float_t  kDfltFloat  = -9999.99;
    static const Int_t    kDfltInt    = -9999;
    static const UInt_t   kDfltUInt   =  9999;
    static const Double_t kDfltDouble = -9999.99;
    static const Bool_t   kDfltBool   = false;
}
