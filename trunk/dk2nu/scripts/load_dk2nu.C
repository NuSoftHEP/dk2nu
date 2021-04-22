int load_dk2nu(bool genie_too = false, bool cnvtinc = false,
                bool verbose = false)
{
  // load the dk2nuTree library
  // if "genie_too" load the library that had the GENIE flux driver
  // if "cvntinc" add include paths to convert areas

  const char* path = gSystem->ExpandPathName("$(DK2NU)");
  if ( ! path ) {
    cout << "DK2NU not defined?!" << endl;
  }

  // add to the library path
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());
  }

  /// Add dk2nu libraries
  const char* libname[2] = { "libdk2nuTree", "libdk2nuGenie" };
  size_t nlib = ( genie_too ? 2 : 1 );
  for ( size_t ilib = 0; ilib < nlib; ++ilib ) {
    if ( verbose ) cout << "Load: " << libname[ilib] << endl;
    gSystem->Load(libname[ilib]);
  }

  // none of the rest makes sense unless we know the DK2NU path
  if ( ! path ) return -1;

  /// Add locations for include files
  TString ip = gSystem->GetIncludePath();
  ip += " -I";
  ip += path;

  ip += " -I";
  ip += path;
  ip += "/include";

  ip += " -I";
  ip += path;
  ip += "/include/dk2nu";

  ip += " -I";
  ip += path;
  ip += "/scripts";

  const char* cnvt[5] = { "g3numi", "flugg", "g4numi", "g4minerva", "g4lbne" };
  size_t ncnvt = sizeof(cnvt)/sizeof(const char*);
  if ( cnvtinc ) {
    for (size_t i = 0; i < ncnvt; ++i ) {
      ip += " -I";
      ip += path;
      ip += "/scripts/convert/";
      ip += cnvt[i];
    }
  }

  gSystem->SetIncludePath(ip);
  if ( verbose ) cout << "SetIncludePath:  " << ip << endl;

  // additions to .include must be done individually or CINT will
  // try to quote all the spaces as a single path
  TString dip;

  dip = ".include ";
  dip += path;
  gROOT->ProcessLine(dip.Data());

  dip = ".include ";
  dip += path;
  dip += "/include";
  gROOT->ProcessLine(dip.Data());

  dip = ".include ";
  dip += path;
  dip += "/include/dk2nu";
  gROOT->ProcessLine(dip.Data());

  dip = ".include ";
  dip += path;
  dip += "/scripts";
  gROOT->ProcessLine(dip.Data());

  if ( cnvtinc ) {
    for (size_t i = 0; i < ncnvt; ++i ) {
      dip = ".include ";
      dip += path;
      dip += "/scripts/convert/";
      dip += cnvt[i];
      gROOT->ProcessLine(dip.Data());
    }
  }

  return 0;
}
