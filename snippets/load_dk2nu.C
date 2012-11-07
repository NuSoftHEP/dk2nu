void load_dk2nu(bool genie_too = false, bool verbose = false)
{

  //TString libs0 = gSystem->GetDynamicPath();
  //TString libs  = libs0 + ":/usr/lib:/usr/local/lib:/opt/lib:/opt/local/lib";
  //gSystem->SetDynamicPath(libs.Data());

  /// Add dk2nu libraries
  const char* libname[2] = { "libdk2nuTree", "libdk2nuGenie" };
  size_t nlib = ( genie_too ? 2 : 1 );
  for ( size_t ilib = 0; ilib < nlib; ++ilib ) {
    if ( verbose ) cout << "Load: " << libname[ilib] << endl;
    gSystem->Load(libname[ilib]);
  }

  /// Add locations for include files
  const char* path = gSystem->ExpandPathName("$(DK2NU)");
  if ( ! path ) {
    cout << "DK2NU not defined?!" << endl;
    return;
  }
  TString ip = gSystem->GetIncludePath();
  ip += " -I";
  ip += path;
  ip += " -I";
  ip += path;
  ip += "/include";
  gSystem->SetIncludePath(ip);
  if ( verbose ) cout << "SetIncludePath:  " << ip << endl;

}
