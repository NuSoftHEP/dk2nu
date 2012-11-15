#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include "TChain.h"

#ifndef __CINT__
// hide header stuff from CINT, assume load_dk2nu.C run first

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#endif  // ifndef __CINT__

void test_read_dk2nu(string pattern="generic_*_to_dk2nu.root")
{
  TChain* cflux = new TChain("dk2nuTree");
  TChain* cmeta = new TChain("dkmetaTree");

  cflux->Add(pattern.c_str());
  cmeta->Add(pattern.c_str());
  
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
  cflux->SetBranchAddress("dk2nu",&dk2nu);
  cmeta->SetBranchAddress("dkmeta",&dkmeta);

  cout << "before reading any entries" << endl;
  cout << *dk2nu << endl << endl;
  cout << *dkmeta << endl << endl;

  Long64_t nflux = cflux->GetEntries();
  Long64_t nmeta = cmeta->GetEntries();
  cout << "nentries:  " << nflux << " " << nmeta << endl;
  
  for (Long64_t i=0; i < nflux; ++i ) {
    cflux->GetEntry(i);
    if ( i < 5 ) cout << "ntype " << dk2nu->decay.ntype << endl;
  }
  cout << endl << *dk2nu << endl << endl;

  cmeta->GetEntry(0);
  cout << *dkmeta << endl;

}
