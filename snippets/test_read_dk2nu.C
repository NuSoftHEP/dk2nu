#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include "TChain.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

void test_read_dk2nu(string pattern="generic_g4minerva*.root")
{
  TChain* cflux = new TChain("dk2nu");
  TChain* cmeta = new TChain("dkmeta");

  cflux->AddFile(pattern.c_str());
  cmeta->AddFile(pattern.c_str());
  
  dk2nu*  dk2nuObj  = new dk2nu;
  dkmeta* dkmetaObj = new dkmeta;
  cflux->SetBranchAddress("dk2nu",&dk2nuObj);
  cmeta->SetBranchAddress("dkmeta",&dkmetaObj);

  Long64_t nflux = cflux->GetEntries();
  Long64_t nmeta = cmeta->GetEntries();
  cout << "nentries:  " << nflux << " " << nmeta << endl;
  
  for (Long64_t i=0; i < nflux; ++i ) {
    cflux->GetEntry(i);
    if ( i < 50 ) cout << "ntype " << dk2nuObj->ntype << endl;
  }

  cmeta->GetEntry(0);
  cout << *dkmetaObj << endl;
}
