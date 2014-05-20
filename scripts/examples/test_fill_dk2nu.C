/// \file  test_fill_dk2nu.C
/// \brief Test creating and filling a TTree based on:
///     dk2nu.h  (dk2nu.C)  - decays of beam particles to a neutrino
///     dkmeta.h (dkmeta.C) - metadata for the file
///
/// also show the use of reading location information, generating
/// weights at those locations, and how to extend the tree (for one-off
/// tests) without modifying the standard class. 
///
/// This script can be run using:
///       root -b -q test_fill_dk2nu.C+
///
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \created   2012-04-03
/// \modified  2012-10-03
/// \version $Id: test_fill_dk2nu.C,v 1.3 2012-11-15 21:56:53 rhatcher Exp $
///==========================================================================

#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#ifndef __CINT__
// hide header stuff from CINT, assume load_dk2nu.C run first

#include "tree/dk2nu.h"
#include "tree/dkmeta.h"

/// include standardized code for reading location text file
#include "tree/readWeightLocations.h"

/// include standardized code for getting energy/weight vectors for locations
#include "tree/calcLocationWeights.h"

#endif  // ifndef __CINT__

// if running bare CINT don't try adding the NonStd addition it won't work
#define ADD_NONSTD

#ifdef ADD_NONSTD
/// example class for extending the tree with non-standard extras
/// that doesn't require modifying the basic "dk2nu" class
class NonStd {
  public:
    NonStd() { }
    virtual ~NonStd() { }
    void Clear() { }
    double foo;  ///< data for my one-off test
    double bar;  ///< more data  
  ClassDef(NonStd,1)
};

/// make a dictionary for classes used in the tree
/// again do this because we have no external linkages to libraries
#pragma link C++ class NonStd+;
#endif

// flugg 500K POT lowth files seem to have 510000 as an upper limit on 
// # of entries.   So to test for estimate of file size one needs to have 
// that many entries _and_ semi-sensible values for all branches (so 
// compression isn't better than it would be in real life).
void test_fill_dk2nu(unsigned int nentries=1000)
{

  // stuff...
  TRandom3* rndm = new TRandom3();

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::NumiAnalysis() in g4numi
  ///
  ///-----------------------------------------------------------------------

  // create objects
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
#ifdef ADD_NONSTD
  NonStd* nonstd = new NonStd;
#endif

  // read the text file for locations, fill the dkmeta object
  std::string locfilename = "$DK2NU/etc/locations.txt";
  bsim::readWeightLocations(locfilename,dkmeta);

  // print out what we have for locations
  size_t nloc = dkmeta->location.size();
  std::cout << "Read " << nloc << " locations read from \"" 
            << locfilename << "\"" << std::endl;
  for (size_t iloc = 0; iloc < nloc; ++iloc ) {
    std::cout << "[ " << std::setw(2) << iloc << "] "
              << dkmeta->location[iloc] << std::endl;
  }

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::book() in g4numi
  ///
  ///-----------------------------------------------------------------------

  // create file, book tree, set branch address to created object 
  TFile* treeFile = new TFile("test_dk2nu.root","RECREATE");

  TTree* dk2nuTree = new TTree("dk2nuTree","neutrino ntuple");
  dk2nuTree->Branch("dk2nu","bsim::Dk2Nu",&dk2nu,32000,1);
#ifdef ADD_NONSTD
  // extend the tree with additional branches without modifying std class
  dk2nuTree->Branch("nonstdb","NonStd",&nonstd,32000,1);
#endif

  TTree* dkmetaTree  = new TTree("dkmetaTree","neutrino ntuple metadata");
  dkmetaTree->Branch("dkmeta","bsim::DkMeta",&dkmeta,32000,1);

  int myjob = 42;  // unique identifying job # for this series

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::? in g4numi
  ///  this is the main loop, making entries as the come about
  ///
  ///-----------------------------------------------------------------------
  // fill a few element of a few entries
  for (unsigned int ipot=1; ipot <= nentries; ++ipot) {

    ///
    ///  equivalent to NumiAnalysis::FillNeutrinoNtuple() in g4numi
    ///  (only the part within the loop over ipot)
    ///

    // clear the object in preparation for filling an entry
    dk2nu->clear();

    // fill with info ... only a few elements, just for test purposes
    dk2nu->job    = myjob;
    dk2nu->potnum = ipot;

    // pick a bogus particle type to decay, and a neutrino flavor
    int ptype = 211;  // pi+
    if ( ipot %  5 == 0 ) ptype = 321;  // k+
    if ( ipot % 50 == 0 ) ptype = 13;   // mu-
    int ntype = ( ( ptype == 321 ) ? 12 : 14 );
    TVector3 p3nu(1,2,3); // bogus random neutrino decay vector

    // calcLocationWeights needs these filled if it isn't going assert()
    // really need to fill the other bits at this point as well:
    //   ntype, ptype, vx, vy, vz, pdpx, pdpy, pdpz, necm, 
    //   ppenergy, ppdxdz, ppdydz, pppz, 
    //   muparpx, muparpy, muparpz, mupare
    dk2nu->decay.ptype  = ptype;  
    dk2nu->decay.ntype  = ntype;

    // fill nupx, nupy, nupz, nuenergy, nuwgt(=1) for random decay
    // should be the 0-th entry
    if ( dkmeta->location[0].name == "random decay" ) {
      bsim::NuRay nurndm(p3nu.x(),p3nu.y(),p3nu.z(),p3nu.Mag(),1.0);
      dk2nu->nuray.push_back(nurndm);
    }
    // fill location specific p3, energy and weights; locations in metadata
    bsim::calcLocationWeights(dkmeta,dk2nu);

    // test the filling of vector where entries vary in length
    // ... really need to fill whole dk2nu object
    unsigned int nancestors = rndm->Integer(12) + 1;  // at least one entry
    for (unsigned int janc = 0; janc < nancestors; ++janc ) {
      int xpdg = rndm->Integer(100);
      bsim::Ancestor anc;
      anc.pdg = janc*10000+xpdg;
      dk2nu->ancestor.push_back(anc);
    }

    // push a couple of user defined values for each entry
    dk2nu->vint.push_back(42);
    dk2nu->vint.push_back(ipot);

#ifdef ADD_NONSTD
    // fill non-standard extension to tree with user additions
    nonstd->foo = ptype + 1000000;
    nonstd->bar = ipot + ptype;
#endif

    // push entry out to tree
    dk2nuTree->Fill();

  } // end of fill loop

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::finish() in g4numi
  ///
  ///-----------------------------------------------------------------------

  /// fill the rest of the metadata (locations filled above)
  //no! would clear location info // dkmeta->Clear();
  dkmeta->job  = myjob;  // needs to match the value in each dk2nu entry
  dkmeta->pots = 50000;  // ntuple represents this many protons-on-target 
  dkmeta->beamsim = "test_fill_dk2nu.C";
  dkmeta->physics = "bogus";
  dkmeta->vintnames.push_back("mytemp_42");
  dkmeta->vintnames.push_back("mytemp_ipot");
  // push entry out to meta-data tree
  dkmetaTree->Fill();

  // finish and clean-up
  treeFile->cd();
  dk2nuTree->Write();
  dkmetaTree->Write();
  treeFile->Close();
  delete treeFile; treeFile=0;
  dk2nuTree=0;
  dkmetaTree=0;

}
