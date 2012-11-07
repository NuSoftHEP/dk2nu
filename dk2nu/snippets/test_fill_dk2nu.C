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
/// \version $Id: test_fill_dk2nu.C,v 1.1 2012-11-07 05:21:56 rhatcher Exp $
///==========================================================================

#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "dk2nu.h"
#include "dkmeta.h"

/// include these because we're not linking to anything external
/// so we need to include the source for dk2nu::Clear() and dkmeta::Clear()
#include "dk2nu.cc"
#include "dkmeta.cc"

/// example class for extending the tree with non-standard extras
/// that doesn't require modifying the basic "dk2nu" class
class nonstd {
  public:
    nonstd() { }
    virtual ~nonstd() { }
    void Clear() { }
    double foo;  ///< data for my one-off test
    double bar;  ///< more data  
  ClassDef(nonstd,1)
};

/// make a dictionary for classes used in the tree
/// again do this because we have no external linkages to libraries
#ifdef __CINT__
#pragma link C++ class dk2nu+;
#pragma link C++ class dkmeta+;
#pragma link C++ class nonstd+;
#endif

/// include standardized code for reading location text file
#include "readWeightLocations.C"

/// include standardized code for getting energy/weight vectors for locations
#include "calcLocationWeights.C"

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
  dk2nu*  dk2nuObj  = new dk2nu;
  dkmeta* dkmetaObj = new dkmeta;
  nonstd* nonstdObj = new nonstd;

  // read the text file for locations, fill the dkmeta object
  std::string locfilename = "locfile.txt";
  readWeightLocations(locfilename,dkmetaObj);

  // print out what we have for locations
  size_t nloc = dkmetaObj->nameloc.size();
  std::cout << "Read " << nloc << " locations read from \"" 
            << locfilename << "\"" << std::endl;
  for (size_t iloc = 0; iloc < nloc; ++iloc ) {
    std::cout << "{" << std::setw(10) << dkmetaObj->xloc[iloc]
              << "," << std::setw(10) << dkmetaObj->yloc[iloc]
              << "," << std::setw(10) << dkmetaObj->zloc[iloc]
              << " } \"" << dkmetaObj->nameloc[iloc] << "\""
              << std::endl;
  }

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::book() in g4numi
  ///
  ///-----------------------------------------------------------------------

  // create file, book tree, set branch address to created object 
  TFile* treeFile = new TFile("test_dk2nu.root","RECREATE");

  TTree* dk2nu_tree = new TTree("dk2nu","FNAL neutrino ntuple");
  dk2nu_tree->Branch("dk2nu","dk2nu",&dk2nuObj,32000,1);
  // extend the tree with additional branches without modifying std class
  dk2nu_tree->Branch("nonstd","nonstd",&nonstdObj,32000,1);

  TTree* dkmeta_tree  = new TTree("dkmeta","FNAL neutrino ntuple metadata");
  dkmeta_tree->Branch("dkmeta","dkmeta",&dkmetaObj,32000,1);

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
    dk2nuObj->Clear();

    // fill with info ... only a few elements, just for test purposes
    dk2nuObj->job    = myjob;
    dk2nuObj->potnum = ipot;

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
    dk2nuObj->ptype  = ptype;  
    dk2nuObj->ntype  = ntype;

    // fill nupx, nupy, nupz, nuenergy, nuwgt(=1) for random decay
    // should be the 0-th entry
    if ( dkmetaObj->nameloc[0] == "random decay" ) {
      dk2nuObj->nupx.push_back(p3nu.x());
      dk2nuObj->nupy.push_back(p3nu.y());
      dk2nuObj->nupz.push_back(p3nu.z());
      dk2nuObj->nuenergy.push_back(p3nu.Mag());
      dk2nuObj->nuwgt.push_back(1.0);
    }
    // fill location specific p3, energy and weights; locations in metadata
    calcLocationWeights(dkmetaObj,dk2nuObj);

    // test the filling of vector where entries vary in length
    // ... really need to fill whole dk2nu object
    unsigned int nancestors = rndm->Integer(12) + 1;  // at least one entry
    for (unsigned int janc = 0; janc < nancestors; ++janc ) {
      int xpdg = rndm->Integer(100);
      dk2nuObj->apdg.push_back(janc*10000+xpdg);
    }

    // push a couple of user defined values for each entry
    dk2nuObj->vint.push_back(42);
    dk2nuObj->vint.push_back(ipot);

    // fill non-standard extension to tree with user additions
    nonstdObj->foo = ptype + 1000000;
    nonstdObj->bar = ipot + ptype;

    // push entry out to tree
    dk2nu_tree->Fill();

  } // end of fill loop

  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::finish() in g4numi
  ///
  ///-----------------------------------------------------------------------

  /// fill the rest of the metadata (locations filled above)
  //no! would clear location info // dkmetaObj->Clear();
  dkmetaObj->job  = myjob;  // needs to match the value in each dk2nu entry
  dkmetaObj->pots = 50000;  // ntuple represents this many protons-on-target 
  dkmetaObj->beamsim = "test_fill_dk2nu.C";
  dkmetaObj->physics = "bogus";
  dkmetaObj->vintnames.push_back("mytemp_42");
  dkmetaObj->vintnames.push_back("mytemp_ipot");
  // push entry out to meta-data tree
  dkmeta_tree->Fill();

  // finish and clean-up
  treeFile->cd();
  dk2nu_tree->Write();
  dkmeta_tree->Write();
  treeFile->Close();
  delete treeFile; treeFile=0;
  dk2nu_tree=0;
  dkmeta_tree=0;
}
