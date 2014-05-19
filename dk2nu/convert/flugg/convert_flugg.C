//
// Convert a flugg ntuple to the new format
//
// this script can be run using:
//       $ root $(DK2NU)/snippets/load_dk2nu.C convert_flugg.C("file.root")
//
// rhatcher@fnal.gov  2012-11-06
//====================================================================

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include "TMath.h"

#include "convert/common_convert.C"
#include "convert/flugg/flugg.C"

void copy_flugg_to_dk2nu(const flugg& fluggObj);
void fluggCrossChecks(const flugg& fluggObj, string inputloc);

void convert_flugg(string ifname="../fluxfiles/generic_flugg.root",
                   int jobnum=42,
                   string inputloc="MINOS",  // "MINOS"or "NOvA" for xcheck
                   Long64_t maxentries=-1,
                   Long64_t moddump=-1) // modulo for dump
{
  // set globals
  myjob = jobnum; // allow override because flugg files forgot to set this
  pots  = 0;
  // allowance in location energy/weight cross-check
  frac_diff_tolerance = 2.5e-4;

  int highest_potnum = 0;

  // open input ntuple
  TFile* fin = TFile::Open(ifname.c_str());
  if ( ! fin ) {
    cout << "couldn't open input file: " << ifname << endl;
    return;
  }
  const char* objName = "h10";
  TTree* tin = 0;
  fin->GetObject(objName,tin);

  if ( ! tin ) {
    cout << "couldn't find " << objName << endl;
    return;
  }
  cout << endl << "Input file:  " << ifname << endl;

  // use makeclass interface for flugg
  flugg fluggObj(tin);

  // construct output filename
  string ofname = construct_outfilename(ifname);
  cout << "Output file: " << ofname << endl << endl;

  ConvertInit();
  ConvertBookNtuple(ofname);

  ///
  /// loop over entries
  ///
  Long64_t nentries = fluggObj.fChain->GetEntriesFast();
  if ( maxentries > 0 ) {
    cout << "limit # of entries to " << nentries 
         << " / " << maxentries << endl;
    nentries = TMath::Min(nentries,maxentries);
  }
  cout << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    Long64_t ientry = fluggObj.LoadTree(jentry);
    if (ientry < 0) break;
    nb = fluggObj.fChain->GetEntry(jentry);
    nbytes += nb;

    // always clear the dk2nu object 
    dk2nu->clear(); //  !!! important !!! always do this

    // fill the dk2nu object from the flugg entry
    copy_flugg_to_dk2nu(fluggObj);

    // fill location specific p3, energy and weights
    // locations to fill are in the metadata
    // assumes that prior copying filled the first entry w/ random decay
    calcLocationWeights(dkmeta,dk2nu);

    // keep track of potnum
    if ( dk2nu->potnum > highest_potnum ) 
      highest_potnum = dk2nu->potnum;

    // push entry out to tree
    dk2nuTree->Fill();

    // just for fun print every n entries
    if ( moddump > 0 && jentry%moddump == 0 ) cout << endl << *dk2nu << endl;

    fluggCrossChecks(fluggObj,inputloc);

  }
  cout << endl;

  pots = estimate_pots(highest_potnum);
  cout << "estimated pots " << pots << " (" << highest_potnum << ")" << endl;
  cout << "maximum observed fractional difference in cross checks: "
       << obs_frac_diff_max << endl
       << "print ## message when exceeded " << frac_diff_tolerance << endl;

  dkmeta->job     = myjob;
  dkmeta->pots    = pots;
  dkmeta->beamsim = "convert_flugg.C";
  dkmeta->physics = "bogus";

  // print last entries
  cout << endl << "========== last entries" << endl << endl;
  cout << *dk2nu << endl << endl;
  cout << *dkmeta << endl;

  // finish off ntuples (including meta data)
  ConvertFinish();

  // close input file
  fin->Close();
}

void copy_flugg_to_dk2nu(const flugg& fluggObj)
{
  // fill the global dk2nu object

  dk2nu->job    = myjob;
  dk2nu->potnum = fluggObj.evtno;

  // calcLocationWeights needs random decay entries filled first
  // but don't copy any other (i.e. near/far) values as they'll
  // be recalculated for the whole set of locations
  double pzrndm = fluggObj.Npz;
  double pxrndm = fluggObj.Ndxdz * pzrndm;
  double pyrndm = fluggObj.Ndydz * pzrndm;
  bsim::NuRay nuray(pxrndm,pyrndm,pzrndm,fluggObj.Nenergy,1.0);
  dk2nu->nuray.push_back(nuray);

  // calcLocationWeights needs these filled if it isn't going assert()
  // really need to fill the other bits at this point as well:
  //   ntype, ptype, vx, vy, vz, pdpx, pdpy, pdpz, necm, 
  //   ppenergy, ppdxdz, ppdydz, pppz, 
  //   muparpx, muparpy, muparpz, mupare

  dk2nu->decay.norig    = fluggObj.Norig;
  dk2nu->decay.ndecay   = fluggObj.Ndecay;
  dk2nu->decay.ntype    = Convert5xToPdg(fluggObj.Ntype);
  dk2nu->decay.vx       = fluggObj.Vx;
  dk2nu->decay.vy       = fluggObj.Vy;
  dk2nu->decay.vz       = fluggObj.Vz;
  dk2nu->decay.pdpx     = fluggObj.pdPx;
  dk2nu->decay.pdpy     = fluggObj.pdPy;
  dk2nu->decay.pdpz     = fluggObj.pdPz;
  dk2nu->decay.ppdxdz   = fluggObj.ppdxdz;
  dk2nu->decay.ppdydz   = fluggObj.ppdydz;
  dk2nu->decay.pppz     = fluggObj.pppz;
  dk2nu->decay.ppenergy = fluggObj.ppenergy;
  dk2nu->decay.ppmedium = fluggObj.ppmedium;
  dk2nu->decay.ptype    = ConvertGeantToPdg(fluggObj.ptype,"ptype");
  dk2nu->decay.muparpx  = fluggObj.muparpx;
  dk2nu->decay.muparpy  = fluggObj.muparpy;
  dk2nu->decay.muparpz  = fluggObj.muparpz;
  dk2nu->decay.mupare   = fluggObj.mupare;

  dk2nu->decay.necm     = fluggObj.Necm;
  dk2nu->decay.nimpwt   = fluggObj.Nimpwt;

  dk2nu->ppvx     = fluggObj.ppvx;
  dk2nu->ppvy     = fluggObj.ppvy;
  dk2nu->ppvz     = fluggObj.ppvz;

  //not-in-new//dk2nu->xpoint   = fluggObj.xpoint;
  //not-in-new//dk2nu->ypoint   = fluggObj.ypoint;
  //not-in-new//dk2nu->zpoint   = fluggObj.zpoint;

  dk2nu->tgtexit.tvx      = fluggObj.tvx;
  dk2nu->tgtexit.tvy      = fluggObj.tvy;
  dk2nu->tgtexit.tvz      = fluggObj.tvz;
  dk2nu->tgtexit.tpx      = fluggObj.tpx;
  dk2nu->tgtexit.tpy      = fluggObj.tpy;
  dk2nu->tgtexit.tpz      = fluggObj.tpz;
  dk2nu->tgtexit.tptype   = ConvertGeantToPdg(fluggObj.tptype,"tptype");
  dk2nu->tgtexit.tgen     = fluggObj.tgen;

  //no-equiv//  dk2nu->tgptype  = ConvertGeantToPdg(fluggObj.tgptype,"tgptype");
  //no-equiv//  dk2nu->tgppx    = fluggObj.tgppx;
  //no-equiv//  dk2nu->tgppy    = fluggObj.tgppy;
  //no-equiv//  dk2nu->tgppz    = fluggObj.tgppz;
  //no-equiv//  dk2nu->tprivx   = fluggObj.tprivx;
  //no-equiv//  dk2nu->tprivy   = fluggObj.tprivy;
  //no-equiv//  dk2nu->tprivz   = fluggObj.tprivz;

  //not-in-new//dk2nu->beamx    = fluggObj.protonX;
  //not-in-new//dk2nu->beamy    = fluggObj.protonY;
  //not-in-new//dk2nu->beamz    = fluggObj.protonZ;
  //not-in-new//dk2nu->beampx   = fluggObj.protonPx;
  //not-in-new//dk2nu->beampy   = fluggObj.protonPy;
  //not-in-new//dk2nu->beampz   = fluggObj.protonPz;

  //no-equiv//  dk2nu->tgptype  = fluggObj.tgptype;

}

void fluggCrossChecks(const flugg& fluggObj, string inputloc)
{

  static bool first = true;
  static size_t indxNear = -1, indxFar = -1;
  if ( first ) {
    first = false;
    /// find indices for new weights so we can compare later
    indxNear = find_loc_index(inputloc+" NearDet");
    indxFar  = find_loc_index(inputloc+" FarDet");
    cout << "indx " << inputloc << " Near " << indxNear 
         << " Far " << indxFar << endl;
  }
  // cross check location energy/weights
  // 
  histCompare(fluggObj.Nenergyn,dk2nu->nuray[indxNear].E,  "near energy");
  histCompare(fluggObj.Nwtnear, dk2nu->nuray[indxNear].wgt,"near wgt");
  histCompare(fluggObj.Nenergyf,dk2nu->nuray[indxFar].E,   "far energy");
  histCompare(fluggObj.Nwtfar,  dk2nu->nuray[indxFar].wgt, "far wgt");

}

