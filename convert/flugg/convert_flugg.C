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

#include "dk2nu/convert/common_convert.C"

#include "flugg.C"

void copy_flugg_to_dk2nu(const flugg& fluggObj);

void convert_flugg(string ifname="../generic_flugg.root",
                   string inputloc="MINOS",  // "MINOS" or "NOvA"
                   int jobnum=42,
                   Long64_t maxentries=200000000,
                   Long64_t moddump=-1) // modulo for dump
{
  // set globals
  myjob = jobnum; // allow override because flugg files forgot to set this
  pots  = 0;
  // allowance in location energy/weight cross-check
  frac_diff_tolerance = 5.0e-4;

  int highest_potnum = 0;

  // open input ntuple
  TFile* fin = TFile::Open(ifname.c_str());
  const char* objName = "h10";
  TTree* tin = 0;
  fin->GetObject(objName,tin);

  if ( ! tin ) {
    cout << "couldn't find " << objName << endl;
    return;
  }
  cout << "Input file:  " << ifname << endl;

  // use makeclass interface for flugg
  flugg fluggObj(tin);

  // construct output filename
  string ofname = construct_outfilename(ifname);
  cout << "Output file: " << ofname << endl;

  ConvertInit();
  ConvertBookNtuple(ofname);

  /// find indices for new weights so we can compare later
  size_t indxN = find_loc_index(inputloc+" NearDet");
  size_t indxF = find_loc_index(inputloc+" FarDet");
  cout << "indx Near " << indxN << " Far " << indxF << endl;

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
    dk2nuObj->Clear(); //  !!! important always to this
    // fill the dk2nu object from the flugg entry
    copy_flugg_to_dk2nu(fluggObj);
    // convert Geant to PDG codes (ntype,ptype,tptype,tgptype) in dk2nu
    ConvertPartCodes();
    // fill location specific p3, energy and weights
    // locations to fill are in the metadata
    // assumes that prior copying filled the first entry w/ random decay
    calcLocationWeights(dkmetaObj,dk2nuObj);
    // keep track of potnum
    if ( dk2nuObj->potnum > highest_potnum ) 
      highest_potnum = dk2nuObj->potnum;
    // push entry out to tree
    dk2nu_tree->Fill();

    // just for fun print every n entries
    if ( moddump > 0 && jentry%moddump == 0 ) cout << endl << *dk2nuObj << endl;

    // cross check location energy/weights
    compare(fluggObj.Nenergyn,dk2nuObj->nuenergy[indxN],"near energy");
    compare(fluggObj.Nwtnear,dk2nuObj->nuwgt[indxN],"near wgt");
    compare(fluggObj.Nenergyf,dk2nuObj->nuenergy[indxF],"far energy");
    compare(fluggObj.Nwtfar,dk2nuObj->nuwgt[indxF],"far wgt");

  }
  cout << endl;

  pots = estimate_pots(highest_potnum);
  cout << "estimated pots " << pots << " (" << highest_potnum << ")" << endl;
  cout << "maximum observed fractional difference in cross checks: "
       << obs_frac_diff_max << endl
       << "print ## message when exceeded " << frac_diff_tolerance << endl;

  // finish off ntuples (including meta data)
  ConvertFinish();

  // close input file
  fin->Close();
}

void copy_flugg_to_dk2nu(const flugg& fluggObj)
{
  // fill the global dk2nuObj

  dk2nuObj->norig    = fluggObj.Norig;
  dk2nuObj->ndecay   = fluggObj.Ndecay;
  dk2nuObj->ntype    = fluggObj.Ntype;
  dk2nuObj->vx       = fluggObj.Vx;
  dk2nuObj->vy       = fluggObj.Vy;
  dk2nuObj->vz       = fluggObj.Vz;
  dk2nuObj->pdpx     = fluggObj.pdPx;
  dk2nuObj->pdpy     = fluggObj.pdPy;
  dk2nuObj->pdpz     = fluggObj.pdPz;
  dk2nuObj->ppdxdz   = fluggObj.ppdxdz;
  dk2nuObj->ppdydz   = fluggObj.ppdydz;
  dk2nuObj->pppz     = fluggObj.pppz;
  dk2nuObj->ppenergy = fluggObj.ppenergy;
  dk2nuObj->ppmedium = fluggObj.ppmedium;
  dk2nuObj->ptype    = fluggObj.ptype;
  dk2nuObj->ppvx     = fluggObj.ppvx;
  dk2nuObj->ppvy     = fluggObj.ppvy;
  dk2nuObj->ppvz     = fluggObj.ppvz;
  dk2nuObj->muparpx  = fluggObj.muparpx;
  dk2nuObj->muparpy  = fluggObj.muparpy;
  dk2nuObj->muparpz  = fluggObj.muparpz;
  dk2nuObj->mupare   = fluggObj.mupare;

  dk2nuObj->necm     = fluggObj.Necm;
  dk2nuObj->nimpwt   = fluggObj.Nimpwt;
  dk2nuObj->xpoint   = fluggObj.xpoint;
  dk2nuObj->ypoint   = fluggObj.ypoint;
  dk2nuObj->zpoint   = fluggObj.zpoint;

  dk2nuObj->tvx      = fluggObj.tvx;
  dk2nuObj->tvy      = fluggObj.tvy;
  dk2nuObj->tvz      = fluggObj.tvz;
  dk2nuObj->tpx      = fluggObj.tpx;
  dk2nuObj->tpy      = fluggObj.tpy;
  dk2nuObj->tpz      = fluggObj.tpz;
  dk2nuObj->tptype   = fluggObj.tptype;
  dk2nuObj->tgen     = fluggObj.tgen;
  dk2nuObj->tgptype  = fluggObj.tgptype;
  dk2nuObj->tgppx    = fluggObj.tgppx;
  dk2nuObj->tgppy    = fluggObj.tgppy;
  dk2nuObj->tgppz    = fluggObj.tgppz;
  dk2nuObj->tprivx   = fluggObj.tprivx;
  dk2nuObj->tprivy   = fluggObj.tprivy;
  dk2nuObj->tprivz   = fluggObj.tprivz;
  dk2nuObj->beamx    = fluggObj.beamx;
  dk2nuObj->beamy    = fluggObj.beamy;
  dk2nuObj->beamz    = fluggObj.beamz;
  dk2nuObj->beampx   = fluggObj.beampx;
  dk2nuObj->beampy   = fluggObj.beampy;
  dk2nuObj->beampz   = fluggObj.beampz;

  //cout << "pot # " << fluggObj.evtno << endl;

  dk2nuObj->job    = myjob;
  dk2nuObj->potnum = fluggObj.evtno;

  // calcLocationWeights needs these filled if it isn't going assert()
  // really need to fill the other bits at this point as well:
  //   ntype, ptype, vx, vy, vz, pdpx, pdpy, pdpz, necm, 
  //   ppenergy, ppdxdz, ppdydz, pppz, 
  //   muparpx, muparpy, muparpz, mupare
  dk2nuObj->ptype    = fluggObj.ptype;  
  dk2nuObj->ntype    = fluggObj.Ntype;
  dk2nuObj->vx       = fluggObj.Vx;
  dk2nuObj->vy       = fluggObj.Vy;
  dk2nuObj->vz       = fluggObj.Vz;
  dk2nuObj->pdpx     = fluggObj.pdPx;
  dk2nuObj->pdpy     = fluggObj.pdPy;
  dk2nuObj->pdpz     = fluggObj.pdPz;
  dk2nuObj->necm     = fluggObj.Necm;
  dk2nuObj->ppenergy = fluggObj.ppenergy;
  dk2nuObj->ppdxdz   = fluggObj.ppdxdz;
  dk2nuObj->ppdydz   = fluggObj.ppdydz;
  dk2nuObj->pppz     = fluggObj.pppz;
  dk2nuObj->muparpx  = fluggObj.muparpx;
  dk2nuObj->muparpy  = fluggObj.muparpy;
  dk2nuObj->muparpz  = fluggObj.muparpz;
  dk2nuObj->mupare   = fluggObj.mupare;

  dk2nuObj->nimpwt   = fluggObj.Nimpwt;

  // calcLocationWeights needs random decay entries filled first
  // but don't copy any other (i.e. near/far) values as they'll
  // be recalculated for the whole set of locations
  double pzrndm = fluggObj.Npz;
  double pxrndm = fluggObj.Ndxdz * pzrndm;
  double pyrndm = fluggObj.Ndydz * pzrndm;
  dk2nuObj->nupx.push_back(pxrndm);
  dk2nuObj->nupy.push_back(pyrndm);
  dk2nuObj->nupz.push_back(pzrndm);
  dk2nuObj->nuenergy.push_back(fluggObj.Nenergy);
  dk2nuObj->nuwgt.push_back(1.0);

  dk2nuObj->tptype   = fluggObj.tptype;
  dk2nuObj->tgptype  = fluggObj.tgptype;

}

