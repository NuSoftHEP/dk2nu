//
// Convert a g4minerva ntuple to the new format
//
// this script can be run using:
//       $ root $(DK2NU)/snippets/load_dk2nu.C convert_g4minerva.C("file.root")
//
// rhatcher@fnal.gov  2012-11-06
//====================================================================

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include "TMath.h"

#include "dk2nu/convert/common_convert.C"

#include "g4minerva.C"

void copy_g4minerva_to_dk2nu(const g4minerva& g4minervaObj);

void convert_g4minerva(string ifname="../generic_g4minerva.root",
                       string inputloc="MINOS",  // "MINOS" or "NOvA"
                       int jobnum=42,
                       Long64_t maxentries=200000000,
                       Long64_t moddump=-1) // modulo for dump
{
  // set globals
  myjob = jobnum; // allow override because g4minerva files forgot to set this
  pots  = 0;
  // allowance in location energy/weight cross-check
  frac_diff_tolerance = 5.0e-4;

  int highest_potnum = 0;

  // open input ntuple
  TFile* fin = TFile::Open(ifname.c_str());
  const char* objName = "nudata";
  TTree* tin = 0;
  fin->GetObject(objName,tin);

  if ( ! tin ) {
    cout << "couldn't find " << objName << endl;
    return;
  }
  cout << "Input file:  " << ifname << endl;

  // use makeclass interface for g4minerva
  g4minerva g4minervaObj(tin);

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
  Long64_t nentries = g4minervaObj.fChain->GetEntriesFast();
  if ( maxentries > 0 ) {
    cout << "limit # of entries to " << nentries 
         << " / " << maxentries << endl;
    nentries = TMath::Min(nentries,maxentries);
  }
  cout << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    Long64_t ientry = g4minervaObj.LoadTree(jentry);
    if (ientry < 0) break;
    nb = g4minervaObj.fChain->GetEntry(jentry);
    nbytes += nb;

    // always clear the dk2nu object 
    dk2nuObj->Clear(); //  !!! important always to this
    // fill the dk2nu object from the g4minerva entry
    copy_g4minerva_to_dk2nu(g4minervaObj);
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
    // don't currently know mapping for minerva weights to dk2nu list
    //compare(g4minervaObj.Nenergyn,dk2nuObj->nuenergy[indxN],"near energy");
    //compare(g4minervaObj.Nwtnear,dk2nuObj->nuwgt[indxN],"near wgt");
    //compare(g4minervaObj.Nenergyf,dk2nuObj->nuenergy[indxF],"far energy");
    //compare(g4minervaObj.Nwtfar,dk2nuObj->nuwgt[indxF],"far wgt");

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

void copy_g4minerva_to_dk2nu(const g4minerva& g4minervaObj)
{
  // fill the global dk2nuObj

  dk2nuObj->norig    = g4minervaObj.Norig;
  dk2nuObj->ndecay   = g4minervaObj.Ndecay;
  dk2nuObj->ntype    = g4minervaObj.Ntype;
  dk2nuObj->vx       = g4minervaObj.Vx;
  dk2nuObj->vy       = g4minervaObj.Vy;
  dk2nuObj->vz       = g4minervaObj.Vz;
  dk2nuObj->pdpx     = g4minervaObj.pdPx;
  dk2nuObj->pdpy     = g4minervaObj.pdPy;
  dk2nuObj->pdpz     = g4minervaObj.pdPz;
  dk2nuObj->ppdxdz   = g4minervaObj.ppdxdz;
  dk2nuObj->ppdydz   = g4minervaObj.ppdydz;
  dk2nuObj->pppz     = g4minervaObj.pppz;
  dk2nuObj->ppenergy = g4minervaObj.ppenergy;
  dk2nuObj->ppmedium = g4minervaObj.ppmedium;
  dk2nuObj->ptype    = g4minervaObj.ptype;
  dk2nuObj->ppvx     = g4minervaObj.ppvx;
  dk2nuObj->ppvy     = g4minervaObj.ppvy;
  dk2nuObj->ppvz     = g4minervaObj.ppvz;
  dk2nuObj->muparpx  = g4minervaObj.muparpx;
  dk2nuObj->muparpy  = g4minervaObj.muparpy;
  dk2nuObj->muparpz  = g4minervaObj.muparpz;
  dk2nuObj->mupare   = g4minervaObj.mupare;

  dk2nuObj->necm     = g4minervaObj.Necm;
  dk2nuObj->nimpwt   = g4minervaObj.Nimpwt;
  dk2nuObj->xpoint   = g4minervaObj.xpoint;
  dk2nuObj->ypoint   = g4minervaObj.ypoint;
  dk2nuObj->zpoint   = g4minervaObj.zpoint;

  dk2nuObj->tvx      = g4minervaObj.tvx;
  dk2nuObj->tvy      = g4minervaObj.tvy;
  dk2nuObj->tvz      = g4minervaObj.tvz;
  dk2nuObj->tpx      = g4minervaObj.tpx;
  dk2nuObj->tpy      = g4minervaObj.tpy;
  dk2nuObj->tpz      = g4minervaObj.tpz;
  dk2nuObj->tptype   = g4minervaObj.tptype;
  dk2nuObj->tgen     = g4minervaObj.tgen;
  //no-equiv//  dk2nuObj->tgptype  = g4minervaObj.tgptype;
  //no-equiv//  dk2nuObj->tgppx    = g4minervaObj.tgppx;
  //no-equiv//  dk2nuObj->tgppy    = g4minervaObj.tgppy;
  //no-equiv//  dk2nuObj->tgppz    = g4minervaObj.tgppz;
  //no-equiv//  dk2nuObj->tprivx   = g4minervaObj.tprivx;
  //no-equiv//  dk2nuObj->tprivy   = g4minervaObj.tprivy;
  //no-equiv//  dk2nuObj->tprivz   = g4minervaObj.tprivz;
  dk2nuObj->beamx    = g4minervaObj.protonX;
  dk2nuObj->beamy    = g4minervaObj.protonY;
  dk2nuObj->beamz    = g4minervaObj.protonZ;
  dk2nuObj->beampx   = g4minervaObj.protonPx;
  dk2nuObj->beampy   = g4minervaObj.protonPy;
  dk2nuObj->beampz   = g4minervaObj.protonPz;

  //cout << "pot # " << g4minervaObj.evtno << endl;

  dk2nuObj->job    = myjob;
  dk2nuObj->potnum = g4minervaObj.evtno;

  // calcLocationWeights needs these filled if it isn't going assert()
  // really need to fill the other bits at this point as well:
  //   ntype, ptype, vx, vy, vz, pdpx, pdpy, pdpz, necm, 
  //   ppenergy, ppdxdz, ppdydz, pppz, 
  //   muparpx, muparpy, muparpz, mupare
  dk2nuObj->ptype    = g4minervaObj.ptype;  
  dk2nuObj->ntype    = g4minervaObj.Ntype;
  dk2nuObj->vx       = g4minervaObj.Vx;
  dk2nuObj->vy       = g4minervaObj.Vy;
  dk2nuObj->vz       = g4minervaObj.Vz;
  dk2nuObj->pdpx     = g4minervaObj.pdPx;
  dk2nuObj->pdpy     = g4minervaObj.pdPy;
  dk2nuObj->pdpz     = g4minervaObj.pdPz;
  dk2nuObj->necm     = g4minervaObj.Necm;
  dk2nuObj->ppenergy = g4minervaObj.ppenergy;
  dk2nuObj->ppdxdz   = g4minervaObj.ppdxdz;
  dk2nuObj->ppdydz   = g4minervaObj.ppdydz;
  dk2nuObj->pppz     = g4minervaObj.pppz;
  dk2nuObj->muparpx  = g4minervaObj.muparpx;
  dk2nuObj->muparpy  = g4minervaObj.muparpy;
  dk2nuObj->muparpz  = g4minervaObj.muparpz;
  dk2nuObj->mupare   = g4minervaObj.mupare;

  dk2nuObj->nimpwt   = g4minervaObj.Nimpwt;

  // calcLocationWeights needs random decay entries filled first
  // but don't copy any other (i.e. near/far) values as they'll
  // be recalculated for the whole set of locations
  double pzrndm = g4minervaObj.Npz;
  double pxrndm = g4minervaObj.Ndxdz * pzrndm;
  double pyrndm = g4minervaObj.Ndydz * pzrndm;
  dk2nuObj->nupx.push_back(pxrndm);
  dk2nuObj->nupy.push_back(pyrndm);
  dk2nuObj->nupz.push_back(pzrndm);
  dk2nuObj->nuenergy.push_back(g4minervaObj.Nenergy);
  dk2nuObj->nuwgt.push_back(1.0);

  dk2nuObj->tptype   = g4minervaObj.tptype;
  //no-equiv//  dk2nuObj->tgptype  = g4minervaObj.tgptype;

  // now copy ancestor history
  int nmx = TMath::Min(g4minervaObj.ntrajectory,10);
  for (int ian=0; ian < nmx; ++ian) {

    const double mm2cm = 0.1;
    const double mev2gev = 0.001;

    dk2nuObj->apdg.push_back(g4minervaObj.pdg[ian]);
    dk2nuObj->trackid.push_back(g4minervaObj.trackId[ian]);

    dk2nuObj->parentid.push_back(g4minervaObj.parentId[ian]);
    dk2nuObj->startx.push_back(g4minervaObj.startx[ian]*mm2cm);
    dk2nuObj->starty.push_back(g4minervaObj.starty[ian]*mm2cm);
    dk2nuObj->startz.push_back(g4minervaObj.startz[ian]*mm2cm);
    dk2nuObj->stopx.push_back(g4minervaObj.stopx[ian]*mm2cm);
    dk2nuObj->stopy.push_back(g4minervaObj.stopy[ian]*mm2cm);
    dk2nuObj->stopz.push_back(g4minervaObj.stopz[ian]*mm2cm);

    dk2nuObj->startpx.push_back(g4minervaObj.startpx[ian]*mev2gev);
    dk2nuObj->startpy.push_back(g4minervaObj.startpy[ian]*mev2gev);
    dk2nuObj->startpz.push_back(g4minervaObj.startpz[ian]*mev2gev);
    dk2nuObj->stoppx.push_back(g4minervaObj.stoppx[ian]*mev2gev);
    dk2nuObj->stoppy.push_back(g4minervaObj.stoppy[ian]*mev2gev);
    dk2nuObj->stoppz.push_back(g4minervaObj.stoppz[ian]*mev2gev);

    dk2nuObj->pprodpx.push_back(g4minervaObj.pprodpx[ian]*mev2gev);
    dk2nuObj->pprodpy.push_back(g4minervaObj.pprodpy[ian]*mev2gev);
    dk2nuObj->pprodpz.push_back(g4minervaObj.pprodpz[ian]*mev2gev);

    // TString -> std::string via char*
    dk2nuObj->proc.push_back(g4minervaObj.proc[ian].Data());
    dk2nuObj->ivol.push_back(g4minervaObj.ivol[ian].Data());
    dk2nuObj->fvol.push_back(g4minervaObj.fvol[ian].Data());

    if ( ian > 0 ) {
      // except for first particle start[i] should == stop[i-1]
      bool match = true;
      match = match && ( close_enough(g4minervaObj.startx[ian],g4minervaObj.stopx[ian-1]) );
      match = match && ( close_enough(g4minervaObj.starty[ian],g4minervaObj.stopy[ian-1]) );
      match = match && ( close_enough(g4minervaObj.startz[ian],g4minervaObj.stopz[ian-1]) );
      if ( ! match ) {
        cout << "## ancestor " << ian << " didn't start where previous entry stopped" << endl;
      }
    }
  }

}

