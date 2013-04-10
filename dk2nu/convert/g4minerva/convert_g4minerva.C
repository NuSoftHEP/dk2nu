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

#include "dk2nu/convert/g4minerva/g4minerva.C"

void copy_g4minerva_to_dk2nu(const g4minerva& g4minervaObj);
void g4minervaCrossChecks(const g4minerva& g4minervaObj);

void convert_g4minerva(string ifname="../fluxfiles/generic_g4minerva.root",
                       int jobnum=42,
                       Long64_t maxentries=-1,
                       Long64_t moddump=-1) // modulo for dump
{
  // set globals
  myjob = jobnum; // allow override because g4minerva files forgot to set this
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
  const char* objName = "nudata";
  TTree* tin = 0;
  fin->GetObject(objName,tin);

  if ( ! tin ) {
    cout << "couldn't find " << objName << endl;
    return;
  }
  cout << endl << "Input file:  " << ifname << endl;

  // use makeclass interface for g4minerva
  g4minerva g4minervaObj(tin);

  // construct output filename
  string ofname = construct_outfilename(ifname);
  cout << "Output file: " << ofname << endl << endl;

  ConvertInit();
  ConvertBookNtuple(ofname);

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
    dk2nu->clear(); //  !!! important !!! always do this

    // fill the dk2nu object from the g4minerva entry
    copy_g4minerva_to_dk2nu(g4minervaObj);

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

    g4minervaCrossChecks(g4minervaObj);

  }
  cout << endl;

  pots = estimate_pots(highest_potnum);
  cout << "estimated pots " << pots << " (" << highest_potnum << ")" << endl;
  cout << "maximum observed fractional difference in cross checks: "
       << obs_frac_diff_max << endl
       << "print ## message when exceeded " << frac_diff_tolerance << endl;

  dkmeta->job     = myjob;
  dkmeta->pots    = pots;
  dkmeta->beamsim = "convert_g4minerva.C";
  dkmeta->physics = "bogus";

  dkmeta->tgtcfg  = TString::Format("Z=%6.1fcm",g4minervaObj.nuTarZ);
  dkmeta->horncfg = TString::Format("I=%6.1fkA",g4minervaObj.hornCurrent);

  dkmeta->beam0x = g4minervaObj.beamX;
  dkmeta->beam0y = g4minervaObj.beamY;
  dkmeta->beamhwidth = g4minervaObj.beamHWidth;
  dkmeta->beamvwidth = g4minervaObj.beamVWidth;

  // print last entries
  cout << endl << "========== last entries" << endl << endl;
  cout << *dk2nu << endl << endl;
  cout << *dkmeta << endl;

  // finish off ntuples (including meta data)
  ConvertFinish();

  // close input file
  fin->Close();
}

void copy_g4minerva_to_dk2nu(const g4minerva& g4minervaObj)
{
  // fill the global dk2nu object

  dk2nu->job    = myjob;
  dk2nu->potnum = g4minervaObj.evtno;

  // calcLocationWeights needs random decay entries filled first
  // but don't copy any other (i.e. near/far) values as they'll
  // be recalculated for the whole set of locations
  double pzrndm = g4minervaObj.Npz;
  double pxrndm = g4minervaObj.Ndxdz * pzrndm;
  double pyrndm = g4minervaObj.Ndydz * pzrndm;
  bsim::NuRay nuray(pxrndm,pyrndm,pzrndm,g4minervaObj.Nenergy,1.0);
  dk2nu->nuray.push_back(nuray);

  // calcLocationWeights needs these filled if it isn't going assert()
  // really need to fill the other bits at this point as well:
  //   ntype, ptype, vx, vy, vz, pdpx, pdpy, pdpz, necm, 
  //   ppenergy, ppdxdz, ppdydz, pppz, 
  //   muparpx, muparpy, muparpz, mupare

  dk2nu->decay.norig    = g4minervaObj.Norig;
  dk2nu->decay.ndecay   = g4minervaObj.Ndecay;
  dk2nu->decay.ntype    = Convert5xToPdg(g4minervaObj.Ntype);
  dk2nu->decay.vx       = g4minervaObj.Vx;
  dk2nu->decay.vy       = g4minervaObj.Vy;
  dk2nu->decay.vz       = g4minervaObj.Vz;
  dk2nu->decay.pdpx     = g4minervaObj.pdPx;
  dk2nu->decay.pdpy     = g4minervaObj.pdPy;
  dk2nu->decay.pdpz     = g4minervaObj.pdPz;
  dk2nu->decay.ppdxdz   = g4minervaObj.ppdxdz;
  dk2nu->decay.ppdydz   = g4minervaObj.ppdydz;
  dk2nu->decay.pppz     = g4minervaObj.pppz;
  dk2nu->decay.ppenergy = g4minervaObj.ppenergy;
  dk2nu->decay.ppmedium = g4minervaObj.ppmedium;
  dk2nu->decay.ptype    = ConvertGeantToPdg(g4minervaObj.ptype,"ptype");
  dk2nu->decay.muparpx  = g4minervaObj.muparpx;
  dk2nu->decay.muparpy  = g4minervaObj.muparpy;
  dk2nu->decay.muparpz  = g4minervaObj.muparpz;
  dk2nu->decay.mupare   = g4minervaObj.mupare;

  dk2nu->decay.necm     = g4minervaObj.Necm;
  dk2nu->decay.nimpwt   = g4minervaObj.Nimpwt;

  dk2nu->ppvx     = g4minervaObj.ppvx;
  dk2nu->ppvy     = g4minervaObj.ppvy;
  dk2nu->ppvz     = g4minervaObj.ppvz;

  //not-in-new//dk2nu->xpoint   = g4minervaObj.xpoint;
  //not-in-new//dk2nu->ypoint   = g4minervaObj.ypoint;
  //not-in-new//dk2nu->zpoint   = g4minervaObj.zpoint;

  dk2nu->tgtexit.tvx      = g4minervaObj.tvx;
  dk2nu->tgtexit.tvy      = g4minervaObj.tvy;
  dk2nu->tgtexit.tvz      = g4minervaObj.tvz;
  dk2nu->tgtexit.tpx      = g4minervaObj.tpx;
  dk2nu->tgtexit.tpy      = g4minervaObj.tpy;
  dk2nu->tgtexit.tpz      = g4minervaObj.tpz;
  dk2nu->tgtexit.tptype   = ConvertGeantToPdg(g4minervaObj.tptype,"tptype");
  dk2nu->tgtexit.tgen     = g4minervaObj.tgen;

  //no-equiv//  dk2nu->tgptype  = ConvertGeantToPdg(g4minervaObj.tgptype,"tgptype");
  //no-equiv//  dk2nu->tgppx    = g4minervaObj.tgppx;
  //no-equiv//  dk2nu->tgppy    = g4minervaObj.tgppy;
  //no-equiv//  dk2nu->tgppz    = g4minervaObj.tgppz;
  //no-equiv//  dk2nu->tprivx   = g4minervaObj.tprivx;
  //no-equiv//  dk2nu->tprivy   = g4minervaObj.tprivy;
  //no-equiv//  dk2nu->tprivz   = g4minervaObj.tprivz;

  //not-in-new//dk2nu->beamx    = g4minervaObj.protonX;
  //not-in-new//dk2nu->beamy    = g4minervaObj.protonY;
  //not-in-new//dk2nu->beamz    = g4minervaObj.protonZ;
  //not-in-new//dk2nu->beampx   = g4minervaObj.protonPx;
  //not-in-new//dk2nu->beampy   = g4minervaObj.protonPy;
  //not-in-new//dk2nu->beampz   = g4minervaObj.protonPz;

  //no-equiv//  dk2nu->tgptype  = g4minervaObj.tgptype;

  // now copy ancestor history
  if ( g4minervaObj.overflow ) dk2nu->flagbits |= bsim::kFlgOverflow;  // flag overflow
  int nmx = TMath::Min(g4minervaObj.ntrajectory,10);
  for (int ian=0; ian < nmx; ++ian) {

    const double mm2cm = 0.1;
    const double mev2gev = 0.001;

    bsim::Ancestor ancestor;
    /*
    dk2nu->apdg.push_back(g4minervaObj.pdg[ian]);
    dk2nu->trackid.push_back(g4minervaObj.trackId[ian]);
    dk2nu->parentid.push_back(g4minervaObj.parentId[ian]);

    dk2nu->startx.push_back(g4minervaObj.startx[ian]*mm2cm);
    dk2nu->starty.push_back(g4minervaObj.starty[ian]*mm2cm);
    dk2nu->startz.push_back(g4minervaObj.startz[ian]*mm2cm);
    dk2nu->stopx.push_back(g4minervaObj.stopx[ian]*mm2cm);
    dk2nu->stopy.push_back(g4minervaObj.stopy[ian]*mm2cm);
    dk2nu->stopz.push_back(g4minervaObj.stopz[ian]*mm2cm);

    dk2nu->startpx.push_back(g4minervaObj.startpx[ian]*mev2gev);
    dk2nu->startpy.push_back(g4minervaObj.startpy[ian]*mev2gev);
    dk2nu->startpz.push_back(g4minervaObj.startpz[ian]*mev2gev);
    dk2nu->stoppx.push_back(g4minervaObj.stoppx[ian]*mev2gev);
    dk2nu->stoppy.push_back(g4minervaObj.stoppy[ian]*mev2gev);
    dk2nu->stoppz.push_back(g4minervaObj.stoppz[ian]*mev2gev);

    dk2nu->pprodpx.push_back(g4minervaObj.pprodpx[ian]*mev2gev);
    dk2nu->pprodpy.push_back(g4minervaObj.pprodpy[ian]*mev2gev);
    dk2nu->pprodpz.push_back(g4minervaObj.pprodpz[ian]*mev2gev);

    // TString -> std::string via char*
    dk2nu->proc.push_back(g4minervaObj.proc[ian].Data());
    dk2nu->ivol.push_back(g4minervaObj.ivol[ian].Data());
    dk2nu->fvol.push_back(g4minervaObj.fvol[ian].Data());
    */

    ancestor.pdg = g4minervaObj.pdg[ian];
    ancestor.SetStartXYZT(g4minervaObj.startx[ian]*mm2cm,
                          g4minervaObj.starty[ian]*mm2cm,
                          g4minervaObj.startz[ian]*mm2cm,
                          0);  // don't have a time
    ancestor.SetStartP(g4minervaObj.startpx[ian]*mev2gev,
                       g4minervaObj.startpy[ian]*mev2gev,
                       g4minervaObj.startpz[ian]*mev2gev);
    ancestor.SetStopP(g4minervaObj.stoppx[ian]*mev2gev,
                      g4minervaObj.stoppy[ian]*mev2gev,
                      g4minervaObj.stoppz[ian]*mev2gev);
    ancestor.SetPProdP(g4minervaObj.pprodpx[ian]*mev2gev,
                       g4minervaObj.pprodpy[ian]*mev2gev,
                       g4minervaObj.pprodpz[ian]*mev2gev);

    // ancestor.polx = ancestor.poly = ancestor.polz = ?
    // ancestor.nucleau = ?
    ancestor.proc = g4minervaObj.proc[ian].Data();
    ancestor.ivol = g4minervaObj.ivol[ian].Data();
    ancestor.imat = g4minervaObj.imat[ian].Data();

    dk2nu->ancestor.push_back(ancestor);
  }

}

void g4minervaCrossChecks(const g4minerva& g4minervaObj)
{

  static bool first = true;
  static size_t indxMinosN = -1, indxMinosF = -1, indxNovaF = -1, indxMini = -1;
  static size_t ioldMinosN =  0, ioldMinosF =  0, ioldNovaF =  1, ioldMini =  8;
  if ( first ) {
    first = false;
    /// find indices for new weights so we can compare later
    indxMinosN = find_loc_index("MINOS NearDet");
    indxMinosF = find_loc_index("MINOS FarDet");
    //indxNovaF  = find_loc_index("NOvA FarDet");
    indxNovaF  = find_loc_index("Ash River per Minerva");
    indxMini   = find_loc_index("MiniBooNE");
    cout << "indx MINOS Near " << indxMinosN << " Far " << indxMinosF << endl;
    cout << "indx NOvA Far " << indxNovaF << endl;
    cout << "indx MiniBooNE " << indxMini << endl;

    /****  minerva g4numi/src/NumiDataInput.cc
     *  RWH: these appear to be in meters ... not sure about what
     *       the "NovaX" are meant to represent
     *
     //Near & Far Detector location
     nNear=11;//was 9 without the different energy for the ND positions.
     nFar=2;
     G4double xdetNear[]    = {   0     , 0.     , 7.     , 11.    ,
                                 14.    , 14.    , 14.    , 0.     , 
                                 25.84  , 4.8/2. , -4.8/2.         };
     G4double ydetNear[]    = {  0     , -3.    , -5.    , -5.    ,
                                -6.    , -3.    ,  0.    , 71.    ,
                                78.42  , 3.8/2. , -3.8/2.          };
     G4double zdetNear[]    = {  1040  , 1010.  , 975.   , 958.   ,
                                 940.  ,  840.  ,  740.  , 940.   ,
                                745.25 , 1040+16.6/2. , 1040-16.6/2.  };
     G4String detNameNear[] = {  "Near", "Nova1a","Nova1b","Nova1c",
                               "Nova2a", "Nova2b", "Nova3","MSB",
                               "MiniBooNE","Near +x +y +z","Near -x -y -z"};
     G4double xdetFar[]     = {0     , 28.81258   };
     G4double ydetFar[]     = {0     , 81.39258   };
     G4double zdetFar[]     = {735000, 811400     };
     G4String detNameFar[]  = {"Far" , "Ash River"};
     *
     */
  }
  // cross check location energy/weights
  // 
  histCompare(g4minervaObj.NenergyN[ioldMinosN],dk2nu->nuray[indxMinosN].E,   "MINOS near energy");
  histCompare(g4minervaObj.NWtNear[ioldMinosN], dk2nu->nuray[indxMinosN].wgt, "MINOS near wgt");
  histCompare(g4minervaObj.NenergyF[ioldMinosF],dk2nu->nuray[indxMinosF].E,   "MINOS far energy");
  histCompare(g4minervaObj.NWtFar[ioldMinosF],  dk2nu->nuray[indxMinosF].wgt, "MINOS far wgt");
  if ( indxNovaF > 0 ) {
    histCompare(g4minervaObj.NenergyF[ioldNovaF], dk2nu->nuray[indxNovaF].E,    "NOvA far energy");
    histCompare(g4minervaObj.NWtFar[ioldNovaF],   dk2nu->nuray[indxNovaF].wgt,  "NOvA far wgt");
  } else {
    static int nmsgnovaf = 10;
    if ( nmsgnovaf-- > 0 ) {
      cout << "## no relevant index for NOvA FarDet that will match minerva version of location" << endl;
    }
  }
  histCompare(g4minervaObj.NenergyN[ioldMini],  dk2nu->nuray[indxMini].E,     "MiniBooNE energy");
  histCompare(g4minervaObj.NWtNear[ioldMini],   dk2nu->nuray[indxMini].wgt,   "MiniBooNE wgt");

  static int nmsg = 0;
  double a = g4minervaObj.NWtNear[ioldMini];
  double b = dk2nu->nuray[indxMini].wgt;
  double myfdiff = ( (a+b) > 1.0e-30 ) ? TMath::Abs(2.0*(a-b)/(a+b)) : 0;
  if ( (nmsg < 1) &&  (myfdiff > 1.8) ) {
    ++nmsg;
    cout << "#### extreme case:  mini wtg " << a << " " << b << " fdiff " << myfdiff << endl;
    cout << "## g4 MinosN=" << ioldMinosN << " MinosF=" << ioldMinosF << " NovaF=" << ioldNovaF << " MiniBooNE=" << ioldMini << endl;
    cout << "g4minerva NenergyN ";
    for (int i=0; i<11; ++i) cout << ((i!=0&&i%4==0)?"\n                    ":" ") 
                                  << std::setw(12) << g4minervaObj.NenergyN[i];
    cout << endl;
    cout << "g4minerva NWtNear  ";
    for (int i=0; i<11; ++i) cout << ((i!=0&&i%4==0)?"\n                    ":" ")
                                  << std::setw(12) << g4minervaObj.NWtNear[i];
    cout << endl;
    cout << "g4minerva NenergyF ";
    for (int i=0; i<2; ++i) cout << " " << std::setw(12) << g4minervaObj.NenergyF[i];
    cout << endl;
    cout << "g4minerva NWtFar   ";
    for (int i=0; i<2; ++i) cout << " " << std::setw(12) << g4minervaObj.NWtFar[i];
    cout << endl;
    cout << "## dk MinosN=" << indxMinosN << " MinosF=" << indxMinosF << " NovaF=" << indxNovaF << " MiniBooNE=" << indxMini << endl;
    cout << *dk2nu << endl;
  }

  int nmx = TMath::Min(g4minervaObj.ntrajectory,10);
  for (int ian=0; ian < nmx; ++ian) {
    if ( ian <= 0 ) continue;
    // except for first particle start[i] should == stop[i-1]
    bool match = true;
    match = match && ( histCompare(g4minervaObj.startx[ian],g4minervaObj.stopx[ian-1],"startx") );
    match = match && ( histCompare(g4minervaObj.starty[ian],g4minervaObj.stopy[ian-1],"starty") );
    match = match && ( histCompare(g4minervaObj.startz[ian],g4minervaObj.stopz[ian-1],"startz") );
    match = match && ( g4minervaObj.ivol[ian] == g4minervaObj.fvol[ian-1] );
    if ( ! match ) {
      cout << "## ancestor " << ian << " didn't start where previous entry stopped" << endl;
    }

  }

}

