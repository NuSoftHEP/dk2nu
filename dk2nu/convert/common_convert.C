/// common code for use in the ntuple converter code
///
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \created   2012-11-07
/// \version $Id: common_convert.C,v 1.3 2012-11-15 09:09:26 rhatcher Exp $
///==========================================================================

#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>
  
#include <float.h> // FLT_EPSILON and DBL_EPSILON definitions used for
                   // floating point comparisons

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1D.h"

// dk2nu headers
#include "tree/dk2nu.h"
#include "tree/dkmeta.h"
/// functions for reading location text file and calculating
/// energy and weight vectors for locations
#include "tree/readWeightLocations.h"
#include "tree/calcLocationWeights.h"

// some globals
TRandom3* rndm            = 0;
bsim::Dk2Nu*  dk2nu       = 0;
bsim::DkMeta* dkmeta      = 0;
TFile*        treeFile    = 0;
TTree*        dk2nuTree   = 0;
TTree*        dkmetaTree  = 0;
int           myjob       = 0;
int           pots        = 0;

double    frac_diff_tolerance = 5.0e-4;
double    obs_frac_diff_max  = 0;

bool histCompare(double a, double b, string s);  // forward declaration

//____________________________________________________________________________

void ConvertInit(std::string locfilename = "$(DK2NU)/etc/locations.txt")
{
  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::NumiAnalysis() in g4numi
  ///
  ///-----------------------------------------------------------------------

  // create the objects used for writing entries
  dk2nu  = new bsim::Dk2Nu;
  dkmeta = new bsim::DkMeta;

  rndm = new TRandom3();

  // read the text file for locations, fill the dkmeta object
  /*bsim::*/readWeightLocations(locfilename,dkmeta);

  // print out what we have for locations
  /*bsim::*/printWeightLocations(dkmeta);

}

void ConvertBookNtuple(std::string ofname = "test_dk2nu.root")
{
  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::book() in g4numi
  ///
  ///-----------------------------------------------------------------------

  // create file, book tree, set branch address to created object 
  treeFile = new TFile(ofname.c_str(),"RECREATE");

  dk2nuTree = new TTree("dk2nuTree","neutrino ntuple");
  dk2nuTree->Branch("dk2nu","bsim::Dk2Nu",&dk2nu,32000,1);

  dkmetaTree  = new TTree("dkmetaTree","neutrino ntuple metadata");
  dkmetaTree->Branch("dkmeta","bsim::DkMeta",&dkmeta,32000,1);
}

void ConvertFinish()
{
  ///-----------------------------------------------------------------------
  ///
  ///  equivalent to NumiAnalysis::finish() in g4numi
  ///
  ///-----------------------------------------------------------------------

  /// fill the rest of the metadata (locations filled above)
  //no! would clear location info // dkmeta->clear();
  dkmeta->job  = myjob;  // needs to match the value in each dk2nu entry
  dkmeta->pots = pots;   // ntuple represents this many protons-on-target 
  //leave to user// dkmeta->beamsim = "common_convert.C";
  //leave to user// dkmeta->physics = "bogus";

  //dkmeta->vintnames.push_back("mytemp_42");
  //dkmeta->vintnames.push_back("mytemp_ipot");
  // push entry out to meta-data tree
  dkmetaTree->Fill();

  // finish and clean-up
  treeFile->cd();
  dk2nuTree->Write();
  dkmetaTree->Write();

  treeFile->cd();  // be here so any booked histograms get created inside output file
  treeFile->mkdir("zzz_diff_hists");
  treeFile->cd("zzz_diff_hists");
  histCompare(0,0,"WriteHist");

  treeFile->Close();

  delete treeFile; treeFile=0;
  dk2nuTree=0;
  dkmetaTree=0;
}
//____________________________________________________________________________

const int kPdgNuE              =  12;   //
const int kPdgAntiNuE          = -12;   //
const int kPdgNuMu             =  14;   //
const int kPdgAntiNuMu         = -14;   // 
const int kPdgNuTau            =  16;   //
const int kPdgAntiNuTau        = -16;   //

const int kPdgElectron         =  11;   //
const int kPdgPositron         = -11;   //
const int kPdgMuon             =  13;   //
const int kPdgAntiMuon         = -13;   //
const int kPdgTau              =  15;   //
const int kPdgAntiTau          = -15;   //

const int kPdgUQuark           =   2;   //
const int kPdgAntiUQuark       =  -2;   //
const int kPdgDQuark           =   1;   //
const int kPdgAntiDQuark       =  -1;   //
const int kPdgSQuark           =   3;   //
const int kPdgAntiSQuark       =  -3;   //
const int kPdgCQuark           =   4;   //
const int kPdgAntiCQuark       =  -4;   //
const int kPdgBQuark           =   5;   //
const int kPdgAntiBQuark       =  -5;   //
const int kPdgTQuark           =   6;   //
const int kPdgAntiTQuark       =  -6;   //

const int kPdgDDDiquarkS1      =  1103; // dd, spin = 1 
const int kPdgUDDiquarkS0      =  2101; // ud, spin = 0 
const int kPdgUDDiquarkS1      =  2103; // ud, spin = 1 
const int kPdgUUDiquarkS1      =  2203; // uu, spin = 1 
const int kPdgSDDiquarkS0      =  3101; // sd, spin = 0
const int kPdgSDDiquarkS1      =  3103; // sd, spin = 1
const int kPdgSUDiquarkS0      =  3201; // su, spin = 0
const int kPdgSUDiquarkS1      =  3203; // su, spin = 1
const int kPdgSSDiquarkS1      =  3303; // ss, spin = 1

const int kPdgProton           =  2212; //
const int kPdgAntiProton       = -2212; //
const int kPdgNeutron          =  2112; //
const int kPdgAntiNeutron      = -2112; //
const int kPdgLambda           =  3122; // Lambda
const int kPdgAntiLambda       = -3122; // \bar{Lambda}
const int kPdgSigmaP           =  3222; // Sigma+
const int kPdgSigma0           =  3212; // Sigma0
const int kPdgSigmaM           =  3112; // Sigma-
const int kPdgAntiSigmaP       = -3222; // \bar{Sigma+}
const int kPdgAntiSigma0       = -3212; // \bar{Sigma0}
const int kPdgAntiSigmaM       = -3112; // \bar{Sigma-}
const int kPdgXi0              =  3322; // Xi0
const int kPdgXiM              =  3312; // Xi-
const int kPdgAntiXi0          = -3322; // \bar{Xi0}
const int kPdgAntiXiP          = -3312; // \bar{Xi+}
const int kPdgOmegaM           =  3334; // Omega-
const int kPdgAntiOmegaP       = -3334; // \bar{Omega+}
const int kPdgLambdaPc         =  4122; // Lambda+_{c}
const int kPdgSigma0c          =  4112; // Sigma0_{c}
const int kPdgSigmaPc          =  4212; // Sigma+_{c}
const int kPdgSigmaPPc         =  4222; // Sigma++_{c}

const int kPdgP33m1232_DeltaM  =  1114; // P33(1232) Delta-
const int kPdgP33m1232_Delta0  =  2114; // P33(1232) Delta0
const int kPdgP33m1232_DeltaP  =  2214; // P33(1232) Delta+
const int kPdgP33m1232_DeltaPP =  2224; // P33(1232) Delta++
const int kPdgS11m1535_N0      = 22112; // S11(1535) N0
const int kPdgS11m1535_NP      = 22212; // S11(1535) N+
const int kPdgD13m1520_N0      =  1214; // D13(1520) N0
const int kPdgD13m1520_NP      =  2124; // D13(1520) N+
const int kPdgS11m1650_N0      = 32112; // S11(1650) N0
const int kPdgS11m1650_NP      = 32212; // S11(1650) N+
const int kPdgD13m1700_N0      = 21214; // D13(1700) N0
const int kPdgD13m1700_NP      = 22124; // D13(1700) N+
const int kPdgD15m1675_N0      =  2116; // D15(1675) N0
const int kPdgD15m1675_NP      =  2216; // D15(1675) N+
const int kPdgS31m1620_DeltaM  =  1112; // S31(1620) Delta-
const int kPdgS31m1620_Delta0  =  1212; // S31(1620) Delta0
const int kPdgS31m1620_DeltaP  =  2122; // S31(1620) Delta+
const int kPdgS31m1620_DeltaPP =  2222; // S31(1620) Delta++
const int kPdgD33m1700_DeltaM  = 11114; // D33(1700) Delta-
const int kPdgD33m1700_Delta0  = 12114; // D33(1700) Delta0
const int kPdgD33m1700_DeltaP  = 12214; // D33(1700) Delta+
const int kPdgD33m1700_DeltaPP = 12224; // D33(1700) Delta++
const int kPdgP11m1440_N0      = 12112; // P11(1440) N0
const int kPdgP11m1440_NP      = 12212; // P11(1440) N+
const int kPdgP13m1720_N0      = 31214; // P13(1720) N0
const int kPdgP13m1720_NP      = 32124; // P13(1720) N+
const int kPdgF15m1680_N0      = 12116; // F15(1680) N0
const int kPdgF15m1680_NP      = 12216; // F15(1680) N+
const int kPdgP31m1910_DeltaM  = 21112; // P31(1910) Delta-
const int kPdgP31m1910_Delta0  = 21212; // P31(1910) Delta0
const int kPdgP31m1910_DeltaP  = 22122; // P31(1910) Delta+
const int kPdgP31m1910_DeltaPP = 22222; // P31(1910) Delta++
const int kPdgP33m1920_DeltaM  = 21114; // P33(1920) Delta-
const int kPdgP33m1920_Delta0  = 22114; // P33(1920) Delta0
const int kPdgP33m1920_DeltaP  = 22214; // P33(1920) Delta+
const int kPdgP33m1920_DeltaPP = 22224; // P33(1920) Delta++
const int kPdgF35m1905_DeltaM  =  1116; // F35(1905) Delta-
const int kPdgF35m1905_Delta0  =  1216; // F35(1905) Delta0
const int kPdgF35m1905_DeltaP  =  2126; // F35(1905) Delta+
const int kPdgF35m1905_DeltaPP =  2226; // F35(1905) Delta++
const int kPdgF37m1950_DeltaM  =  1118; // F37(1950) Delta-
const int kPdgF37m1950_Delta0  =  2118; // F37(1950) Delta0
const int kPdgF37m1950_DeltaP  =  2218; // F37(1950) Delta+
const int kPdgF37m1950_DeltaPP =  2228; // F37(1950) Delta++
const int kPdgP11m1710_N0      = 42112; // P11(1710) N0
const int kPdgP11m1710_NP      = 42212; // P11(1710) N+

const int kPdgPiP              =   211; // pi+
const int kPdgPiM              =  -211; // pi-
const int kPdgPi0              =   111; // pi0
const int kPdgEta              =   221; // eta
const int kPdgEtaPrm           =   331; // eta' (prime)
const int kPdgEtac             =   441; // eta_{c}
const int kPdgEtab             =   551; // eta_{b}
const int kPdgRhoP             =   213; // rho+
const int kPdgRhoM             =  -213; // rho-
const int kPdgRho0             =   113; // rho0
const int kPdgomega            =   223; // omega (the meson, not Omega the baryon)
const int kPdgPhi              =   333; // phi
const int kPdgJpsi             =   443; // J/psi
const int kPdgY                =   553; // Y
const int kPdgKP               =   321; // K+
const int kPdgKM               =  -321; // K-
const int kPdgK0               =   311; // K0
const int kPdgAntiK0           =  -311; // \bar{K0}
const int kPdgK0L              =   130; // K0_{long}
const int kPdgK0S              =   310; // K0_{short}
const int kPdgDP               =   411; // D+
const int kPdgDM               =  -411; // D-
const int kPdgD0               =   421; // D0
const int kPdgAntiD0           =  -421; // \bar{D0}
const int kPdgDPs              =   431; // D+_{s}
const int kPdgDMs              =  -431; // D-_{s}

const int kPdgGluon            =    21; // gluon
const int kPdgGamma            =    22; // photon
const int kPdgZ0               =    23; // Z
const int kPdgWP               =    24; // W+
const int kPdgWM               =   -24; // W-

// Note: PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z)
// PDG2006 convention: 10LZZZAAAI 
// Define names for some commonly used nuclear PDG codes:
const int kPdgTgtFreeP     = 1000010010;
const int kPdgTgtFreeN     = 1000000010;
const int kPdgTgtDeuterium = 1000010020;
const int kPdgTgtC12       = 1000060120;
const int kPdgTgtO16       = 1000080160;
const int kPdgTgtFe56      = 1000260560;

// PDG codes for GENIE special particles
const int kPdgHadronicSyst = 2000000001; // dis hadronic system before hadronization
const int kPdgHadronicBlob = 2000000002; // unmodelled fraction of the hadronic system
const int kPdgBindino      = 2000000101; // binding energy subtracted from f/s nucleons
const int kPdgCoulobtron   = 2000000102; // coulomb energy subtracted from f/s leptons
const int kPdgClusterNN    = 2000000200; // a nn cluster within a nucleus
const int kPdgClusterNP    = 2000000201; // a np cluster within a nucleus
const int kPdgClusterPP    = 2000000202; // a pp cluster within a nucleus

// PDG codes for PYTHIA/JETSET special particles
const int kPdgCluster      = 91; 
const int kPdgString       = 92; 
const int kPdgIndep        = 93; 

//____________________________________________________________________________
int ConvertGeantToPdg(int geant_code, std::string tag="?")
{
  if (geant_code ==  3) return kPdgElectron;     //    11 / e-
  if (geant_code ==  2) return kPdgPositron;     //   -11 / e+
  if (geant_code ==  6) return kPdgMuon;         //    13 / mu-
  if (geant_code ==  5) return kPdgAntiMuon;     //   -13 / mu+             
  if (geant_code == 34) return kPdgTau;          //    15 / tau-
  if (geant_code == 33) return kPdgAntiTau;      //   -15 / tau+              
  if (geant_code ==  8) return kPdgPiP;          //   211 / pi+
  if (geant_code ==  9) return kPdgPiM;          //  -211 / pi-
  if (geant_code ==  7) return kPdgPi0;          //   111 / pi0
  if (geant_code == 17) return kPdgEta;          //   221 / eta
  if (geant_code == 11) return kPdgKP;           //   321 / K+
  if (geant_code == 12) return kPdgKM;           //  -321 / K-
  if (geant_code == 10) return kPdgK0L;          //   130 / K0_{long}
  if (geant_code == 16) return kPdgK0S;          //   310 / K0_{short}
  if (geant_code == 35) return kPdgDP;           //   411 / D+
  if (geant_code == 36) return kPdgDM;           //  -411 / D-
  if (geant_code == 37) return kPdgD0;           //   421 / D0
  if (geant_code == 38) return kPdgAntiD0;       //  -421 / \bar{D0}
  if (geant_code == 39) return kPdgDPs;          //   431 / D+_{s}
  if (geant_code == 40) return kPdgDMs;          //  -431 / D-_{s}
  if (geant_code ==  1) return kPdgGamma;        //    22 / photon
  if (geant_code == 44) return kPdgZ0;           //    23 / Z
  if (geant_code == 42) return kPdgWP;           //    24 / W+
  if (geant_code == 43) return kPdgWM;           //   -24 / W-
  if (geant_code == 14) return kPdgProton;       //  2212
  if (geant_code == 15) return kPdgAntiProton;   // -2212
  if (geant_code == 13) return kPdgNeutron;      //  2112
  if (geant_code == 25) return kPdgAntiNeutron;  // -2112
  if (geant_code == 18) return kPdgLambda;       //  3122 / Lambda
  if (geant_code == 26) return kPdgAntiLambda;   // -3122 / \bar{Lambda}
  if (geant_code == 19) return kPdgSigmaP;       //  3222 / Sigma+
  if (geant_code == 20) return kPdgSigma0;       //  3212 / Sigma0
  if (geant_code == 21) return kPdgSigmaM;       //  3112 / Sigma-
  if (geant_code == 29) return kPdgAntiSigmaP;   // -3112 / \bar{Sigma+}
  if (geant_code == 28) return kPdgAntiSigma0;   // -3212 / \bar{Sigma0}
  if (geant_code == 27) return kPdgAntiSigmaM;   // -3112 / \bar{Sigma-}
  if (geant_code == 22) return kPdgXi0;          //  3322 / Xi0
  if (geant_code == 23) return kPdgXiM;          //  3312 / Xi-
  if (geant_code == 30) return kPdgAntiXi0;      // -3322 / \bar{Xi0}
  if (geant_code == 31) return kPdgAntiXiP;      // -3312 / \bar{Xi+}
  if (geant_code == 24) return kPdgOmegaM;       //  3334 / Omega-
  if (geant_code == 32) return kPdgAntiOmegaP;   // -3334 / \bar{Omega+}

  // some rare Geant3 codes that don't really need definitions in PDGCodes.h
  const int kPdgDeuteron = 1000010020; // pdg::IonPdgCode(2,1);
  const int kPdgTritium  = 1000010030; // pdg::IonPdgCode(3,1);
  const int kPdgAlpha    = 1000020040; // pdg::IonPdgCode(4,2);
  const int kPdgHe3      = 1000020030; // pdg::IonPdgCode(3,2);
  if (geant_code == 45) return kPdgDeuteron;
  if (geant_code == 46) return kPdgTritium;
  if (geant_code == 47) return kPdgAlpha;
  if (geant_code == 49) return kPdgHe3;

  if (geant_code != 0 ) 
    std::cerr << "## Can not convert geant code: " << geant_code
              << " to PDG  (" << tag << ")" << std::endl;
  return 0;
}
//____________________________________________________________________________

int Convert5xToPdg(int old_ntype)
{
  // NuMI ntuples have an odd neutrino "geant" convention
  switch ( old_ntype ) {
  case 56: return  14; break;  // kPdgNuMu;     break;
  case 55: return -14; break;  // kPdgAntiNuMu; break;
  case 53: return  12; break;  // kPdgNuE;      break;
  case 52: return -12; break;  // kPdgAntiNuE;  break;
  default:
    std::cerr << "Convert5xToPdg saw ntype " << old_ntype 
              << " -- unknown " << std::endl;
    assert(0);
  }
  return 0;
}

//____________________________________________________________________________

size_t find_loc_index(string match)
{
  // find the index in the array of locations
  for ( size_t i=0; i < dkmeta->location.size(); ++i ) {
    //cout << "looking for \"" << match << "\" vs \""
    //     << dkmeta->nameloc[i] << "\"" << endl;
    if ( dkmeta->location[i].name == match ) return i;
  }
  //cout << "failed to find match" << endl;
  return -1;
}
//____________________________________________________________________________

string construct_outfilename(string infilename)
{
  string ofname = infilename;
  size_t dot = ofname.find(".root");
  ofname.insert(dot,"_to_dk2nu");
  size_t slash = ofname.find_last_of("/");
  if ( slash != std::string::npos ) {
    ofname.erase(0,slash+1);
  }
  return ofname;
}

//____________________________________________________________________________

double estimate_pots(int highest_potnum)
{
  // looks like low counts are due to "evtno" not including                     
  // protons that miss the actual target (hit baffle, etc)                      
  // this will vary with beam conditions parameters                             
  // so we should round way up, those generating flugg files                    
  // aren't going to quantize less than 1000                                    
  // though 500 would probably be okay, 100 would not.                          
  // also can't use _last_ potnum because muons decays don't
  // have theirs set
  const Int_t    nquant = 1000; // 500;  // 100                                 
  const Double_t rquant = nquant;

  Int_t estimate = (TMath::FloorNint((highest_potnum-1)/rquant)+1)*nquant;
  return estimate;
}

//____________________________________________________________________________
bool isCloseEnough(double a, double b, double& difference, double& fdiff,
                   int& nmsg, int mxmsg, string s)
{
  // For float or double, need to determine equality while 
  // accounting for some degree of imprecision
  // The method applied here is adopted from Theodore C. Belding's
  // fcmp implementation of Donald Knuth's algorithm. 
  double value[2] = {a,b};
  //double epsilon = 10.*FLT_EPSILON;
  double epsilon = 10.*FLT_EPSILON;
  //      if ( isDouble[0] || isDouble[1] ) epsilon = DBL_EPSILON;
  int exponent;
  frexp(fabs(value[0]) > fabs(value[1]) ? value[0] : value[1], &exponent);
  double delta = ldexp(epsilon,exponent);
  //double 
  difference = value[0] - value[1];
  if ( difference > delta || difference < -delta ) {
    double sum  = a+b;
    if ( TMath::Abs(sum) < 1.0e-30 ) sum = 1.0e-30;  // protect divide-by-0
    //double
    fdiff = TMath::Abs(2.0*(a-b)/(a+b));
    if ( fdiff > obs_frac_diff_max ) obs_frac_diff_max = fdiff;
    if ( fdiff < frac_diff_tolerance ) return true; // we'll allow it anyway
    ++nmsg;
    if ( nmsg <= mxmsg ) {
      cout << "## \"" << s << "\" " << a << " " << b 
        //<< " delta " << delta
           << " diff " << difference
           << " fdiff " << fdiff << endl;
      //cout << "## ndecay (dkproc) " << dk2nu->decay.ndecay << endl;
      if ( nmsg == mxmsg ) 
        cout << "## last message for tag \"" << s << "\"" << endl;
    }
    return false;  // nope
  }
  return true; // these are close
}

bool histCompare(double a, double b, string s)
{
  // booking histograms should end up in the output file

  static map<std::string,int> nmsg;
  double diff, frac;
  const int mxmsg = 10, nbins = 100;
  static map<std::string,TH1D*> diffmap;
  static map<std::string,TH1D*> fracmap;
  ///  static bool first = true;
  //if ( first ) {
  //  first = false;
  //  TH1::SetDefaultSumw2(true);
  //}
  if ( s == "WriteHist" ) { 
    map<std::string,TH1D*>::iterator itr;
    for ( itr=diffmap.begin(); itr != diffmap.end(); ++itr ) (*itr).second->Write();
    for ( itr=fracmap.begin(); itr != fracmap.end(); ++itr ) (*itr).second->Write();
    return true;
  }

  bool close = isCloseEnough(a,b,diff,frac,nmsg[s],mxmsg,s);
  TH1D* diffh = diffmap[s];
  TH1D* frach = fracmap[s];
  if ( ! diffh  ) {
    string sdh = "difference " + s;
    diffh = new TH1D(sdh.c_str(),sdh.c_str(),nbins,1,-1); // autoadjust limits by reversed values
    diffh->Sumw2();
    diffmap[s] = diffh;
  }
  if ( ! frach  ) {
    string sfh = "frac diff " + s;
    frach = new TH1D(sfh.c_str(),sfh.c_str(),nbins,1,-1); // autoadjust limits by reversed values
    frach->Sumw2();
    fracmap[s] = frach;
  }
  diffh->Fill(diff);
  frach->Fill(frac);

  return close;
}

//____________________________________________________________________________

