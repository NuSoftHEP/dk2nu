//____________________________________________________________________________
/*!

\class    genie::flux::GDk2NuFlux

\brief    An implementation of the GENIE GFluxI interface ("flux driver")
          encapsulating reading/processing the "dk2nu" tree structure.

\author   Robert Hatcher <rhatcher \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Nov 6, 2012

\cpright  Copyright (c) 2012, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GDK2NUFLUX_H_
#define _GDK2NUFLUX_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

// post-R3 this needs -I ${GENIE}/src/Framework as well as -I ${GENIE}/src
// #include "Conventions/GVersion.h"
// can't rely on that ... so we'll require -DGENIE_PRE_R3


#ifdef GENIE_PRE_R3
  #include "Conventions/GVersion.h"
  #include "EVGDrivers/GFluxI.h"
  #include "PDG/PDGUtils.h"
  #if __GENIE_RELEASE_CODE__ >= GRELCODE(2,9,0)
    #include "FluxDrivers/GFluxExposureI.h"
    #include "FluxDrivers/GFluxFileConfigI.h"
   #endif
#else
  #include "Framework/Conventions/GVersion.h"
  #include "Framework/EventGen/GFluxI.h"
  #include "Framework/ParticleData/PDGUtils.h"
  #include "Tools/Flux/GFluxExposureI.h"
  #include "Tools/Flux/GFluxFileConfigI.h"
#endif

class TFile;
class TChain;
class TTree;
class TBranch;

namespace bsim {
  class Dk2Nu;
  class DkMeta;
  class NuChoice;
}

namespace genie {
namespace flux  {

// #if __GENIE_RELEASE_CODE__ >= GRELCODE(2,9,0)
class GDk2NuFlux : public GFluxI,
    public GFluxExposureI, public GFluxFileConfigI {
// #else
// class GDk2NuFlux: public GFluxI {
// #endif

public :
  GDk2NuFlux();
 ~GDk2NuFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the NuMI neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void);
  double                 Weight        (void) { return  fWeight;              }
  const TLorentzVector & Momentum      (void);
  const TLorentzVector & Position      (void);
  bool                   End           (void) { return  fEnd;                 }
  long int               Index         (void) { return  fIEntry;              }
  void                   Clear            (Option_t * opt);
  void                   GenerateWeighted (bool gen_weighted);

  // Methods specific to this flux driver,
  // for configuration/initialization of the flux & event generation drivers
  // and and for passing-through flux information (e.g. neutrino parent decay
  // kinematics) not used by the generator but required by analyses/processing
  // further downstream

  //
  // information about or actions on current entry
  //
  const bsim::NuChoice &  GetNuChoice(void) { return *fCurNuChoice; };
  const bsim::Dk2Nu &     GetDk2Nu(void)    { return *fCurDk2Nu; };
  const bsim::DkMeta &    GetDkMeta(void)   { LoadDkMeta(); return *fCurDkMeta; };

  Long64_t GetEntryNumber() { return fIEntry; }   ///< index in chain

  double    GetDecayDist() const; ///< dist (user units) from dk to current pos
  void      MoveToZ0(double z0);  ///< move ray origin to user coord Z0

  //
  // information about the current state
  //
// #if __GENIE_RELEASE_CODE__ >= GRELCODE(2,9,0)
  virtual double GetTotalExposure() const;  // GFluxExposureI interface
// #else
//   // these are elminated with the introduction of GFluxFileConfigI
//   void                   SetXMLFile(std::string xmlbasename="GNuMIFlux.xml")
//                               { fXMLbasename = xmlbasename; }
//   std::string            GetXMLFile() const { return fXMLbasename; }
// #endif

  double    POT_curr(void);             ///< current average POT (RWH?)
  double    UsedPOTs(void) const;       ///< # of protons-on-target used
  long int  NFluxNeutrinos(void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double    SumWeight(void) const { return fSumWeight;  } ///< integrated weight for flux neutrinos looped so far

  void      PrintCurrent(void);         ///< print current entry from leaves
  void      PrintConfig();              ///< print the current configuration

  std::vector<std::string> GetFileList();  ///< list of files currently part of chain

  //
  // GFluxFileConfigI interface
  //
  virtual void  LoadBeamSimData(const std::vector<std::string>& filenames,
                                const std::string&         det_loc);
// #if __GENIE_RELEASE_CODE__ >= GRELCODE(2,9,0)
  using GFluxFileConfigI::LoadBeamSimData; // inherit the rest
  virtual void GetBranchInfo(std::vector<std::string>& branchNames,
                             std::vector<std::string>& branchClassNames,
                             std::vector<void**>&      branchObjPointers);
  virtual TTree* GetMetaDataTree();
// #else
//   void      LoadBeamSimData(std::set<std::string>  filenames, std::string det_loc);     ///< load root flux ntuple files and config
//   void      LoadBeamSimData(std::string filename, std::string det_loc);     ///< older (obsolete) single file version
//
// #endif

  //
  // configuration of GDk2NuFlux
  //
  void      SetTreeNames(std::string fname = "dk2nuTree",
                         std::string mname = "dkmetaTree")
  { fTreeNames[0] = fname; fTreeNames[1] = mname; }


  bool      LoadConfig(std::string cfg);      ///< load a named configuration

  void      SetMaxEnergy(double Ev);          ///< specify max flx nu energy

  void      SetGenWeighted(bool genwgt=false) { fGenWeighted = genwgt; } ///< toggle whether GenerateNext() returns weight=1 flux (initial default false)

  void      SetEntryReuse(long int nuse=1);   ///<  # of times to use entry before moving to next

  void      ScanForMaxWeight(void);           ///< scan for max flux weight (before generating unweighted flux neutrinos)
  void      SetMaxWgtScan(double fudge = 1.05, long int nentries = 2500000)      ///< configuration when estimating max weight
            { fMaxWgtFudge = fudge; fMaxWgtEntries = nentries; }
  void      SetMaxEFudge(double fudge = 1.05) ///< extra fudge factor in estimating maximum energy
            { fMaxEFudge = fudge; }

  void      SetMinMaxWeight(double minwgt)    ///< user floor on estimating maxweight (to avoid bumping)
            { fMinMaxWeight = minwgt; }
  void      SetMaxWeightFailModel(int i=0)    ///< 0=bump, 1=leave frozen, 2=abort
            { fMaxWgtFailModel = i; }


  void      SetApplyWindowTiltWeight(bool apply = true)           ///< apply wgt due to tilt of flux window relative to beam
            { fApplyTiltWeight = apply; }

  // GDk2NuFlux uses "cm" as the length unit consistently internally (this is
  // the length units used by the dk2nu ntuples).  User interactions
  // setting the beam-to-detector coordinate transform, flux window, and the
  // returned position might need to be in other units.  Use:
  //     double scale = genie::utils::units::UnitFromString("cm");
  // ( #include "Utils/UnitUtils.h for declaration )
  // to get the correct scale factor to pass in.  This should get set
  // FIRST before setting detector position/rotation

  void   SetLengthUnits(double user_units);  ///< Set units assumed by user
  double    LengthUnits(void) const;         ///< Return user units
  void   SetTimeUnits(double user_units);    ///< Set units assumed by user
  double    TimeUnits(void) const;           ///< Return user units

  // set relative orientation of user coords vs. beam system, i.e.
  //  x3_user = ( beamrot * x3_beam ) + x0beam_user
  //  p3_user =   beamrot * p3_beam

  ///< beam (0,0,0) relative to user frame, beam direction in user frame
  void      SetBeamRotation(TRotation beamrot);
  void      SetBeamCenter(TVector3 beam0);
  TRotation GetBeamRotation() const; ///< rotation to apply from beam->user
  TVector3  GetBeamCenter() const;   ///< beam origin in user frame

  // configure a flux window (or point) where E_nu and weight are evaluated

  // rwh: potential upgrade: allow flux window set/get in beam coords
  // as optional flag to *etFluxWindow
  bool      IsFluxSphere() const { return fIsSphere; }                     ///< flat window or sphere
  void      SetFluxWindow(TVector3  p1, TVector3  p2, TVector3  p3);       ///< 3 points define a plane (by default in user coordinates)
  void      GetFluxWindow(TVector3& p1, TVector3& p2, TVector3& p3) const; ///< 3 points define a plane in beam coordinate
  void      SetFluxSphere(TVector3  center, double  radius, bool inDetCoord=true);       ///< specification of a sphere
  void      GetFluxSphere(TVector3& center, double& radius, bool inDetCoord=true) const; ///< specification of a sphere

// #if __GENIE_RELEASE_CODE__ < GRELCODE(2,9,0)
//   // migrated to GFluxFileConfigI
//   void      SetUpstreamZ(double z0);
//   void      SetNumOfCycles(long int ncycle);
//   void      SetFluxParticles(const PDGCodeList & particles);
// #endif

  //
  // Actual coordinate transformations  b=beam, u=user (e.g. detector)
  //
  void      Beam2UserPos(const TLorentzVector& beamxyz,
                               TLorentzVector& usrxyz  ) const;
  void      Beam2UserDir(const TLorentzVector& beamdir,
                               TLorentzVector& usrdir  ) const;
  void      Beam2UserP4 (const TLorentzVector& beamp4,
                               TLorentzVector& usrp4   ) const;
  void      User2BeamPos(const TLorentzVector& usrxyz,
                               TLorentzVector& beamxyz ) const;
  void      User2BeamDir(const TLorentzVector& usrdir,
                               TLorentzVector& beamdir ) const;
  void      User2BeamP4 (const TLorentzVector& usrp4,
                               TLorentzVector& beamp4  ) const;

  TVector3  FluxWindowNormal() const { return (fIsSphere) ? SphereNormal() : fFluxWindowNormal; }
  TVector3  SphereNormal() const;

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);
  void AddFile               (TTree* fluxtree, TTree* metatree, std::string fname);
  void CalcEffPOTsPerNu      (void);
  void LoadDkMeta            (void);

  // Private data members
  //
  double         fMaxEv;          ///< maximum energy

  bool           fEnd;            ///< end condition reached

// #if __GENIE_RELEASE_CODE__ < GRELCODE(2,9,0)
//   // incorporated into GFluxFileConfigI
//   PDGCodeList *  fPdgCList;       ///< list of neutrino pdg-codes to generate
//   PDGCodeList *  fPdgCListRej;    ///< list of neutrino pdg-codes seen but rejected
//   std::string    fXMLbasename;         ///< XML filename for config data
//   double    fZ0;                  ///< configurable starting z position for each flux neutrino (in detector coord system)
//   long int  fNCycles;             ///< # times to cycle through the flux ntuple
//   long int  fICycle;              ///< current file cycle
// #endif

  std::vector<std::string> fNuFluxFilePatterns;   ///< (potentially wildcarded) path(s)

  std::string fTreeNames[2];      ///< pair of names "dk2nuTree", "dkmetaTree"
  TChain*   fNuFluxTree;          ///< TTree // REF ONLY!
  TChain*   fNuMetaTree;          ///< TTree // REF ONLY!

  bsim::Dk2Nu*     fCurDk2Nu;
  bsim::DkMeta*    fCurDkMeta;
  bsim::NuChoice*  fCurNuChoice;

  int       fNFiles;              ///< number of files in chain
  Long64_t  fNEntries;            ///< number of flux ntuple entries
  Long64_t  fIEntry;              ///< current flux ntuple entry
  Long64_t  fNuTot;               ///< cummulative # of entries (=fNEntries)
  Long64_t  fFilePOTs;            ///< # of protons-on-target represented by all files

  std:: map<int,int>  fJobToMetaIndex;  ///< quick lookup from job# to meta chain

  double    fWeight;              ///< current neutrino weight, =1 if generating unweighted entries
  double    fMaxWeight;           ///< max flux neutrino weight in input file
  double    fMinMaxWeight;        ///< user set lower limit on estimate
  double    fMaxWeightScan;       ///< initial estimate from scan
  double    fMaxWeightInit;       ///< max of scan & minmaxweight
  double    fMaxWeightMax;        ///< if "frozen" this is what bump would given
  double    fMaxWgtFudge;         ///< fudge factor for estimating max wgt
  long int  fMaxWgtEntries;       ///< # of entries in estimating max wgt
  double    fMaxEFudge;           ///< fudge factor for estmating max enu (0=> use fixed 120GeV)
  long int  fMaxWgtExceeded;      ///< track failures of estimate
  int       fMaxWgtFailModel;     ///< what to do ... 0=bump, 1=frozen, 2=abort

  long int  fNUse;                ///< how often to use same entry in a row
  long int  fIUse;                ///< current # of times an entry has been used
  double    fSumWeight;           ///< sum of weights for nus thrown so far
  long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far
  double    fEffPOTsPerNu;        ///< what a entry is worth ...
  double    fAccumPOTs;           ///< POTs used so far

  bool      fGenWeighted;         ///< does GenerateNext() give weights?
  bool      fApplyTiltWeight;     ///< wgt due to window normal not || beam
  bool      fDetLocIsSet;         ///< is a flux location (near/far) set?

  double           fLengthUnits;    ///< units for coord in user exchanges
  double           fTimeUnits;      ///< units for coord in user exchanges
  double           fLengthScaleB2U; ///< scale factor beam (cm) --> user
  double           fLengthScaleU2B; ///< scale factor beam user --> (cm)
  double           fTimeScaleB2U;   ///< scale factor beam (ns) --> user
  double           fTimeScaleU2B;   ///< scale factor beam user --> (ns)

  TLorentzVector   fBeamZero;       ///< beam origin in user coords
  TLorentzRotation fBeamRot;        ///< rotation applied beam --> user coord
  TLorentzRotation fBeamRotInv;

  bool             fIsSphere;             ///< doing this on a sphere rather than a flat window?
  TVector3         fFluxWindowPtUser[3];  ///<  user points of flux window
  TLorentzVector   fFluxWindowBase;       ///< base point for flux window - beam coord
  TLorentzVector   fFluxWindowDir1;       ///< extent for flux window (direction 1)
  TLorentzVector   fFluxWindowDir2;       ///< extent for flux window (direction 2)
  double           fFluxWindowLen1;
  double           fFluxWindowLen2;
  TVector3         fFluxWindowNormal;     ///< normal direction for flux window -- beam coord
  TVector3         fFluxSphereCenterUser; ///< center for flux sphere - user coords
  TVector3         fFluxSphereCenterBeam; ///< center for flux sphere - beam coords
  double           fFluxSphereRadius;     ///< radius for flux sphere

  TLorentzVector   fgX4dkvtx;             ///< decay 4-position beam coord

};

} // flux namespace
} // genie namespace

#endif // _GDK2NUFLUX_H_
