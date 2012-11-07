//____________________________________________________________________________
/*!

\class    genie::flux::GDk2NuFlux

\brief    A GENIE flux driver encapsulating the "dk2nu" neutrino flux.
          It reads-in the official unified neutrino flux ntuples.

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

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#include "EVGDrivers/GFluxI.h"
#include "PDG/PDGUtils.h"

class TFile;
class TChain;
class TTree;
class TBranch;

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

namespace genie {
namespace flux  {

/// GDk2NuFluxPassThroughInfo:
/// =========================
/// A small persistable C-struct -like class that mirrors (some of) the 
/// structure of the gnumi ntuples.  This can then be stored as an extra 
/// branch of the output event tree -alongside with the generated event 
/// branch- for use further upstream in the analysis chain - e.g. beam 
/// reweighting etc.
/// To do future x-y reweighting users must retain the info found in:
//     Ntype   Vx      Vy      Vz      
//     pdPx    pdPy    pdPz    
//     ppdxdz  ppdydz  pppz    ppenergy ptype
//     muparpx muparpy muparpz mupare   Necm
//     Nimpwt  
///
class GDk2NuFluxPassThroughInfo: public TObject {
public:
   GDk2NuFluxPassThroughInfo();
   GDk2NuFluxPassThroughInfo(const dk2nu*);
   /* allow default copy constructor ... for now nothing special
   GDk2NuFluxPassThroughInfo(const GDk2NuFluxPassThroughInfo & info);
   */
   virtual ~GDk2NuFluxPassThroughInfo() { };

   void ResetCopy();     // reset portion copied from ntuple
   void ResetCurrent();  // reset generated xy positioned info
   void Print(const Option_t* opt = "") const;

   friend ostream & operator << (ostream & stream, const GDk2NuFluxPassThroughInfo & info);

   // values from ntuple
   dk2nu          dk2nuObj;

   // Values for GDk2NuFlux chosen x-y-z position, not from flux ntuple
   int            fgPdgC;   ///< generated nu pdg-code
   double         fgXYWgt;  ///< generated nu x-y weight
                            ///   not the same as GDk2NuFlux::Weight()
                            ///   which include importance wgt and deweighting
   TLorentzVector fgP4;     ///< generated nu 4-momentum beam coord
   TLorentzVector fgX4;     ///< generated nu 4-position beam coord
   TLorentzVector fgP4User; ///< generated nu 4-momentum user coord
   TLorentzVector fgX4User; ///< generated nu 4-position user coord

ClassDef(GDk2NuFluxPassThroughInfo,1)
};

/// GDk2NuFlux:
/// ==========
/// An implementation of the GFluxI interface that provides flux
/// from the unified dk2nu ntuples
///
class GDk2NuFlux: public GFluxI {

public :
  GDk2NuFlux();
 ~GDk2NuFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the NuMI neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fCurEntry->fgPdgC;    }
  double                 Weight        (void) { return  fWeight;              }
  const TLorentzVector & Momentum      (void) { return  fCurEntry->fgP4User;  }
  const TLorentzVector & Position      (void) { return  fCurEntry->fgX4User;  }
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
  const GDk2NuFluxPassThroughInfo &
     PassThroughInfo(void) { return *fCurEntry; } ///< GDk2NuFluxPassThroughInfo
  Long64_t GetEntryNumber() { return fIEntry; }   ///< index in chain

  double    GetDecayDist() const; ///< dist (user units) from dk to current pos
  void      MoveToZ0(double z0);  ///< move ray origin to user coord Z0

  //
  // information about the current state
  //
  double    POT_curr(void);             ///< current average POT (RWH?)
  double    UsedPOTs(void) const;       ///< # of protons-on-target used
  long int  NFluxNeutrinos(void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double    SumWeight(void) const { return fSumWeight;  } ///< integrated weight for flux neutrinos looped so far

  void      PrintCurrent(void);         ///< print current entry from leaves
  void      PrintConfig();              ///< print the current configuration

  std::vector<std::string> GetFileList();  ///< list of files currently part of chain

  //
  // configuration of GDk2NuFlux
  //
  void      SetXMLFile(string xmlbasename="GDk2NuFlux.xml") { fXMLbasename = xmlbasename; }  ///< set the name of the file that hold XML config param_sets
  std::string GetXMLFile() const { return fXMLbasename; }  ///< return the name of the file that hold XML config param_sets

  void      LoadBeamSimData(std::vector<string> filenames, string det_loc);     ///< load root flux ntuple files and config
  void      LoadBeamSimData(std::set<string>    filenames, string det_loc);     ///< load root flux ntuple files and config
  void      LoadBeamSimData(string filename, string det_loc);     ///< older (obsolete) single file version

  bool      LoadConfig(string cfg);                               ///< load a named configuration
  void      SetFluxParticles(const PDGCodeList & particles);      ///< specify list of flux neutrino species
  void      SetMaxEnergy(double Ev);                              ///< specify maximum flx neutrino energy

  void      SetGenWeighted(bool genwgt=false) { fGenWeighted = genwgt; } ///< toggle whether GenerateNext() returns weight=1 flux (initial default false)

  void      SetNumOfCycles(long int ncycle);                      ///< set how many times to cycle through the ntuple (default: 1 / n=0 means 'infinite')
  void      SetEntryReuse(long int nuse=1);                       ///<  # of times to use entry before moving to next

  void      ScanForMaxWeight(void);                               ///< scan for max flux weight (before generating unweighted flux neutrinos)
  void      SetMaxWgtScan(double fudge = 1.05, long int nentries = 2500000)      ///< configuration when estimating max weight
            { fMaxWgtFudge = fudge; fMaxWgtEntries = nentries; }
  void      SetMaxEFudge(double fudge = 1.05)                  ///< extra fudge factor in estimating maximum energy
            { fMaxEFudge = fudge; }

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
  void      SetFluxWindow(TVector3  p1, TVector3  p2, TVector3  p3); ///< 3 points define a plane (by default in user coordinates)
  void      GetFluxWindow(TVector3& p1, TVector3& p2, TVector3& p3) const; ///< 3 points define a plane in beam coordinate 

  void      SetUpstreamZ(double z0);                           ///< set flux neutrino initial z position (upstream of the detector) pushed back from the flux window

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

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);
  void AddFile               (string fname);
  void AddTreeFile           (TTree* tree, string fname);
  void CalcEffPOTsPerNu      (void);
  
  // Private data members
  //
  double         fMaxEv;          ///< maximum energy
  PDGCodeList *  fPdgCList;       ///< list of neutrino pdg-codes to generate
  PDGCodeList *  fPdgCListRej;    ///< list of neutrino pdg-codes seen but rejected
  bool           fEnd;            ///< end condition reached

  string    fXMLbasename;         ///< XML filename for config data
  std::vector<string> fNuFluxFilePatterns;   ///< (potentially wildcarded) path(s)

  TChain*   fNuFluxTree;          ///< TTree // REF ONLY!

  int       fNFiles;              ///< number of files in chain
  Long64_t  fNEntries;            ///< number of flux ntuple entries
  Long64_t  fIEntry;              ///< current flux ntuple entry
  Long64_t  fNuTot;               ///< cummulative # of entries (=fNEntries)
  Long64_t  fFilePOTs;            ///< # of protons-on-target represented by all files

  double    fWeight;              ///< current neutrino weight, =1 if generating unweighted entries
  double    fMaxWeight;           ///< max flux neutrino weight in input file
  double    fMaxWgtFudge;         ///< fudge factor for estimating max wgt
  long int  fMaxWgtEntries;       ///< # of entries in estimating max wgt
  double    fMaxEFudge;           ///< fudge factor for estmating max enu (0=> use fixed 120GeV)

  long int  fNCycles;             ///< # times to cycle through the flux ntuple
  long int  fICycle;              ///< current file cycle
  long int  fNUse;                ///< how often to use same entry in a row
  long int  fIUse;                ///< current # of times an entry has been used
  double    fSumWeight;           ///< sum of weights for nus thrown so far
  long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far
  double    fEffPOTsPerNu;        ///< what a entry is worth ...
  double    fAccumPOTs;           ///< POTs used so far

  bool      fGenWeighted;         ///< does GenerateNext() give weights?
  bool      fDetLocIsSet;         ///< is a flux location (near/far) set?
  
  double           fLengthUnits;    ///< units for coord in user exchanges
  double           fLengthScaleB2U; ///< scale factor beam (cm) --> user
  double           fLengthScaleU2B; ///< scale factor beam user --> (cm)

  TLorentzVector   fBeamZero;       ///< beam origin in user coords
  TLorentzRotation fBeamRot;        ///< rotation applied beam --> user coord
  TLorentzRotation fBeamRotInv;
  double           fZ0;             ///< configurable starting z position for each flux neutrino (in detector coord system)

  TVector3         fFluxWindowPtUser[3]; ///<  user points of flux window
  TLorentzVector   fFluxWindowBase; ///< base point for flux window - beam coord
  TLorentzVector   fFluxWindowDir1; ///< extent for flux window (direction 1)
  TLorentzVector   fFluxWindowDir2; ///< extent for flux window (direction 2)
  double           fFluxWindowLen1;
  double           fFluxWindowLen2;

  TLorentzVector   fgX4dkvtx;       ///< decay 4-position beam coord

  GDk2NuFluxPassThroughInfo* fCurEntry;  ///< copy of current ntuple entry info (owned structure)

};

} // flux namespace
} // genie namespace

#endif // _GDK2NUFLUX_H_
