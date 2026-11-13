//____________________________________________________________________________
/*
 Copyright (c) 2012, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.


*/
//____________________________________________________________________________

//#define __GENIE_LOW_LEVEL_MESG_ENABLED__

#include <cstdlib>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>
#include <climits>
#include <iostream>
#include <iomanip>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

// ROOT headers
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TBranch.h>
#include <TTree.h>
#include <TSystem.h>
#include <TStopwatch.h>

// GENIE R-3 (or later) headers
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/XmlParserUtils.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Numerical/MathUtils.h"

// Dk2Nu headers
#include "genie/GDk2NuFlux.h"
#include "tree/dk2nu.h"
#include "tree/dkmeta.h"
#include "tree/NuChoice.h"
#include "tree/calcLocationWeights.h"

// need GDk2NuFlux header to register w/ factory
#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GDk2NuFlux,genie::flux::GDk2NuFlux)

#include <vector>
#include <algorithm>
#include <iomanip>
#include "TRegexp.h"
#include "TString.h"

using namespace genie;
using namespace genie::flux;

// declaration of helper class
namespace genie {
  namespace flux  {
    class GDk2NuFluxXMLHelper {
    public:
      GDk2NuFluxXMLHelper(GDk2NuFlux* dk2nuFlux)
        : fVerbose(0), fGDk2NuFlux(dk2nuFlux) { ; }
      ~GDk2NuFluxXMLHelper() { ; }
      bool LoadConfig(std::string cfg);

      // these should go in a more general package
      std::vector<double>   GetDoubleVector(std::string str);
      std::vector<long int> GetIntVector(std::string str);

    private:
      bool     LoadParamSet(xmlDocPtr&, std::string cfg);
      void     ParseParamSet(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamDir(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamPos(std::string);
      void     ParseRotSeries(xmlDocPtr&, xmlNodePtr&);
      void     ParseWindowSeries(xmlDocPtr&, xmlNodePtr&);
      void     ParseEnuMax(std::string);
      void     ParseMaxWgtFail(std::string);
      TVector3 AnglesToAxis(double theta, double phi, std::string units = "deg");
      TVector3 ParseTV3(const std::string& );

      int         fVerbose;  ///< how noisy to be when parsing XML
      // who to apply these changes to
      GDk2NuFlux*  fGDk2NuFlux;

      // derived offsets/rotations
      TVector3    fBeamPosXML;
      TRotation   fBeamRotXML;
      TVector3    fFluxWindowPtXML[3];
    };
  }
}

//____________________________________________________________________________
namespace {
  std::string formatLorentzRotation(const TLorentzRotation& lr) {
    std::stringstream ss;
    std::ios_base::fmtflags original_flags = ss.flags();
    std::streamsize original_precision = ss.precision();

    ss << std::fixed << std::setprecision(6);
    ss << "  [ " << std::setw(9) << lr.XX()
       << ", " << std::setw(9) << lr.XY()
       << ", " << std::setw(9) << lr.XZ()
       << ", " << std::setw(9) << lr.XT() << " ]\n";
    ss << "  [ " << std::setw(9) << lr.YX()
       << ", " << std::setw(9) << lr.YY()
       << ", " << std::setw(9) << lr.YZ()
       << ", " << std::setw(9) << lr.YT() << " ]\n";
    ss << "  [ " << std::setw(9) << lr.ZX()
       << ", " << std::setw(9) << lr.ZY()
       << ", " << std::setw(9) << lr.ZZ()
       << ", " << std::setw(9) << lr.ZT() << " ]\n";
    ss << "  [ " << std::setw(9) << lr.TX()
       << ", " << std::setw(9) << lr.TY()
       << ", " << std::setw(9) << lr.TZ()
       << ", " << std::setw(9) << lr.TT() << " ]\n";

    ss.flags(original_flags);
    ss.precision(original_precision);

    return ss.str();
  }
}


//____________________________________________________________________________
GDk2NuFlux::GDk2NuFlux()
  : GFluxExposureI(genie::flux::kPOTs)
{
  this->Initialize();
}
//___________________________________________________________________________
GDk2NuFlux::~GDk2NuFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
double GDk2NuFlux::GetTotalExposure() const
{
  // complete the GFluxExposureI interface
  return UsedPOTs();
}
//___________________________________________________________________________
int GDk2NuFlux::PdgCode(void)
{
  return fCurNuChoice->pdgNu;
}
//___________________________________________________________________________
const TLorentzVector& GDk2NuFlux::Momentum(void)
{
  return fCurNuChoice->p4NuUser;
}
//___________________________________________________________________________
const TLorentzVector& GDk2NuFlux::Position(void)
{
  return fCurNuChoice->x4NuUser;
}
//___________________________________________________________________________
bool GDk2NuFlux::GenerateNext(void)
{
// Get next (unweighted) flux ntuple entry on the specified detector location
//
  RandomGen* rnd = RandomGen::Instance();
  while ( true ) {
     // Check for end of flux ntuple
     bool end = this->End();
     if ( end ) {
       LOG("Flux", pNOTICE) << "GenerateNext signaled End() ";
       return false;
     }

     // Get next weighted flux ntuple entry
     bool nextok = this->GenerateNext_weighted();
     if ( fGenWeighted ) {
       if ( nextok ) this->FinishCurrent();
       return nextok;
     }
     if ( ! nextok ) continue;

     /* RWH - debug purposes
     if ( fNCycles == 0 ) {
       LOG("Flux", pNOTICE)
          << "Got flux entry: " << fIEntry
          << " - Cycle: "<< fICycle << "/ infinite";
     } else {
       LOG("Flux", pNOTICE)
          << "Got flux entry: "<< fIEntry
          << " - Cycle: "<< fICycle << "/"<< fNCycles;
     }
     */

     // Get fractional weight & decide whether to accept curr flux neutrino
     double f = this->Weight() / fMaxWeight;
     //LOG("Flux", pNOTICE)
     //   << "Curr flux neutrino fractional weight = " << f;
     bool finished = false;
     if (f > 1.) {
       // complete the ray (full entry, times, move to z0) for the printout
       this->FinishCurrent();
       finished = true;
       // Oh,dear me ... what the weight max _should_ have been
       double bumped = this->Weight() * fMaxWgtFudge;
       if ( bumped > fMaxWeightMax ) fMaxWeightMax = bumped;
       ++fMaxWgtExceeded;

       // what to do ...
       if ( fMaxWgtFailModel <= 0 ) {
         fMaxWeight = bumped;
         LOG("Flux", pERROR)
           << "** Fractional weight = " << f
           << " > 1 !! Bump fMaxWeight estimate to " << fMaxWeight << "\n"
           << fCurDk2Nu->AsString() << "\n" << fCurNuChoice->AsString();
         std::cout << std::flush;
       } else {
         LOG("Flux", pERROR)
           << "** Fractional weight = " << f
           << " > 1 !! Leave fMaxWeight frozen at " << fMaxWeight << " "
           << " not bumped to " << bumped << ", max so far: "
           << fMaxWeightMax << "\n"
           << fCurDk2Nu->AsString() << "\n" << fCurNuChoice->AsString();
         std::cout << std::flush;
         if ( fMaxWgtFailModel >= 2 ) {
           PrintConfig();
           std::cout << std::flush;
           LOG("Flux", pFATAL)
             << "User requested not to continue";
           abort();
         }
       }
     }
     double r = rnd->RndFlux().Rndm(); // (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = ( r < f );
     if ( accept ) {

       // read the rest of the entry, set times, move to upstream z
       // (already done if the ray bumped fMaxWeight)
       if ( ! finished ) this->FinishCurrent();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("Flux", pNOTICE)
         << "Generated beam neutrino: "
         << "\n pdg-code: " << fCurNuChoice->pdgNu
         << "\n p4: " << utils::print::P4AsShortString(&(fCurNuChoice->p4NuBeam))
         << "\n x4: " << utils::print::X4AsString(&(fCurNuChoice->x4NuBeam))
         << "\n p4: " << utils::print::P4AsShortString(&(fCurNuChoice->p4NuUser))
         << "\n x4: " << utils::print::X4AsString(&(fCurNuChoice->x4NuUser));
#endif

       fWeight = 1.;
       return true;
     }

     //LOG("Flux", pNOTICE)
     //  << "** Rejecting current flux neutrino based on the flux weight only";
  }
  return false;
}
//___________________________________________________________________________
bool GDk2NuFlux::GenerateNext_weighted(void)
{
// Get next (weighted) flux ntuple entry on the specified detector location
//

  // Check whether a flux ntuple has been loaded
  if ( ! fNuFluxTree ) {
     LOG("Flux", pERROR)
          << "The flux driver has not been properly configured";
     return false;
  }

  // Reuse an entry?
  //std::cout << " ***** iuse " << fIUse << " nuse " << fNUse
  //          << " ientry " << fIEntry << " nentry " << fNEntries
  //          << " icycle " << fICycle << " ncycle " << fNCycles << std::endl;
  if ( fIUse < fNUse && fIEntry >= 0 ) {
    // Reuse this entry
    fIUse++;
  } else {
    // Reset previously generated neutrino code / 4-p / 4-x
    this->ResetCurrent();
    // Move on, read next flux ntuple entry
    fIEntry++;
    if ( fIEntry >= fNEntries ) {
      // Ran out of entries @ the current cycle of this flux file
      // Check whether more (or infinite) number of cycles is requested
      if (fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
      } else {
        LOG("Flux", pWARN)
          << "No more entries in input flux neutrino ntuple, cycle "
          << fICycle << " of " << fNCycles;
        fEnd = true;
        //assert(0);
        return false;
      }
    }

    // only bsim::Decay if partial reading is possible; see FinishCurrent()
    this->ReadCurrentEntry();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux",pDEBUG)
    << "got " << fNNeutrinos << " new fIEntry " << fIEntry
    << ( fIEntryFull == fIEntry ? " (full)" : " (decay only)" );
#endif

    fIUse = 1;

    // here we might want to do flavor oscillations or simple mappings
    fCurNuChoice->pdgNu  = fCurDk2Nu->decay.ntype;
    fCurNuChoice->impWgt = fCurDk2Nu->decay.nimpwt;
    fCurNuChoice->impWgtSum2 = fCurDk2Nu->decay.sumnimpwt2;
  }

  // update the # POTs & number of neutrinos
  // Do this HERE (before rejecting flavors that users might be weeding out)
  // in order to keep the POT accounting correct.  This allows one to get
  // the right normalization for generating only events from the intrinsic
  // nu_e entries.
  fAccumPOTs += fEffPOTsPerNu / fMaxWeight;
  fNNeutrinos++;

  // Check neutrino pdg against declared list of neutrino species declared
  // by the current instance of the neutrino flux driver.
  // No undeclared neutrino species will be accepted at this point as GENIE
  // has already been configured to handle the specified list.
  // Make sure that the appropriate list of flux neutrino species was set at
  // initialization via GDk2NuFlux::SetFluxParticles(const PDGCodeList &)
  // The result for the last flavour checked is cached (consecutive rays are
  // usually the same flavour); SetFluxParticles() resets the cache.

  if ( fCurNuChoice->pdgNu != fLastPdg ) {
    fLastPdg   = fCurNuChoice->pdgNu;
    fLastPdgOK = fPdgCList->ExistsInPDGCodeList(fLastPdg);
  }
  if ( ! fLastPdgOK ) {
     /// user might modify list via SetFluxParticles() in order to reject certain
     /// flavors, even if they're found in the file.  So don't make a big fuss.
     /// Spit out a single message and then stop reporting that flavor as problematic.
     int badpdg = fCurNuChoice->pdgNu;
     if ( ! fPdgCListRej->ExistsInPDGCodeList(badpdg) ) {
       fPdgCListRej->push_back(badpdg);
       LOG("Flux", pWARN)
         << "Encountered neutrino specie (" << badpdg
         << " that wasn't in SetFluxParticles() list, "
         << "\nDeclared list of neutrino species: " << *fPdgCList;
     }
     return false;
  }

  // Update the curr neutrino weight and energy

  // Check current neutrino energy against the maximum flux neutrino energy
  // declared by the current instance of the neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this
  // point as that maximum energy has already been used for normalizing the
  // interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at
  // initialization via GDk2NuFlux::SetMaxEnergy(double Ev)

  // Pick a point on the flux window and weight the ray, in beam coordinates.
  // Plain doubles throughout: each TVector3/TLorentzVector temporary costs a
  // TObject construction and destruction, which dominated this part of the
  // loop.  The products and sums keep the grouping of the vector expressions
  // they replace (base + (r1*dir1 + r2*dir2), TVector3::Unit, TVector3::Dot).
  RandomGen * rnd = RandomGen::Instance();
  // the two random numbers in a fixed order: in the former vector expression
  // r1*dir1 + r2*dir2 the order of the two Rndm() calls was up to the compiler
  const double r1 = rnd->RndFlux().Rndm();
  const double r2 = rnd->RndFlux().Rndm();

  const double wx = fFluxWindowBase.X() +
    ( r1*fFluxWindowDir1.X() + r2*fFluxWindowDir2.X() );
  const double wy = fFluxWindowBase.Y() +
    ( r1*fFluxWindowDir1.Y() + r2*fFluxWindowDir2.Y() );
  const double wz = fFluxWindowBase.Z() +
    ( r1*fFluxWindowDir1.Z() + r2*fFluxWindowDir2.Z() );
  const double wt = fFluxWindowBase.T() +
    ( r1*fFluxWindowDir1.T() + r2*fFluxWindowDir2.T() );
  fCurNuChoice->x4NuBeam.SetXYZT(wx, wy, wz, wt);

  double Ev = 0;
  double& wgt_xy = fCurNuChoice->xyWgt;
  fgWinPtBeam.SetXYZ(wx, wy, wz);   // reused member; calcEnuWgt takes a TVector3
  bsim::calcEnuWgt(fCurDk2Nu->decay, fgWinPtBeam, Ev, wgt_xy);

  if (Ev > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy"
          << "\nEv = " << Ev << "(> Ev{max} = " << fMaxEv << ")";
  }

  // Set the current flux neutrino 4-momentum, in *beam* coordinates, along
  // the direction from the decay vertex to the window point
  // (normalised as TVector3::Unit() does)
  const bsim::Decay& dk = fCurDk2Nu->decay;
  fgX4dkvtx.SetXYZT(dk.vx, dk.vy, dk.vz, 0.);
  double ux = wx - dk.vx;
  double uy = wy - dk.vy;
  double uz = wz - dk.vz;
  const double umag2 = ux*ux + uy*uy + uz*uz;
  const double uinv  = ( umag2 > 0 ) ? 1.0/TMath::Sqrt(umag2) : 1.0;
  ux *= uinv;  uy *= uinv;  uz *= uinv;
  fCurNuChoice->p4NuBeam.SetPxPyPzE(Ev*ux, Ev*uy, Ev*uz, Ev);

  fWeight = fCurNuChoice->impWgt * fCurNuChoice->xyWgt;  // full weight
  if ( fApplyTiltWeight ) {
    // additional weight due to window tilt relative to flux direction
    double nx = fFluxWindowNormal.X();
    double ny = fFluxWindowNormal.Y();
    double nz = fFluxWindowNormal.Z();
    if ( fIsSphere ) {
      // as SphereNormal(): radial direction at the point
      nx = wx - fFluxSphereCenterBeam.X();
      ny = wy - fFluxSphereCenterBeam.Y();
      nz = wz - fFluxSphereCenterBeam.Z();
      const double nmag2 = nx*nx + ny*ny + nz*nz;
      const double ninv  = ( nmag2 > 0 ) ? 1.0/TMath::Sqrt(nmag2) : 1.0;
      nx *= ninv;  ny *= ninv;  nz *= ninv;
    }
    const double tiltwgt = ux*nx + uy*ny + uz*nz;
    fWeight *= TMath::Abs( tiltwgt );
  }

  // the user-coordinate 4-momentum and 4-position, the ray times (which need
  // the ancestor list) and the move to the upstream z are done in
  // FinishCurrent(), once the ray has been accepted; only beam-coordinate
  // quantities are needed for the weight

  if ( Ev > fMaxEv ) {
    LOG("Flux", pFATAL)
      << "Generated neutrino had E_nu = " << Ev << " > " << fMaxEv
      << " maximum ";
    assert(0);
  }

  // update sum of weights
  fSumWeight += fWeight;

  return true;
}
//___________________________________________________________________________
void GDk2NuFlux::FinishCurrent(void)
{
  // Complete the current ray once it has been accepted, or before printing
  // it when it bumps fMaxWeight (or for every ray when generating weighted
  // rays); call at most once per ray: make sure the whole dk2nu entry has been
  // read, not only bsim::Decay, transform the ray to user coordinates, set
  // the ray times (which need the ancestor list) and, if requested, move the
  // ray origin upstream.

  this->ReadFullEntry();

  // the 4-momentum and 4-position in user coordinates
  Beam2UserP4(fCurNuChoice->p4NuBeam,fCurNuChoice->p4NuUser);
  Beam2UserPos(fCurNuChoice->x4NuBeam,fCurNuChoice->x4NuUser);

  // set the time component of the Lorentz vectors
  double dist_dk2start = GetDecayDist();  // dk point to window, user units
  size_t inu = fCurDk2Nu->indxnu();
  // calculate t_dk in user units
  double t_dk = 0;  // time from proton-on-target to where nu generated
  //if ( ! fCurDk2Nu->overflow() && ! fCurDk2Nu->ancestor.empty() ) {
  // overflow bit should signal that the ancestor list was truncated
  // but seems to be set unreliably
  // Assume the list is never truncated ...
  // if it is we'll just have a slightly wrong time (but better than zero)
  if ( ! fCurDk2Nu->ancestor.empty() ) {
    // ancestor times in beam time units (ns); t_dk in user units
    t_dk = fCurDk2Nu->ancestor[inu].startt * fTimeScaleB2U;
  }
  // this ... makes assumptions [blech]
  const double c_mbys  = 299792458;
  //const double c_user  = c_mbys * fLengthScale;
  double tstart = t_dk + ( dist_dk2start / c_mbys ); // user units

  fCurNuChoice->x4NuBeam.SetT(tstart * fTimeScaleU2B);
  fCurNuChoice->x4NuUser.SetT(tstart);
  if ( t_dk == 0 ) {
    // probably wasn't set
    static int nmsg = 5;
    if ( nmsg > 0 ) {
      --nmsg;
      LOG("Flux", pWARN)
        << "Setting time at flux window, \n"
        << "noticed that t_dk from ancestor list was 0, "
        << "this probably means that the calculated time is wrong";
    }
  }

  // if desired, move to user specified user coord z
  if ( TMath::Abs(fZ0) < 1.0e30 ) this->MoveToZ0(fZ0);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
    << "Generated neutrino: " << fIEntry << " " << fCurDk2Nu->potnum
    << "\n pdg-code: " << fCurNuChoice->pdgNu
    << "\n p4 beam: " << utils::print::P4AsShortString(&fCurNuChoice->p4NuBeam)
    << "\n x4 beam: " << utils::print::X4AsString(&fCurNuChoice->x4NuBeam)
    << "\n p4 user: " << utils::print::P4AsShortString(&(fCurNuChoice->p4NuUser))
    << "\n x4 user: " << utils::print::X4AsString(&(fCurNuChoice->x4NuUser));
#endif
}
//___________________________________________________________________________
void GDk2NuFlux::ReadCurrentEntry(void)
{
  // Read entry fIEntry of the flux chain.  With partial reading enabled and
  // a split dk2nu branch, only the "decay" sub-branch (everything needed to
  // weight the ray) is read; ReadFullEntry() reads the rest if the ray is
  // accepted.  Otherwise, or for a file whose dk2nu branch is not split,
  // the whole entry is read.  The check is made for each file in the chain,
  // so split and unsplit files may be mixed.

  fIEntryFull = -1;
  if ( fPartialRead ) {
    Long64_t local = fNuFluxTree->LoadTree(fIEntry);
    int itree = fNuFluxTree->GetTreeNumber();
    if ( itree != fDecayTreeNum ) {
      fDecayTreeNum = itree;
      fDecayBranch  = this->FindDecayBranch();
      if ( ! fDecayBranch ) {
        LOG("Flux", pNOTICE)
          << "dk2nu branch in file "
          << ( fNuFluxTree->GetFile() ? fNuFluxTree->GetFile()->GetName()
                                      : "(unknown)" )
          << " is not split; reading whole entries";
      }
    }
    if ( local >= 0 && fDecayBranch ) {
      fDecayBranch->GetEntry(local);
      return;
    }
  }
  fNuFluxTree->GetEntry(fIEntry);
  fIEntryFull = fIEntry;
}
//___________________________________________________________________________
void GDk2NuFlux::ReadFullEntry(void)
{
  // make sure all of entry fIEntry (not just bsim::Decay) is in fCurDk2Nu
  if ( fIEntry < 0 || fIEntryFull == fIEntry ) return;
  fNuFluxTree->GetEntry(fIEntry);
  fIEntryFull = fIEntry;
}
//___________________________________________________________________________
TBranch* GDk2NuFlux::FindDecayBranch(void)
{
  // the sub-branch holding bsim::Decay in the chain's current tree, if the
  // dk2nu branch was written split (split level >= 1); it is named "decay",
  // or "dk2nu.decay" if the top-level branch was booked as "dk2nu."
  TTree* tree = fNuFluxTree->GetTree();
  if ( ! tree ) return 0;
  TBranch* top = tree->GetBranch("dk2nu");
  if ( ! top || top->GetListOfBranches()->GetEntries() == 0 ) return 0;
  const char* names[] = { "decay", "dk2nu.decay" };
  for (size_t i = 0; i < 2; ++i) {
    TBranch* b = (TBranch*)top->GetListOfBranches()->FindObject(names[i]);
    if ( b ) return b;
  }
  return 0;
}
//___________________________________________________________________________
double GDk2NuFlux::GetDecayDist() const
{
  // return distance (user units) between dk point and start (window) position
  // these are in beam units
  //                position on window              dk2nu nu start
  TVector3 x3diff = fCurNuChoice->x4NuBeam.Vect() - fgX4dkvtx.Vect();
  return x3diff.Mag() * fLengthScaleB2U;
}
//___________________________________________________________________________
void GDk2NuFlux::MoveToZ0(double z0usr)
{
  // move ray origin to specified user z0
  // move beam coord entry correspondingly

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  static bool first = true;
  if (first) {
    first = false;
    LOG("Flux", pDEBUG)
      << "\n"
      << "TLorentzRotation fBeamRot" << "\n"
      << formatLorentzRotation(fBeamRot) << "\n"
      << "TLorentzRotation fBeamRotInv" << "\n"
      << formatLorentzRotation(fBeamRotInv);
  }

  LOG("Flux", pDEBUG)
    << "MoveToZ0 (z0usr=" << z0usr << ") before:"
    << "\n p4 user: " << utils::print::P4AsShortString(&(fCurNuChoice->p4NuUser))
    << "\n x4 user: " << utils::print::X4AsString(&(fCurNuChoice->x4NuUser))
    << "\n x4 beam: " << utils::print::X4AsString(&(fCurNuChoice->x4NuBeam));
#endif

  double pzusr    = fCurNuChoice->p4NuUser.Pz();
  if ( TMath::Abs(pzusr) < 1.0e-30 ) {
    // neutrino is moving almost entirely in x-y plane
    LOG("Flux", pWARN)
      << "MoveToZ0(" << z0usr << ") not possible due to pz_usr (" << pzusr << ")";
    return;
  }

  double dz = z0usr - fCurNuChoice->x4NuUser.Z();
  double scale = dz / pzusr;
  TVector3 dx3( scale*fCurNuChoice->p4NuUser.Vect() );

  // Find the corresponding time shift
  double v = fCurNuChoice->p4NuUser.Beta() * constants::kLightSpeed
    / ( units::meter / units::second );
  double dt = fCurNuChoice->p4NuUser.P() * dz / ( pzusr * v );

  // appl thes shifts to update the neutrino 4-position
  TLorentzVector dx4(dx3,dt);
  fCurNuChoice->x4NuUser += dx4;
  // rotations this a different dx4, so just the the coordinate transform
  User2BeamPos(fCurNuChoice->x4NuUser,fCurNuChoice->x4NuBeam);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pDEBUG)
    << "MoveToZ0 (" << z0usr << ") after:"
    << "\n x4 user: " << utils::print::X4AsString(&(fCurNuChoice->x4NuUser))
    << "\n x4 beam: " << utils::print::X4AsString(&(fCurNuChoice->x4NuBeam));
#endif
}

//___________________________________________________________________________
void GDk2NuFlux::CalcEffPOTsPerNu()
{
  // do this if flux window changes or # of files changes

  if (!fNuFluxTree) return;  // not yet fully configured

  // effpots = mc_pots * (wgtfunction-area) / window-area / wgt-max-est
  //   wgtfunction-area = pi * radius-det-element^2 = pi * (100.cm)^2

  // this should match what is used in the CalcEnuWgt()
  const double kRDET = 100.0;   // set to flux per 100 cm radius
  const double kRDET2 = kRDET * kRDET;
  double flux_area = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Mag();
  LOG("Flux",pNOTICE) << "in CalcEffPOTsPerNu, area = " << flux_area;

  if ( flux_area < 1.0e-30 ) {
    LOG("Flux", pWARN)
          << "CalcEffPOTsPerNu called with flux window area effectively zero";
    flux_area = 1;
  }
  double area_ratio = TMath::Pi() * kRDET2 / flux_area;
  fEffPOTsPerNu = area_ratio * ( (double)fFilePOTs / (double)fNEntries );
}

//___________________________________________________________________________
void GDk2NuFlux::LoadDkMeta(void)
{
  // load the matching bsim::dkmeta entry that goes w/ bsim::dk2nu entry

  if ( fJobToMetaIndex.empty() ) {
    int nmeta = fNuMetaTree->GetEntries();
    for (int imeta =0; imeta < nmeta; ++imeta ) {
      fNuMetaTree->GetEntry(imeta);
      int mjob = fCurDkMeta->job;
      // there shouldn't already be an entry in the map
      // complain if there is
      std::map<int,int>::const_iterator mitr = fJobToMetaIndex.find(mjob);
      if ( mitr == fJobToMetaIndex.end() ) {
        fJobToMetaIndex[mjob] = imeta;  // make an entry
      } else {
        LOG("Flux", pERROR) << "LoadDkMeta already had an entry for job "
                            << mjob << " at " << mitr->second
                            << " which conflicts with new entry at "
                            << imeta;
      }
    }
  }

  int job = fCurDk2Nu->job;
  if ( fCurDkMeta->job == job ) return;  // already loaded

  std::map<int,int>::const_iterator lookat = fJobToMetaIndex.find(job);
  if ( lookat != fJobToMetaIndex.end() ) {
    int indx = lookat->second;
    fNuMetaTree->GetEntry(indx);
    if ( fCurDkMeta->job != job ) {
      LOG("Flux", pFATAL) << "Failed to get right metadata " << job
                          << " => " << fCurDkMeta->job
                          << " indx " << indx;
      assert(0);
    }
  } else {
    // wasn't already indexed
    // (lookat == end() here, so it must not be dereferenced)
    LOG("Flux", pFATAL) << "Failed to find metadata for job " << job
                        << " among " << fJobToMetaIndex.size()
                        << " indexed dkmeta entries";
    assert(0);
  }

}

//___________________________________________________________________________
double GDk2NuFlux::UsedPOTs(void) const
{
// Compute current number of flux POTs

  if (!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;
  }
  return fAccumPOTs;
}

//___________________________________________________________________________
double GDk2NuFlux::POT_curr(void) {
  // RWH: Not sure what POT_curr is supposed to represent I'll guess for
  // now that that it means what I mean by UsedPOTs().
  return UsedPOTs();
}

//___________________________________________________________________________
void GDk2NuFlux::LoadBeamSimData(const std::vector<std::string>& patterns,
                                 const std::string&         config )
{
// Loads in a beam simulation root file into the GDk2NuFlux driver.

  bool found_cfg = this->LoadConfig(config);
  if ( ! found_cfg ) {
    LOG("Flux", pFATAL)
      << "LoadBeamSimData could not find XML config \"" << config << "\"\n";
    exit(1);
  }

  fNuFluxFilePatterns = patterns;
  std::vector<int> nfiles_from_pattern;

  // create a (sorted) set of file names
  // this also helps ensure that the same file isn't listed multiple times
  std::set<std::string> fnames;

  LOG("Flux", pINFO) << "LoadBeamSimData was passed " << patterns.size()
                     << " patterns";

  for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
    std::string pattern = patterns[ipatt];
    nfiles_from_pattern.push_back(0);

    LOG("Flux", pNOTICE)
      << "Loading dk2nu flux tree from ROOT file pattern ["
      << std::setw(3) << ipatt << "] \"" << pattern << "\"";

    // !WILDCARD only works for file name ... NOT directory
    std::string dirname = gSystem->UnixPathName(gSystem->WorkingDirectory());
    size_t slashpos = pattern.find_last_of("/");
    size_t fbegin;
    if ( slashpos != std::string::npos ) {
      dirname = pattern.substr(0,slashpos);
      LOG("Flux", pINFO) << "Look for flux using directory " << dirname;
      fbegin = slashpos + 1;
    } else { fbegin = 0; }
    std::string basename = pattern.substr(fbegin,pattern.size()-fbegin);

    // A plain file name (no wildcard) is used as given, without listing
    // the directory: listing a remote (e.g. root://) directory needs the
    // server to allow it, and can be slow for large directories.
    // TFile::Open below reports a file that can't be opened.
    bool has_wildcard = ( basename.find_first_of("*?[") != std::string::npos );
    if ( ! has_wildcard ) {
      fnames.insert(pattern);
      nfiles_from_pattern[ipatt]++;
      continue;
    }

    char* expdir = gSystem->ExpandPathName(dirname.c_str());
    void* dirp = gSystem->OpenDirectory(expdir);
    delete [] expdir;
    if ( dirp ) {
      TRegexp re(basename.c_str(),kTRUE);
      const char* onefile;
      while ( ( onefile = gSystem->GetDirEntry(dirp) ) ) {
        TString afile = onefile;
        if ( afile=="." || afile==".." ) continue;
        if ( basename!=afile && afile.Index(re) == kNPOS ) continue;
        std::string fullname = dirname + "/" + afile.Data();
        fnames.insert(fullname);
        nfiles_from_pattern[ipatt]++;
      }
      gSystem->FreeDirectory(dirp);
      if ( nfiles_from_pattern[ipatt] == 0 ) {
        LOG("Flux", pWARN)
          << "No files in directory \"" << dirname
          << "\" match \"" << basename << "\"";
      }
    } else {
      LOG("Flux", pERROR)
        << "Could not list directory \"" << dirname
        << "\" to expand the wildcard in \"" << basename << "\""
        << " (for a remote directory the server must allow listing);"
        << " give explicit file names instead";
    } // legal directory
  } // loop over patterns

  // std::cout << "RWH start adding files" << std::endl << std::flush;

  size_t indx = 0;
  std::set<std::string>::const_iterator sitr = fnames.begin();
  for ( ; sitr != fnames.end(); ++sitr, ++indx ) {
    std::string filename = *sitr;
    // std::cout << "RWH  [" << std::setw(3) << indx << "]  \""
    //           << filename << "\"" << std::endl << std::flush;
    // gSystem->AccessPathName() only works for local files, so we can't use
    // it if we want to stream via xrootd; rely on TFile::Open instead
    TFile* tf = TFile::Open(filename.c_str(),"READ");
    if ( ! tf || tf->IsZombie() ) {
      LOG("Flux", pERROR) << "Could not open flux file \"" << filename
                          << "\"; it was not added to the file list";
      delete tf;
      continue;
    }
    {
      TTree* ftree = (TTree*)tf->Get(fTreeNames[0].c_str());
      TTree* mtree = (TTree*)tf->Get(fTreeNames[1].c_str());
      if ( ftree && mtree ) {
        // found both trees
        this->AddFile(ftree,mtree,filename);
        // std::cout << "RWH  AddFile(" << ftree << ","
        //    << mtree << "," << filename << ")" << std::endl << std::flush;
      } else {
        // std::cout << "RWH  no AddFile" << std::endl << std::flush;
        LOG("Flux", pNOTICE) << "File " << filename << " lacked a tree: "
                             << " \"" << fTreeNames[0] << "\" " << ftree
                             << " \"" << fTreeNames[1] << "\" " << mtree
                             << " and was not added to the file list";
      }
      tf->Close();
      delete tf;
    } // loop over tree type
  } // loop over sorted file names

  if ( ! fNuFluxTree || ! fNuMetaTree ) {
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    LOG("Flux", pERROR)
      << "Missing dk2nu (" << fNuFluxTree << ") or dkmeta ("
      << fNuMetaTree << ") tree";
    LOG("Flux", pERROR)
      << "please check that input files are valid";
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
  }

  // this will open all files and read header!!
  fNEntries = fNuFluxTree->GetEntries();

  if ( fNEntries == 0 ) {
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    LOG("Flux", pERROR)
      << "Loaded flux tree contains " <<  fNEntries << " entries";
    LOG("Flux", pERROR)
      << "Was passed the file patterns: ";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      std::string pattern = patterns[ipatt];
      LOG("Flux", pERROR)
        << "  [" << std::setw(3) << ipatt <<"] " << pattern;
    }
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
  } else {
    LOG("Flux", pNOTICE)
      << "Loaded flux tree contains " <<  fNEntries << " entries"
      << " from " << fnames.size() << " unique files";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      std::string pattern = patterns[ipatt];
      LOG("Flux", pINFO)
        << " pattern: " << pattern << " yielded "
        << nfiles_from_pattern[ipatt] << " files";
    }
  }

  TBranch* bdk2nu = fNuFluxTree->GetBranch("dk2nu");
  if (bdk2nu) bdk2nu->SetAutoDelete(false);
  TBranch* bdkmeta = fNuMetaTree->GetBranch("dkmeta");
  if (bdkmeta) bdkmeta->SetAutoDelete(false);

  // chain (re)built: forget any cached "decay" sub-branch / full entry
  fDecayTreeNum = -1;
  fDecayBranch  =  0;
  fIEntryFull   = -1;

  // we have a file we can work with
  if (!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "LoadBeamSimData left detector location unset";
  }
  if (fMaxWeight<=0) {
     LOG("Flux", pINFO)
       << "Run ScanForMaxWeight() as part of LoadBeamSimData";
     this->ScanForMaxWeight();
  }

  if ( fPartialRead ) {
    // The TTreeCache learns which branches to prefetch from the first
    // entries read.  With partial reading those come from the max-weight
    // scan, which reads only "decay", so the cache would fetch only the
    // decay baskets and every other basket needed by an accepted entry
    // would be a separate read (one round trip each over xrootd/dCache).
    // Accepted entries need nearly every basket anyway, so register all
    // branches now that the scan is done: generation then makes the same
    // reads as full reading, and the scan read only "decay".
    // (A chain sizes its cache only once a tree is loaded; with no cache,
    // e.g. SetCacheSize(0), there is nothing to register.)
    if ( fNuFluxTree->GetTreeNumber() < 0 ) fNuFluxTree->LoadTree(0);
    if ( fNuFluxTree->GetCacheSize() > 0 &&
         fNuFluxTree->AddBranchToCache("*", kTRUE) == 0 ) {
      fNuFluxTree->StopCacheLearningPhase();
    } else {
      LOG("Flux", pINFO) << "no TTreeCache in use for the flux chain";
    }
  }

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle =  0;
  // pick a starting entry index [0:fNEntries-1]
  // pretend we just used up the the previous one
  RandomGen* rnd = RandomGen::Instance();
  fIUse   =  9999999;
  fIEntry = rnd->RndFlux().Integer(fNEntries) - 1;

  // don't count things we used to estimate max weight
  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;

  LOG("Flux",pNOTICE) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();

}

//___________________________________________________________________________
void GDk2NuFlux::GetBranchInfo(std::vector<std::string>& branchNames,
                               std::vector<std::string>& branchClassNames,
                               std::vector<void**>&      branchObjPointers)
{
  // allow flux driver to report back current status and/or ntuple entry
  // info for possible recording in the output file by supplying
  // the class name, and a pointer to the object that will be filled
  // as well as a suggested name for the branch.

  branchNames.push_back("dk2nu");
  branchClassNames.push_back("bsim::Dk2Nu");
  branchObjPointers.push_back((void**)&fCurDk2Nu);

  branchNames.push_back("nuchoice");
  branchClassNames.push_back("bsim::NuChoice");
  branchObjPointers.push_back((void**)&fCurNuChoice);

}
TTree* GDk2NuFlux::GetMetaDataTree() { return fNuMetaTree; }

//___________________________________________________________________________
void GDk2NuFlux::SetFluxParticles(const PDGCodeList & particles)
{
  // set the list of flavours to generate; forget the cached result of
  // the last flavour check, which was made against the old list
  GFluxFileConfigI::SetFluxParticles(particles);
  fLastPdg   = 0;
  fLastPdgOK = false;
}


//___________________________________________________________________________
void GDk2NuFlux::ScanForMaxWeight(void)
{
  if (!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "Specify a detector location before scanning for max weight";
     return;
  }

  // scan for the maximum weight

  double wgtgenmx = 0, enumx = 0;
  TStopwatch t;
  t.Start();
  for (int itry=0; itry < fMaxWgtEntries; ++itry) {
    this->GenerateNext_weighted();
    double wgt = this->Weight();
    if ( wgt > wgtgenmx ) wgtgenmx = wgt;
    double enu = fCurNuChoice->p4NuBeam.Energy();
    if ( enu > enumx ) enumx = enu;
  }
  t.Stop();
  LOG("Flux", pNOTICE) << "Maximum flux weight for spin = "
                       << wgtgenmx << ", energy = " << enumx
                       << " (" << fMaxWgtEntries << " rays, "
                       << t.CpuTime() << " s CPU, "
                       << t.RealTime() << " s real)";

  if (wgtgenmx > fMaxWeight ) fMaxWeight = wgtgenmx;
  // apply a fudge factor to estimated weight
  fMaxWeight *= fMaxWgtFudge;
  if ( fMaxWgtFudge < 1 ) fMaxWgtFudge=1; // avoid run-away downward spiral

  fMaxWeightScan = fMaxWeight;
  // apply possible user supplied minimum (floor)
  if (fMinMaxWeight > fMaxWeight) {
    LOG("Flux", pNOTICE) << "Maximum flux weight set by user = "
                         << fMinMaxWeight << " exceeds estimate, use that";
    fMaxWeight = fMinMaxWeight;
  }

  fMaxWeightInit  = fMaxWeight; // what we started with before generation
  fMaxWeightMax   = fMaxWeight;
  fMaxWgtExceeded = 0;          // not yet exceeded

  // adjust max energy?
  if ( enumx*fMaxEFudge > fMaxEv ) {
    LOG("Flux", pNOTICE) << "Adjust max: was=" << fMaxEv
                         << " now " << enumx << "*" << fMaxEFudge
                         << " = " << enumx*fMaxEFudge;
    fMaxEv = enumx * fMaxEFudge;
  }

  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight
                       << ", energy = " << fMaxEv;

}
//___________________________________________________________________________
void GDk2NuFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO)
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GDk2NuFlux::SetEntryReuse(long int nuse)
{
// With nuse > 1 then the same entry in the file is used "nuse" times
// before moving on to the next entry in the ntuple

  fNUse    = TMath::Max(1L, nuse);
}

//___________________________________________________________________________
void GDk2NuFlux::SetFluxWindow(TVector3 p0, TVector3 p1, TVector3 p2)
                             // bool inDetCoord)  future extension
{
  // set flux window
  // NOTE: internally these are in "cm", but user might have set a preference
  fDetLocIsSet         = true;
  fIsSphere            = false;

  fFluxWindowPtUser[0] = p0;
  fFluxWindowPtUser[1] = p1;
  fFluxWindowPtUser[2] = p2;

  // convert from user to beam coord and from 3 points to base + 2 directions
  // apply units conversion
  TLorentzVector ptbm0, ptbm1, ptbm2;
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[0],0),ptbm0);
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[1],0),ptbm1);
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[2],0),ptbm2);

  fFluxWindowBase = ptbm0;
  fFluxWindowDir1 = ptbm1 - ptbm0;
  fFluxWindowDir2 = ptbm2 - ptbm0;

  fFluxWindowLen1   = fFluxWindowDir1.Mag();
  fFluxWindowLen2   = fFluxWindowDir2.Mag();
  fFluxWindowNormal = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Unit();

  double dot = fFluxWindowDir1.Dot(fFluxWindowDir2);
  if ( TMath::Abs(dot) > 1.0e-8 )
    LOG("Flux",pWARN) << "Dot product between window direction vectors was "
                      << dot << "; please check for orthoganality";

  LOG("Flux",pNOTICE) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();
}

//___________________________________________________________________________
void GDk2NuFlux::SetFluxSphere(TVector3 center, double radius,
                               bool inDetCoord)
{
  fIsSphere = true;
  fFluxSphereRadius = radius;
  if ( inDetCoord ) {
    fFluxSphereCenterUser = center;
    TLorentzVector v4beam;
    User2BeamPos(TLorentzVector(center,0),v4beam);
    fFluxSphereCenterBeam = v4beam.Vect();
  } else {
    fFluxSphereCenterBeam = center;
    TLorentzVector v4user;
    Beam2UserPos(TLorentzVector(center,0),v4user);
    fFluxSphereCenterUser = v4user.Vect();
  }

  // not yet supported
  LOG("Flux", pFATAL) << "use of spherical flux window not yet fully supported";
  assert(0);

}

//___________________________________________________________________________
void GDk2NuFlux::GetFluxSphere(TVector3& center, double& radius,
                               bool inDetCoord) const
{
  center = ( inDetCoord ) ? fFluxSphereCenterUser : fFluxSphereCenterBeam;
  radius = fFluxSphereRadius;
}

//___________________________________________________________________________
void GDk2NuFlux::GetFluxWindow(TVector3& p0, TVector3& p1, TVector3& p2) const
{
  // return flux window points
  p0 = fFluxWindowPtUser[0];
  p1 = fFluxWindowPtUser[1];
  p2 = fFluxWindowPtUser[2];

}
//___________________________________________________________________________
void GDk2NuFlux::SetBeamRotation(TRotation beamrot)
{
  // rotation is really only 3-d vector, but we'll be operating on LorentzV's
  // unless we do a boost or something this leaves [XYZ]T and T[XYZ] = 0, TT=1
  fBeamRot    = TLorentzRotation(beamrot);
  fBeamRotInv = fBeamRot.Inverse();
}

void GDk2NuFlux::SetBeamCenter(TVector3 beam0)
{
  // set coord transform between detector and beam
  // NOTE: internally these are in "cm", but user might have set a preference
  fBeamZero = TLorentzVector(beam0,0);  // no time shift
}

//___________________________________________________________________________
TRotation GDk2NuFlux::GetBeamRotation() const
{
  // rotation is really only 3-d vector, but we'll be operating on LorentzV's
  // give people back the original TRotation ... not pretty
  // ... it think this is right
  TRotation rot3;
  const TLorentzRotation& rot4 = fBeamRot;
  TVector3 newX(rot4.XX(),rot4.XY(),rot4.XZ());
  TVector3 newY(rot4.YX(),rot4.YY(),rot4.YZ());
  TVector3 newZ(rot4.ZX(),rot4.ZY(),rot4.ZZ());
  rot3.RotateAxes(newX,newY,newZ);
  return rot3.Inverse();
}
TVector3 GDk2NuFlux::GetBeamCenter() const
{
  TVector3 beam0 = fBeamZero.Vect();
  return beam0;
}

//___________________________________________________________________________
TVector3  GDk2NuFlux::SphereNormal() const
{
  // return the normal for this current location on the sphere
  TVector3 dir = (fCurNuChoice->x4NuBeam.Vect() - fFluxSphereCenterBeam);
  return dir.Unit();
}

//___________________________________________________________________________
//void GDk2NuFlux::SetCoordTransform(TVector3 beam0, TRotation beamrot)
//{
//  // set coord transform between detector and beam
//  // NOTE: internally these are in "cm", but user might have set a preference
//
//  beam0 *= (1./fLengthScaleB2U);
//  fDetectorZero = TLorentzVector(beam0,0);  // no time shift
//  fDetectorRot  = TLorentzRotation(beamrot);
//
//}
//___________________________________________________________________________
//void GDk2NuFlux::GetDetectorCoord(TLorentzVector& det0, TLorentzRotation& detrot) const
//{
//  // get coord transform between detector and beam
//  // NOTE: internally these are in "cm", but user might have set a preference
//
//  det0 = fDetectorZero;
//  det0 *= fLengthScaleB2U;
//  detrot = fDetectorRot;
//
//}
//___________________________________________________________________________

void GDk2NuFlux::Beam2UserPos(const TLorentzVector& beamxyz,
                                   TLorentzVector& usrxyz) const
{
  usrxyz = fLengthScaleB2U*(fBeamRot*beamxyz) + fBeamZero;
  // above assumed fBeamZero T=0
  // fBeamRot having no T components other than TT=1
  usrxyz.SetT(beamxyz.T()*fTimeScaleB2U);
}
void GDk2NuFlux::Beam2UserDir(const TLorentzVector& beamdir,
                                   TLorentzVector& usrdir) const
{
  usrdir = fLengthScaleB2U*(fBeamRot*beamdir);
}
void GDk2NuFlux::Beam2UserP4 (const TLorentzVector& beamp4,
                                   TLorentzVector& usrp4 ) const
{
  usrp4 = fBeamRot*beamp4;
}

void GDk2NuFlux::User2BeamPos(const TLorentzVector& usrxyz,
                                   TLorentzVector& beamxyz) const
{
  beamxyz = fLengthScaleU2B*(fBeamRotInv*(usrxyz-fBeamZero));
  // above assumed fBeamZero T=0
  // fBeamRotInv having no T components other than TT=1
  beamxyz.SetT(usrxyz.T()*fTimeScaleU2B);
}
void GDk2NuFlux::User2BeamDir(const TLorentzVector& usrdir,
                                   TLorentzVector& beamdir) const
{
  beamdir = fLengthScaleU2B*(fBeamRotInv*usrdir);
}
void GDk2NuFlux::User2BeamP4 (const TLorentzVector& usrp4,
                                   TLorentzVector& beamp4) const
{
  beamp4 = fBeamRotInv*usrp4;
}

//___________________________________________________________________________
void GDk2NuFlux::PrintCurrent(void)
{
  LOG("Flux", pNOTICE) << "CurrentEntry:\n"
                       << fCurDk2Nu->AsString() << "\n"
                       << fCurNuChoice->AsString();
}
//___________________________________________________________________________
void GDk2NuFlux::Clear(Option_t * opt)
{
  // Clear the driver state
  //
  LOG("Flux", pWARN) << "GSimpleNtpFlux::Clear(" << opt << ") called";
  // do it in all cases, but EVGDriver/GMCJDriver will pass "CycleHistory"

  fICycle     = 0;

  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;
}
//___________________________________________________________________________
void GDk2NuFlux::GenerateWeighted(bool gen_weighted)
{
  // Set whether to generate weighted rays
  //
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GDk2NuFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GDk2NuFlux driver";

  fMaxEv           =  0;
  fEnd             =  false;


  fTreeNames[0]    = "dk2nuTree";
  fTreeNames[1]    = "dkmetaTree";
  fNuFluxTree      =  0;
  fNuMetaTree      =  0;
  fCurDk2Nu        =  0;
  fCurDkMeta       =  0;
  fCurNuChoice     =  0;
  fNFiles          =  0;

  fNEntries        =  0;
  fIEntry          = -1;
  fIEntryFull      = -1;
  fPartialRead     = true;
  fLastPdg         =  0;
  fLastPdgOK       = false;
  fDecayTreeNum    = -1;
  fDecayBranch     =  0;
  fICycle          =  0;
  fNUse            =  1;
  fIUse            =  999999;

  fNuTot           = 0;
  fFilePOTs        = 0;

  fMaxWeight       = -1;
  fMinMaxWeight    =  0;
  fMaxWgtFudge     =  1.05;
  fMaxWgtEntries   = 2500000;
  fMaxEFudge       =  1.05;

  fMaxWgtExceeded  =  0;
  fMaxWgtFailModel =  0;  // default: 0=bump, 1=leave frozen, 2=abort

  fSumWeight       =  0;
  fNNeutrinos      =  0;
  fEffPOTsPerNu    =  0;
  fAccumPOTs       =  0;

  fGenWeighted     = false;
  fApplyTiltWeight = true;
  fIsSphere        = false;
  fDetLocIsSet     = false;
  // by default assume user length is meters, time is seconds
  SetLengthUnits(genie::utils::units::UnitFromString("m"));
  SetTimeUnits(genie::utils::units::UnitFromString("s"));



  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GDk2NuFlux::SetDefaults(void)
{
// - Set default neutrino species list (nue, nuebar, numu, numubar) and
//   maximum energy (120 GeV).
//   These defaults can be overwritten by user calls (at the driver init) to
//   GDk2NuFlux::SetMaxEnergy(double Ev) and
//   GDk2NuFlux::SetFluxParticles(const PDGCodeList & particles)
// - Set the default file normalization to 1E+21 POT
// - Set the default flux neutrino start z position at -5m (z=0 is the
//   detector centre).
// - Set number of cycles to 1

  LOG("Flux", pNOTICE) << "Setting default GDk2NuFlux driver options";

  PDGCodeList particles;
  particles.push_back(kPdgNuMu);
  particles.push_back(kPdgAntiNuMu);
  particles.push_back(kPdgNuE);
  particles.push_back(kPdgAntiNuE);

  this->SetFluxParticles (particles);
  this->SetMaxEnergy     (120./*GeV*/);  // was 200, but that would be wasteful
  this->SetUpstreamZ     (-3.4e38); // way upstream ==> use flux window
  this->SetNumOfCycles   (0);
  this->SetEntryReuse    (1);

  std::string xmlfile = "GNuMIFlux.xml";
  const char* altxml = gSystem->Getenv("GDK2NUFLUXXML");
  if ( altxml ) xmlfile = altxml;
  this->SetXMLFileBase(xmlfile);

}
//___________________________________________________________________________
void GDk2NuFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  if ( fCurDk2Nu )    fCurDk2Nu->clear();
  if ( fCurNuChoice ) fCurNuChoice->clear();
}
//___________________________________________________________________________
void GDk2NuFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if ( fPdgCList )    delete fPdgCList;
  if ( fPdgCListRej ) delete fPdgCListRej;
  if ( fCurNuChoice ) delete fCurNuChoice;

  LOG("Flux", pNOTICE)
    << " flux file cycles: " << fICycle << " of " << fNCycles
    << ", entry " << fIEntry << " use: " << fIUse << " of " << fNUse;
}

//___________________________________________________________________________
void GDk2NuFlux::AddFile(TTree* ftree, TTree* mtree, std::string fname)
{
  // Add a file to the chain

  ULong64_t nentries = ftree->GetEntries();

  // generally files will only have one meta data entry, but if they've
  // been combined (i.e. "hadd") there might be more than one.
  int nmeta = mtree->GetEntries();
  double potsum = 0;
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
  mtree->SetBranchAddress("dkmeta",&dkmeta);
  for (int imeta = 0; imeta < nmeta; ++imeta ) {
    mtree->GetEntry(imeta);
    double potentry = dkmeta->pots;
    potsum += potentry;
  }
  delete dkmeta;

  // don't need these anymore
  delete ftree;
  delete mtree;

  // make sure the chains are defined
  bool firstFile = ( fNuFluxTree == 0 );
  if ( firstFile ) {
    fNuFluxTree  = new TChain(fTreeNames[0].c_str());
    fNuMetaTree  = new TChain(fTreeNames[1].c_str());
    fCurDk2Nu    = new bsim::Dk2Nu;
    fCurDkMeta   = new bsim::DkMeta;
    fCurNuChoice = new bsim::NuChoice;
  }

  // add the file to the chains
  int stat0 = fNuFluxTree->AddFile(fname.c_str());
  int stat1 = fNuMetaTree->AddFile(fname.c_str());

  // attach the branch objects once the chains have a file; ROOT warns
  // ("no dataset schema available") if this is done on an empty chain
  if ( firstFile ) {
    fNuFluxTree->SetBranchAddress("dk2nu",&fCurDk2Nu);
    fNuMetaTree->SetBranchAddress("dkmeta",&fCurDkMeta);
  }

  LOG("Flux",pINFO)
    << "flux->AddFile() of " << nentries
    << " " << ((mtree)?"[+meta]":"[no-meta]")
    << " [status=" << stat0 << "," << stat1 << "]"
    << nentries << " (" << nmeta << ")"
    << " entries in file: " << fname;

  if ( stat0 != 1 || stat1 != 1 ) {
    SLOG("GDk2NuFlux", pFATAL) << "Add: \"" << fname << "\" failed";
  }

  fNuTot    += nentries;
  fFilePOTs += potsum;
  fNFiles++;

}

//___________________________________________________________________________
void GDk2NuFlux::SetLengthUnits(double user_units_dist)
{
  // Set the scale factor for lengths going from beam (cm) to user coordinates

  // GDk2NuFlux uses "cm" as the length unit consistently internally (this is
  // the length units used by both the g3 and g4 ntuples).  User interactions
  // setting the beam-to-detector coordinate transform, flux window, and the
  // returned position might need to be in other units.  Use:
  //     double scale = genie::utils::units::UnitFromString("cm");
  // ( #include "Utils/UnitUtils.h for declaration )
  // to get the correct scale factor to pass in.

  double rescale = fLengthUnits / user_units_dist;
  fLengthUnits = user_units_dist;
  double cm = genie::utils::units::UnitFromString("cm");
  fLengthScaleB2U = cm / user_units_dist;
  fLengthScaleU2B = user_units_dist / cm;

  // in case we're changing units without resetting transform/window
  // not recommended, but should work
  if (fCurNuChoice) fCurNuChoice->x4NuUser  *= rescale;
  fBeamZero               *= rescale;
  fFluxWindowPtUser[0]    *= rescale;
  fFluxWindowPtUser[1]    *= rescale;
  fFluxWindowPtUser[2]    *= rescale;
}

//___________________________________________________________________________
void GDk2NuFlux::SetTimeUnits(double user_units_time)
{
  // Set the scale factor for time going from beam (ns) to user coordinates

  // GDk2NuFlux uses "ns" as the length unit consistently internally (this is
  // the length units used by both the g3 and g4 ntuples).  User interactions
  // setting the beam-to-detector coordinate transform, flux window, and the
  // returned position might need to be in other units.  Use:
  //     double scale = genie::utils::units::UnitFromString("s");
  // ( #include "Utils/UnitUtils.h for declaration )
  // to get the correct scale factor to pass in.

  fTimeUnits = user_units_time;
  double ns = genie::utils::units::UnitFromString("ns");
  fTimeScaleB2U = ns / user_units_time;
  fTimeScaleU2B = user_units_time / ns;

}

//___________________________________________________________________________
double GDk2NuFlux::LengthUnits(void) const
{
  // Return the user length unit (as set by SetLengthUnits) relative to
  // metres, e.g. 1 for "m", 0.01 for "cm".  This is the same convention
  // as ROOTGeomAnalyzer::LengthUnits().
  return fLengthUnits / units::meter;
}
//___________________________________________________________________________
double GDk2NuFlux::TimeUnits(void) const
{
  // Return the user time unit (as set by SetTimeUnits) relative to
  // seconds, e.g. 1 for "s", 1.0e-9 for "ns".
  return fTimeUnits / units::second;
}

//___________________________________________________________________________

bool GDk2NuFlux::LoadConfig(std::string cfg)
{
  genie::flux::GDk2NuFluxXMLHelper helper(this);
  return helper.LoadConfig(cfg);
}

//___________________________________________________________________________

void GDk2NuFlux::PrintConfig()
{

  std::ostringstream s;
  PDGCodeList::const_iterator itr = fPdgCList->begin();
  for ( ; itr != fPdgCList->end(); ++itr) s << (*itr) << " ";
  s << "[rejected: ";
  itr = fPdgCListRej->begin();
  for ( ; itr != fPdgCListRej->end(); ++itr) s << (*itr) << " ";
  s << " ] ";

  std::ostringstream fpattout;
  for (size_t i = 0; i < fNuFluxFilePatterns.size(); ++i)
    fpattout << "\n [" << std::setw(3) << i << "] " << fNuFluxFilePatterns[i];

  std::ostringstream flistout;
  std::vector<std::string> flist = GetFileList();
  for (size_t i = 0; i < flist.size(); ++i)
    flistout << "\n [" << std::setw(3) << i << "] " << flist[i];

  TLorentzVector usr0(0,0,0,0);
  TLorentzVector usr0asbeam;
  User2BeamPos(usr0,usr0asbeam);

  const int w=10, p=6;
  std::ostringstream beamrot_str, beamrotinv_str;
  beamrot_str
    << "fBeamRot: " << std::setprecision(p) << "\n"
    << "  [ "
    << std::setw(w) << fBeamRot.XX() << " "
    << std::setw(w) << fBeamRot.XY() << " "
    << std::setw(w) << fBeamRot.XZ() << " ]\n"
    << "  [ "
    << std::setw(w) << fBeamRot.YX() << " "
    << std::setw(w) << fBeamRot.YY() << " "
    << std::setw(w) << fBeamRot.YZ() << " ]\n"
    << "  [ "
    << std::setw(w) << fBeamRot.ZX() << " "
    << std::setw(w) << fBeamRot.ZY() << " "
    << std::setw(w) << fBeamRot.ZZ() << " ]";
  beamrotinv_str
    << "fBeamRotInv: " << std::setprecision(p) << "\n"
    << "  [ "
    << std::setw(w) << fBeamRotInv.XX() << " "
    << std::setw(w) << fBeamRotInv.XY() << " "
    << std::setw(w) << fBeamRotInv.XZ() << " ]\n"
    << "  [ "
    << std::setw(w) << fBeamRotInv.YX() << " "
    << std::setw(w) << fBeamRotInv.YY() << " "
    << std::setw(w) << fBeamRotInv.YZ() << " ]\n"
    << "  [ "
    << std::setw(w) << fBeamRotInv.ZX() << " "
    << std::setw(w) << fBeamRotInv.ZY() << " "
    << std::setw(w) << fBeamRotInv.ZZ() << " ]";

  std::ostringstream config_str;
  config_str
    << "GDk2NuFlux Config:"
    << "\n Enu_max " << fMaxEv
    << "\n pdg-codes: " << s.str() << "\n "
    << fNEntries << " entries"
    << " (FilePOTs " << fFilePOTs << ") "
    <<  "in " << fNFiles << " files: "
    << flistout.str()
    << "\n from file patterns:"
    << fpattout.str()
    << "\n wgt max=" << fMaxWeight
    << " min=" << fMinMaxWeight << " scan=" << fMaxWeightScan
    << " ( fudge=" << fMaxWgtFudge << " using scan of "
    << fMaxWgtEntries << " rays )"
    << "\n exceeded init=" << fMaxWeightInit << ", N=" << fMaxWgtExceeded
    << " times (fail model " << fMaxWgtFailModel << ") max = " << fMaxWeightMax
    << "\n Z0 pushback " << fZ0
    << "\n used entry " << fIEntry << " " << fIUse << "/" << fNUse
    << " times, in " << fICycle << "/" << fNCycles << " cycles"
    << "\n SumWeight " << fSumWeight << " for " << fNNeutrinos << " neutrino entries"
    << "\n EffPOTsPerNu " << fEffPOTsPerNu << " AccumPOTs " << fAccumPOTs
    << "\n GenWeighted: \"" << (fGenWeighted?"true":"false") << "\", "
    << "ApplyTiltWeight: \"" << (fApplyTiltWeight?"true":"false") << "\", "
    << "PartialRead: \"" << (fPartialRead?"true":"false") << "\", "
    << "Detector location set: \"" << (fDetLocIsSet?"true":"false") << "\", "
    << "\n LengthUnits " << fLengthUnits << ", scale b2u " << fLengthScaleB2U
    << ", scale u2b " << fLengthScaleU2B;

  if ( IsFluxSphere() ) {
    config_str
      << "\n Flux Sphere radius: " << fFluxSphereRadius
      << "\n   User Center: " << utils::print::Vec3AsString(&fFluxSphereCenterUser)
      << "\n   Beam Center: " << utils::print::Vec3AsString(&fFluxSphereCenterBeam);
  } else {
    config_str
      << "\n User Flux Window: "
      << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[0]))
      << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[1]))
      << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[2]))
      << "\n Flux Window (cm, beam coord): "
      << "\n  base " << utils::print::X4AsString(&fFluxWindowBase)
      << "\n  dir1 " << utils::print::X4AsString(&fFluxWindowDir1) << " len " << fFluxWindowLen1
      << "\n  dir2 " << utils::print::X4AsString(&fFluxWindowDir2) << " len " << fFluxWindowLen2;
  }
  config_str
      << "\n User Beam Origin: "
      << "\n  base " << utils::print::X4AsString(&fBeamZero)
      << "\n " << beamrot_str.str() << " "
      << "\n Detector Origin (beam coord): "
      << "\n  base " << utils::print::X4AsString(&usr0asbeam)
      << "\n " << beamrotinv_str.str();

  LOG("Flux", pNOTICE) << config_str.str();
}

//___________________________________________________________________________
std::vector<std::string> GDk2NuFlux::GetFileList()
{
  std::vector<std::string> flist;
  if ( ! fNuFluxTree ) return flist;
  TObjArray *fileElements=fNuFluxTree->GetListOfFiles();
  if ( ! fileElements ) return flist;
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    flist.push_back(chEl->GetTitle());
  }
  return flist;
}

//___________________________________________________________________________

std::vector<double> GDk2NuFluxXMLHelper::GetDoubleVector(std::string str)
{
  // turn string into vector<double>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str," ,;:()[]=");
  std::vector<double> vect;
  size_t ntok = strtokens.size();

  if ( fVerbose > 2 )
    std::cout << "GetDoubleVector \"" << str << "\"" << std::endl;

  for (size_t i=0; i < ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if ( " " == trimmed || "" == trimmed ) continue;  // skip empty strings
    double val = strtod(trimmed.c_str(), (char**)NULL);
    if ( fVerbose > 2 )
      std::cout << "(" << vect.size() << ") = " << val << std::endl;
    vect.push_back(val);
  }

  return vect;
}

std::vector<long int> GDk2NuFluxXMLHelper::GetIntVector(std::string str)
{
  // turn string into vector<long int>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str," ,;:()[]=");
  std::vector<long int> vect;
  size_t ntok = strtokens.size();

  if ( fVerbose > 2 )
    std::cout << "GetIntVector \"" << str << "\"" << std::endl;

  for (size_t i=0; i < ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if ( " " == trimmed || "" == trimmed ) continue;  // skip empty strings
    long int val = strtol(trimmed.c_str(),(char**)NULL,10);
    if ( fVerbose > 2 )
      std::cout << "(" << vect.size() << ") = " << val << std::endl;
    vect.push_back(val);
  }
  return vect;
}

bool GDk2NuFluxXMLHelper::LoadConfig(std::string cfg)
{
  std::string fname = utils::xml::GetXMLFilePath(fGDk2NuFlux->GetXMLFileBase());

  bool is_accessible = ! (gSystem->AccessPathName(fname.c_str()));
  if (!is_accessible) {
    SLOG("GDk2NuFlux", pERROR)
      << "The XML doc doesn't exist! (filename: " << fname << ")";
    return false;
  }

  xmlDocPtr xml_doc = xmlParseFile( fname.c_str() );
  if ( xml_doc == NULL) {
    SLOG("GDk2NuFlux", pERROR)
      << "The XML doc can't be parsed! (filename: " << fname << ")";
    return false;
  }

  xmlNodePtr xml_root = xmlDocGetRootElement( xml_doc );
  if ( xml_root == NULL ) {
    SLOG("GDk2NuFlux", pERROR)
      << "The XML doc is empty! (filename: " << fname << ")";
    return false;
  }

  std::string rootele    = "gdk2nu_config";
  std::string rooteleAlt = "gnumi_config";  // read older GNuMIFlux files
  if ( xmlStrcmp(xml_root->name, (const xmlChar*)rootele.c_str() ) &&
       xmlStrcmp(xml_root->name, (const xmlChar*)rooteleAlt.c_str() ) ) {
    SLOG("GDk2NuFlux", pERROR)
      << "The XML doc has invalid root element! (filename: " << fname << ")"
      << " expected \"" << rootele << "\", saw \"" << xml_root->name << "\"";
    return false;
  }

  SLOG("GDk2NuFlux", pINFO) << "Attempt to load config \"" << cfg
                           << "\" from file: " << fname;

  bool found = this->LoadParamSet(xml_doc,cfg);

  xmlFree(xml_doc);
  return found;

}

bool GDk2NuFluxXMLHelper::LoadParamSet(xmlDocPtr& xml_doc, std::string cfg)
{

  xmlNodePtr xml_root = xmlDocGetRootElement( xml_doc );

  // loop over all xml tree nodes that are children of the root node
  // read the entries looking for "param_set" of the right name

  // loop looking for particular config
  bool found = false;
  xmlNodePtr xml_pset = xml_root->xmlChildrenNode;
  for ( ; xml_pset != NULL ; xml_pset = xml_pset->next ) {
    if ( ! xmlStrEqual(xml_pset->name, (const xmlChar*)"param_set") ) continue;
    // every time there is a 'param_set' tag
    std::string param_set_name =
      utils::str::TrimSpaces(utils::xml::GetAttribute(xml_pset,"name"));

    if ( param_set_name != cfg ) continue;

    SLOG("GDk2NuFlux", pINFO) << "Found config \"" << cfg;

    this->ParseParamSet(xml_doc,xml_pset);
    found = true;

  } // loop over elements of root
  xmlFree(xml_pset);

  return found;
}

void GDk2NuFluxXMLHelper::ParseParamSet(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // handle basic gnumi_config/param_set
    // bad cast away const on next line, but function sig requires it
    std::string pname =
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( pname == "text" || pname == "comment" ) continue;
    std::string pval  =
      utils::xml::TrimSpaces(
              xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));

    if ( fVerbose > 1 )
      SLOG("GDk2NuFlux", pINFO)
        << "   pname \"" << pname << "\", string value \"" << pval << "\"";

    if        ( pname == "verbose" ) {
      fVerbose = atoi(pval.c_str());

    } else if ( pname == "using_param_set" ) {
      SLOG("GDk2NuFlux", pWARN) << "start using_param_set: \"" << pval << "\"";
      bool found = this->LoadParamSet(xml_doc,pval); // recurse
      if ( ! found ) {
        SLOG("GDk2NuFlux", pFATAL) << "using_param_set: \"" << pval << "\" NOT FOUND";
        assert(found);
      }
      SLOG("GDk2NuFlux", pWARN) << "done using_param_set: \"" << pval << "\"";
    } else if ( pname == "units" ) {
      // this is legacy ... equivalent to <lunits> for length
      double scale = genie::utils::units::UnitFromString(pval);
      fGDk2NuFlux->SetLengthUnits(scale);
      SLOG("GDk2NuFlux", pINFO) << "set user length units to \"" << pval << "\"";

    } else if ( pname == "lunits" ) {
      double scale = genie::utils::units::UnitFromString(pval);
      fGDk2NuFlux->SetLengthUnits(scale);
      SLOG("GDk2NuFlux", pINFO) << "set user length units to \"" << pval << "\"";

    } else if ( pname == "tunits" ) {
      double scale = genie::utils::units::UnitFromString(pval);
      fGDk2NuFlux->SetTimeUnits(scale);
      SLOG("GDk2NuFlux", pINFO) << "set user time units to \"" << pval << "\"";

    } else if ( pname == "beamdir" ) {
      ParseBeamDir(xml_doc,xml_child);
      fGDk2NuFlux->SetBeamRotation(fBeamRotXML);

    } else if ( pname == "beampos" ) {
      ParseBeamPos(pval);
      fGDk2NuFlux->SetBeamCenter(fBeamPosXML);

    } else if ( pname == "window" ) {
      ParseWindowSeries(xml_doc,xml_child);
      // RWH  !!!! MEMORY LEAK!!!!
      //std::cout << " flux window " << std::endl
      //          << " [0] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[0],0)) << std::endl
      //          << " [1] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[1],0)) << std::endl
      //          << " [2] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[2],0)) << std::endl;

      fGDk2NuFlux->SetFluxWindow(fFluxWindowPtXML[0],
                            fFluxWindowPtXML[1],
                            fFluxWindowPtXML[2]);

    } else if ( pname == "enumax" ) {
      ParseEnuMax(pval);

    } else if ( pname == "minmaxwgt" ) {
      double minmaxwgt = 0;
      std::vector<double> v = GetDoubleVector(pval);
      if ( v.size() > 0 ) minmaxwgt = v[0];
      fGDk2NuFlux->SetMinMaxWeight(minmaxwgt);
      SLOG("GDk2NuFlux", pINFO) << "set minmaxwgt = " << minmaxwgt;

    } else if ( pname == "maxwgtfail" ) {
      ParseMaxWgtFail(pval);

    } else if ( pname == "upstreamz" ) {
      double z0usr = -3.4e38;
      std::vector<double> v = GetDoubleVector(pval);
      if ( v.size() > 0 ) z0usr = v[0];
      fGDk2NuFlux->SetUpstreamZ(z0usr);
      SLOG("GDk2NuFlux", pINFO) << "set upstreamz = " << z0usr;

    } else if ( pname == "reuse" ) {
      long int nreuse = 1;
      std::vector<long int> v = GetIntVector(pval);
      if ( v.size() > 0 ) nreuse = v[0];
      fGDk2NuFlux->SetEntryReuse(nreuse);
      SLOG("GDk2NuFlux", pINFO) << "set entry reuse = " << nreuse;

    } else if ( pname == "partialread" ) {
      // read only bsim::Decay until a ray is accepted?
      std::string v = utils::str::ToLower(pval);
      if        ( v == "true"  || v == "yes" || v == "on"  || v == "1" ) {
        fGDk2NuFlux->SetPartialRead(true);
      } else if ( v == "false" || v == "no"  || v == "off" || v == "0" ) {
        fGDk2NuFlux->SetPartialRead(false);
      } else {
        SLOG("GDk2NuFlux", pWARN)
          << "unrecognized <partialread> value \"" << pval
          << "\" (use true/false); left as "
          << ( fGDk2NuFlux->GetPartialRead() ? "true" : "false" );
      }
      SLOG("GDk2NuFlux", pINFO) << "set partialread = "
        << ( fGDk2NuFlux->GetPartialRead() ? "true" : "false" );

    } else {
      SLOG("GDk2NuFlux", pWARN)
        << "  NOT HANDLED: pname \"" << pname
        << "\", string value \"" << pval << "\"";

    }

  } // loop over param_set contents
  xmlFree(xml_child);
}

void GDk2NuFluxXMLHelper::ParseBeamDir(xmlDocPtr& xml_doc, xmlNodePtr& xml_beamdir)
{
  fBeamRotXML.SetToIdentity(); // start fresh

  std::string dirtype =
    utils::str::TrimSpaces(
      utils::xml::GetAttribute(xml_beamdir,"type"));

  std::string pval  =
    utils::xml::TrimSpaces(
      xmlNodeListGetString(xml_doc, xml_beamdir->xmlChildrenNode, 1));

  if        ( dirtype == "series" ) {
    // series of rotations around an axis
    ParseRotSeries(xml_doc,xml_beamdir);

  } else if ( dirtype == "thetaphi3") {
    // G3 style triplet of (theta,phi) pairs
    std::vector<double> thetaphi3 = GetDoubleVector(pval);
    std::string units =
      utils::str::TrimSpaces(utils::xml::GetAttribute(xml_beamdir,"units"));
    if ( thetaphi3.size() == 6 ) {
      TRotation fTempRot;
      TVector3 newX = AnglesToAxis(thetaphi3[0],thetaphi3[1],units);
      TVector3 newY = AnglesToAxis(thetaphi3[2],thetaphi3[3],units);
      TVector3 newZ = AnglesToAxis(thetaphi3[4],thetaphi3[5],units);
      fTempRot.RotateAxes(newX,newY,newZ);
      fBeamRotXML = fTempRot;  //.Inverse();
    } else {
      SLOG("GDk2NuFlux", pWARN)
        << " type=\"" << dirtype << "\" within <beamdir> needs 6 values";
    }

  } else if ( dirtype == "newxyz" ) {
    // G4 style new axis values
    std::vector<double> newdir = GetDoubleVector(pval);
    if ( newdir.size() == 9 ) {
      TRotation fTempRot;
      TVector3 newX = TVector3(newdir[0],newdir[1],newdir[2]).Unit();
      TVector3 newY = TVector3(newdir[3],newdir[4],newdir[5]).Unit();
      TVector3 newZ = TVector3(newdir[6],newdir[7],newdir[8]).Unit();
      fTempRot.RotateAxes(newX,newY,newZ);
      fBeamRotXML = fTempRot.Inverse(); // weirdly necessary: frame vs. obj rot
    } else {
      SLOG("GDk2NuFlux", pWARN)
        << " type=\"" << dirtype << "\" within <beamdir> needs 9 values";
    }

  } else {
    // yet something else ... what? 3 choices weren't sufficient?
    SLOG("GDk2NuFlux", pWARN)
      << " UNHANDLED type=\"" << dirtype << "\" within <beamdir>";
  }

  if ( fVerbose > 1 ) {
    int w=10, p=6;
    std::cout << " fBeamRotXML: " << std::setprecision(p) << std::endl;
    std::cout << " [ "
              << std::setw(w) << fBeamRotXML.XX() << " "
              << std::setw(w) << fBeamRotXML.XY() << " "
              << std::setw(w) << fBeamRotXML.XZ() << std::endl
              << "   "
              << std::setw(w) << fBeamRotXML.YX() << " "
              << std::setw(w) << fBeamRotXML.YY() << " "
              << std::setw(w) << fBeamRotXML.YZ() << std::endl
              << "   "
              << std::setw(w) << fBeamRotXML.ZX() << " "
              << std::setw(w) << fBeamRotXML.ZY() << " "
              << std::setw(w) << fBeamRotXML.ZZ() << " ] " << std::endl;
    std::cout << std::endl;
  }

}

void GDk2NuFluxXMLHelper::ParseBeamPos(std::string str)
{
  std::vector<double> xyz = GetDoubleVector(str);
  if ( xyz.size() == 3 ) {
    fBeamPosXML = TVector3(xyz[0],xyz[1],xyz[2]);
  } else if ( xyz.size() == 6 ) {
    // should check for '=' between triplets but we won't be so pedantic
    // ( userx, usery, userz ) = ( beamx, beamy, beamz )
    TVector3 userpos(xyz[0],xyz[1],xyz[2]);
    TVector3 beampos(xyz[3],xyz[4],xyz[5]);
    fBeamPosXML = userpos - fBeamRotXML*beampos;
  } else {
    SLOG("GDk2NuFlux", pWARN)
      << "Unable to parse " << xyz.size() << " values in <beampos>";
    return;
   }
  if ( fVerbose > 1 ) {
    int w=16, p=10;
    std::cout << " fBeamPosXML: [ " << std::setprecision(p)
              << std::setw(w) << fBeamPosXML.X() << " , "
              << std::setw(w) << fBeamPosXML.Y() << " , "
              << std::setw(w) << fBeamPosXML.Z() << " ] "
              << std::endl;
  }
}

void GDk2NuFluxXMLHelper::ParseEnuMax(std::string str)
{
  std::vector<double> v = GetDoubleVector(str);
  size_t n = v.size();
  if ( n > 0 ) {
    fGDk2NuFlux->SetMaxEnergy(v[0]);
    if ( fVerbose > 1 )
      std::cout << "ParseEnuMax SetMaxEnergy(" << v[0] << ") " << std::endl;
  }
  if ( n > 1 ) {
    fGDk2NuFlux->SetMaxEFudge(v[1]);
    if ( fVerbose > 1 )
      std::cout << "ParseEnuMax SetMaxEFudge(" << v[1] << ")" << std::endl;
  }
  if ( n > 2 ) {
    if ( n == 3 ) {
      fGDk2NuFlux->SetMaxWgtScan(v[2]);
      if ( fVerbose > 1 )
        std::cout << "ParseEnuMax SetMaxWgtScan(" << v[2] << ")" << std::endl;
    } else {
      long int nentries = (long int)v[3];
      fGDk2NuFlux->SetMaxWgtScan(v[2],nentries);
      if ( fVerbose > 1 )
        std::cout << "ParseEnuMax SetMaxWgtScan(" << v[2] << "," << nentries << ")" << std::endl;
    }
  }
}

void GDk2NuFluxXMLHelper::ParseMaxWgtFail(std::string str)
{
  // what to do if we exceed estimate maximum weight
  // 0 = "bump"
  // 1 = "frozen"  treat as frozen:  accept and move on
  // 2 = "abort"
  int model = 0;  // default
  if      ( str == "0" || str == "bump"   ) model = 0;
  else if ( str == "1" || str == "frozen" ) model = 1;
  else if ( str == "2" || str == "abort"  ) model = 2;
  else {
    std::cout << "ParseMaxWeightFail \"" << str << "\" unrecognized "
              << std::endl;
  }
  if ( fVerbose > 1 ) {
    std::cout << "ParseMaxWeightFail \"" << str << "\" treat as "
              << model << std::endl;
  }
  fGDk2NuFlux->SetMaxWeightFailModel(model);
}

void GDk2NuFluxXMLHelper::ParseRotSeries(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  TRotation fTempRot; // reset matrix

  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // in a <beamdir> of type "series"
    // should be a sequence of <rotation> entries
    std::string name =
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "rotation" ) {
      std::string val = utils::xml::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      std::string axis =
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"axis"));

      std::string units =
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"units"));

      double rot = atof(val.c_str());
      // assume radians unless given a hint that it's degrees
      if ( 'd' == units[0] || 'D' == units[0] ) rot *= TMath::DegToRad();

      if ( fVerbose > 0 )
        SLOG("GDk2NuFlux", pINFO)
          << " rotate " << rot << " radians around " << axis << " axis";

      if      ( axis[0] == 'x' || axis[0] == 'X' ) fTempRot.RotateX(rot);
      else if ( axis[0] == 'y' || axis[0] == 'Y' ) fTempRot.RotateY(rot);
      else if ( axis[0] == 'z' || axis[0] == 'Z' ) fTempRot.RotateZ(rot);
      else {
        SLOG("GDk2NuFlux", pINFO)
          << " no " << axis << " to rotate around";
      }

    } else {
      SLOG("GDk2NuFlux", pWARN)
        << " found <" << name << "> within <beamdir type=\"series\">";
    }
  }
  // TRotation rotates objects not frames, so we want the inverse
  fBeamRotXML = fTempRot.Inverse();
  xmlFree(xml_child);
}

void GDk2NuFluxXMLHelper::ParseWindowSeries(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  int ientry = -1;

  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // in a <windowr> element
    // should be a sequence of <point> entries
    std::string name =
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "point" ) {
      std::string val  =
        utils::xml::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      std::string coord =
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"coord"));

      std::vector<double> xyz = GetDoubleVector(val);
      if ( xyz.size() != 3 || coord != "det" ) {
        SLOG("GDk2NuFlux", pWARN)
          << "parsing <window> found <point> but size=" << xyz.size()
          << " (expect 3) and coord=\"" << coord << "\" (expect \"det\")"
          << " IGNORE problem";
      }
      ++ientry;
      if ( ientry < 3 && ientry >= 0 ) {
        TVector3 pt(xyz[0],xyz[1],xyz[2]);
        if ( fVerbose > 0 ) {
          int w=16, p=10;
          std::cout << " point[" << ientry <<"] = [ " << std::setprecision(p)
                    << std::setw(w) << pt.X() << " , "
                    << std::setw(w) << pt.Y() << " , "
                    << std::setw(w) << pt.Z() << " ] "
                    << std::endl;
        }
        fFluxWindowPtXML[ientry] = pt;  // save the point
      } else {
        SLOG("GDk2NuFlux", pWARN)
          << " <window><point> ientry " << ientry << " out of range (0-2)";
      }

    } else {
      SLOG("GDk2NuFlux", pWARN)
        << " found <" << name << "> within <window>";
    }
  }
  xmlFree(xml_child);
}

TVector3 GDk2NuFluxXMLHelper::AnglesToAxis(double theta, double phi, std::string units)
{
  double xyz[3];
  // assume radians unless given a hint that it's degrees
  double scale = ('d'==units[0]||'D'==units[0]) ? TMath::DegToRad() : 1.0 ;

  xyz[0] = TMath::Cos(scale*phi)*TMath::Sin(scale*theta);
  xyz[1] = TMath::Sin(scale*phi)*TMath::Sin(scale*theta);
  xyz[2] = TMath::Cos(scale*theta);
  // condition vector to eliminate most floating point errors
  for (int i=0; i<3; ++i) {
    const double eps = 1.0e-15;
    if (TMath::Abs(xyz[i])   < eps ) xyz[i] =  0;
    if (TMath::Abs(xyz[i]-1) < eps ) xyz[i] =  1;
    if (TMath::Abs(xyz[i]+1) < eps ) xyz[i] = -1;
  }
  return TVector3(xyz[0],xyz[1],xyz[2]);
}

TVector3 GDk2NuFluxXMLHelper::ParseTV3(const std::string& str)
{
  std::vector<double> xyz = GetDoubleVector(str);
  if ( xyz.size() != 3 ) {
    return TVector3();
    SLOG("GDk2NuFlux", pWARN)
      << " ParseTV3 \"" << str << "\" had " << xyz.size() << " elements ";
  }
  return TVector3(xyz[0],xyz[1],xyz[2]);

}
//___________________________________________________________________________
