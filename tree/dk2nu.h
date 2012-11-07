/**
 * \class dk2nu
 * \file  dk2nu.h
 *
 * \brief A class that defines the "dk2nu" object used as the primary
 *        branch for a TTree for the output of neutrino flux simulations
 *        such as g4numi, g4numi_flugg, etc.
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2012-11-07 01:35:47 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: dk2nu.h,v 1.1 2012-11-07 01:35:47 rhatcher Exp $
 *
 * Notes tagged with "DK2NU" are questions that should be answered
 */

#ifndef DK2NU_H
#define DK2NU_H

#include "TROOT.h"
#include "TObject.h"

#include <vector>
#include <string>

class dk2nu;
std::ostream& operator<<(std::ostream& os, const dk2nu& entry);

class dk2nu
{
private:
  ClassDef(dk2nu,3) // KEEP THIS UP-TO-DATE!  increment for each change

public:
  /**
   *   Public methods for constructing/destruction and resetting the data
   */
  dk2nu();
  virtual     ~dk2nu();
  void        Clear(const std::string &opt = "");    ///< reset everything
  std::string AsString(const std::string& opt = "") const; ///< output as a string

  /**
   *  All the data members are public as this class is used as a
   *  generalized struct, with just the addition of the Clear() method.
   *  As they will be branches of a TTree no specialized naming 
   *  indicators signifying that they are member data of a class 
   *  will be used, nor will any fancy capitalization schemes.
   */

  /**
   *=======================================================================
   *  General Info
   */
   Int_t job;              ///< identifying job #
   Int_t potnum;           ///< proton # processed by simulation

  /**
   *=======================================================================
   *  Fixed Decays:
   *  A random ray plus those directed at specific points.
   */
   std::vector<Double_t> nupx;     ///< px for nu at location(s)
   std::vector<Double_t> nupy;     ///< py for nu at location(s)
   std::vector<Double_t> nupz;     ///< pz for nu at location(s)
   std::vector<Double_t> nuenergy; ///< E for nu at location(s)
   std::vector<Double_t> nuwgt;    ///< weight for nu at location(s)

  /**
   *=======================================================================
   *  Decay Data:
   *  Core information about the neutrino and the decay that gave rise to it.
   *  % = necessary for reweighting
   */
   Int_t    norig;        ///< not used?
   Int_t    ndecay;       ///< decay process (see dkproc_t)
   Int_t    ntype;        ///< % neutrino flavor (PDG? code)

   Double_t vx;           ///< % neutrino production vertex x
   Double_t vy;           ///< % neutrino production vertex y
   Double_t vz;           ///< % neutrino production vertex z
   Double_t pdpx;         ///< % px momentum of nu parent at (vx,vy,vz)
   Double_t pdpy;         ///< % py momentum of nu parent at (vx,vy,vz)
   Double_t pdpz;         ///< % pz momentum of nu parent at (vx,vy,vz)

   /**  these are used in muon decay case? */
   Double_t ppdxdz;       ///< % direction of nu parent at its production point
   Double_t ppdydz;       ///< % direction of nu parent at its production point
   Double_t pppz;         ///< % z momentum of nu parent at its production point
   Double_t ppenergy;     ///< % energy of nu parent at its production point

   Double_t ppmedium;     ///< material nu parent was produced in
   Int_t    ptype;        ///< % nu parent species (PDG? code)

   /** momentum and energy of nu grandparent at
       muons:    grandparent decay point
       hadrons:  grandparent production point
       Huh?  this needs better documentation
    */
   Double_t muparpx;      ///< %
   Double_t muparpy;      ///< %
   Double_t muparpz;      ///< %
   Double_t mupare;       ///< % energy of nu grandparent

   Double_t necm;         ///< % nu energy in center-of-mass frame
   Double_t nimpwt;       ///< % production vertex z of nu parent

  /**
   *=======================================================================
   *  (Grand)Parent Info:
   *
   */

   /**
    * DK2NU: are these needed for any/all cases?
    */
   Double_t ppvx;         ///< production vertex x of nu parent
   Double_t ppvy;         ///< production vertex y of nu parent
   Double_t ppvz;         ///< production vertex z of nu parent

   /**
    * DK2NU: do we need these?  these aren't filled by flugg, others?
    */
   Double_t xpoint;       ///< ?
   Double_t ypoint;       ///< ?
   Double_t zpoint;       ///< ?

   /**
    * these ancestors are possibly, but not necessarily, the direct nu parent
    * DK2NU: can these be removed in favor of cascade info below?
    */
   Double_t tvx;          ///< x position of nu ancestor as it exits target
   Double_t tvy;          ///< y position of nu ancestor as it exits target
   Double_t tvz;          ///< z position of nu ancestor as it exits target
   Double_t tpx;          ///< x momentum of nu ancestor as it exits target
   Double_t tpy;          ///< y momentum of nu ancestor as it exits target
   Double_t tpz;          ///< z momentum of nu ancestor as it exits target
   Int_t    tptype;       ///< species of ancestor exiting the target
   Int_t    tgen;         ///< nu parent generation in cascade:
                          ///<   1=primary proton
                          ///<   2=particles produced by proton interaction
                          ///<   etc
   /**
    * these are only in g3numi and flugg
    * DK2NU: can these be removed in favor of cascade info below?
    *        for now we'll leave them in place
    */
   Int_t    tgptype;      ///< species of parent of particle exiting the target (PDG code?)

   Double_t tgppx;        ///< x momentum of parent of particle exiting target at the parent production point ...
   Double_t tgppy;        ///< y momentum
   Double_t tgppz;        ///< z momentum
   Double_t tprivx;       ///< primary particle interaction vtx (not used?)
   Double_t tprivy;       ///< primary particle interaction vtx (not used?)
   Double_t tprivz;       ///< primary particle intereaction vtx (not used?)
   Double_t beamx;        ///< primary proton origin
   Double_t beamy;        ///< primary proton origin
   Double_t beamz;        ///< primary proton origin
   Double_t beampx;       ///< primary proton momentum
   Double_t beampy;       ///< primary proton momentum
   Double_t beampz;       ///< primary proton momentum

   /**
    * these are in the g4numi and minerva ntuples
    * DK2NU: but what do they mean and are the duplicative to
    *        the more complete progenitor info below?
    */
   std::vector<Double_t> trkx;
   std::vector<Double_t> trky;
   std::vector<Double_t> trkz;
   std::vector<Double_t> trkpx;
   std::vector<Double_t> trkpy;
   std::vector<Double_t> trkpz;

  /**
   *=======================================================================
   *  Progenitor Info:
   *  Complete ancestral info from primary proton down to decaying particle
   *
   *  DK2NU: this is mainly (based on) the minerva extensions *except*
   *         some names are changed to avoid confusion and
   *         distances will be cm, energies in GeV (unless the whole
   *         record uniformly uses something else and is flagged as such)
   */
   std::vector<Int_t>    apdg;     ///< ancestor species
   std::vector<Int_t>    trackid;  ///< ??? particle trackId
   std::vector<Int_t>    parentid; ///< ??? parentId

   std::vector<Double_t> startx;   ///< particle x initial position
   std::vector<Double_t> starty;   ///< particle y initial position
   std::vector<Double_t> startz;   ///< particle z initial position
   std::vector<Double_t> stopx;    ///< particle x final position
   std::vector<Double_t> stopy;    ///< particle y final position
   std::vector<Double_t> stopz;    ///< particle z final position

   std::vector<Double_t> startpx;  ///< particle x initial momentum
   std::vector<Double_t> startpy;  ///< particle y initial momentum
   std::vector<Double_t> startpz;  ///< particle z initial momentum
   std::vector<Double_t> stoppx;   ///< particle x final momentum
   std::vector<Double_t> stoppy;   ///< particle y final momentum
   std::vector<Double_t> stoppz;   ///< particle z final momentum

   std::vector<Double_t> pprodpx;  ///< parent x momentum when producing this particle, MeV/c
   std::vector<Double_t> pprodpy;  ///< parent y momentum when producing this particle
   std::vector<Double_t> pprodpz;  ///< parent z momentum when producing this particle

   std::vector<std::string> proc;  ///< name of the process that creates this particle

   std::vector<std::string> ivol;  ///< name of the volume where the particle starts
   std::vector<std::string> fvol;  ///< name of the volume where the particle stops

   /**
    *=======================================================================
    *  Special Info:
    */
   Int_t    flagbits;      ///< bits signify non-std setting such as
                           ///< Geant vs. PDG codes, mm vs. cm, Mev vs. GeV
   std::vector<Int_t>    vint;    ///< user defined vector of integers
   std::vector<Double_t> vdbl;    ///< user defined vector of doubles

   /**
    *=======================================================================
    *  Random Info:
    *  blah, blah, blah
    */
   
   Int_t    ptrkid;        ///< lbne addition

   /**
    *=======================================================================
    *  Specialized enumerations
    */

   /**
    *  Proposed flag bits:
    */
   typedef enum flgbitval {
     flg_dist_m           = 0x00000000,  ///< no special bit for meters
     flg_dist_cm          = 0x00020000,  ///< distances in cm (default)
     flg_dist_mm          = 0x00030000,  ///< distances in mm
     flg_e_gev            = 0x00000000,  ///< no special bit for GeV (default)
     flg_e_mev            = 0x00300000,  ///< energies in MeV
     flg_usr_mask         = 0x0000FFFF,
     flg_reserved_mask    = 0xFFFF0000
   } flgbitval_t;

   /**
    *  Enumeration of decay processes, stored in "ndecay"
    *  store as integer; these are for reference
    *  DK2NU:  should there be an associated AsString() method
    *          that returns a text (optionally formatted for latex?)?
    */
   typedef enum dkproc {
     dkp_unknown         =  0,
     dkp_k0l_nuepimep    =  1,  ///< k0long => nu_e + pi- + e+
     dkp_k0l_nuebpipem   =  2,  ///< k0long => nu_e_bar + p+ + e-
     dkp_k0l_numupimmup  =  3,  ///< k0long => nu_mu + pi- + mu+
     dkp_k0l_numubpipmum =  4,  ///< k0long => nu_mu_bar + pi+ + mu-
     dkp_kp_numumup      =  5,  ///< k+ => nu_mu + mu+
     dkp_kp_nuepi0ep     =  6,  ///< k+ => nu_e + pi0 + e+
     dkp_kp_numupi0mup   =  7,  ///< k+ => nu_mu + pi0 + mu+
     dkp_kp_numubmum     =  8,  ///< k- => nu_mu_bar + mu-
     dkp_kp_nuebpi0em    =  9,  ///< k- => nu_e_bar + pi0 + e-
     dkp_kp_numubpi0mum  = 10,  ///< k- => nu_mu_bar + pi0 + mu-
     dkp_mup_nusep       = 11,  ///< mu+ => nu_mu_bar + nu_e + e+
     dkp_mum_nusep       = 12,  ///< mu- => nu_mu + nu_e_bar + e-
     dk_pip_numumup      = 13,  ///< pi+ => nu_mu + mu+
     dk_pim_numubmum     = 14,  ///< pi- => nu_mu_bar + mu-
     dkp_maximum,               ///< one-beyond end for iterating
     dkp_other           = 999  ///< flag for unusual cases
   } dkproc_t;
   
 };

#endif
