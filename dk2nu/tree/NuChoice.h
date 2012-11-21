///----------------------------------------------------------------------------
/**
 * \class bsim::NuChoice
 * \file  NuChoice.h
 *
 * \brief A class that defines the final choice of the neutrino ray
 *        including it position and momentum in both detector and beam
 *        coordinate systems.  This class is generally used in 
 *        conjunction with the "dk2nu" class.
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2012-11-21 04:47:28 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: NuChoice.h,v 1.1 2012-11-21 04:47:28 rhatcher Exp $
 *
 */
///----------------------------------------------------------------------------

#ifndef BSIM_NUCHOICE_H
#define BSIM_NUCHOICE_H

#include "TLorentzVector.h"
#include "TVector3.h"

namespace bsim {
  /**
   *  All the data members are public as these classes are used as
   *  generalized structs.  As they will be branches of a TTree no
   *  specialized naming indicators signifying that they are member data
   *  of a class will be used, nor will any fancy capitalization schemes.
   *
   *  All classes must implement a clear() method that resets their values
   *  to an identifiably invalid state or clears any vectors.  Additionally
   *  classes should provide a AsString() method for formatting themselves
   *  for use output.
   */
  
  ///---------------------------------------------------------------------------
  /**
   *============================================================================
   *  
   */
  class NuChoice
  {
  public:
    int             pdgNu;       ///< generated nu pdg code
    double          xyWgt;       ///< generated nu x-y weight
    double          impWgt;      ///< original importance weight
                                 ///  GDk2NuFlux::Weight() might be the product of these
    TLorentzVector  p4NuBeam;    ///< generated nu 4-momentum in beam coord
    TLorentzVector  x4NuBeam;    ///< generated nu 4-position in beam coord
    TLorentzVector  p4NuUser;    ///< generated nu 4-momentum in user/det coord
    TLorentzVector  x4NuUser;    ///< generated nu 4-position in user/det coord
    
  public:
    NuChoice();
    virtual     ~NuChoice();
    void        clear(const std::string &opt = "");    ///< reset everything
    std::string AsString(const std::string& opt = "") const;
    
  private:
    ClassDef(bsim::NuChoice,1)
  };  // end-of-class bsim::NuChoice

} // end-of-namespace "bsim"

// not part of namespace bsim
std::ostream& operator<<(std::ostream& os, const bsim::NuChoice& nuchoice);

#endif
