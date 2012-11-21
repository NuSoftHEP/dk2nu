/**
 * \class NuChoice
 * \file  NuChoice.cc
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
 * $Id: NuChoice.cc,v 1.1 2012-11-21 04:47:28 rhatcher Exp $
 *
 */

#include "NuChoice.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "TMath.h"

//-----------------------------------------------------------------------------
ClassImp(bsim::NuChoice)
bsim::NuChoice::NuChoice() { clear(); }
bsim::NuChoice::~NuChoice() { ; }
void bsim::NuChoice::clear(const std::string &)
{ 
  pdgNu  = 0;
  xyWgt  = 0;
  impWgt = 0;
  static TLorentzVector NullLV = TLorentzVector(0,0,0,0);
  p4NuBeam = NullLV;
  x4NuBeam = NullLV;
  p4NuUser = NullLV;
  x4NuUser = NullLV;
}
std::string bsim::NuChoice::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "NuChoice: " << pdgNu << " xyWgt " << xyWgt 
    << " impWgt " << impWgt << "\n";
  s << "beam  p3={ " << std::setw(12) << p4NuBeam.Px()
    << "," << std::setw(12) << p4NuBeam.Py()
    << "," << std::setw(12) << p4NuBeam.Pz()
    << "; E=" << std::setw(12) << p4NuBeam.E()
    << "}\n";
  s << "      x3={ " << std::setw(12) << x4NuBeam.X()
    << "," << std::setw(12) << x4NuBeam.Y()
    << "," << std::setw(12) << x4NuBeam.Z()
    << "; t=" << std::setw(12) << x4NuBeam.T()
    << "}\n";
  s << "user  p3={ " << std::setw(12) << p4NuUser.Px()
    << "," << std::setw(12) << p4NuUser.Py()
    << "," << std::setw(12) << p4NuUser.Pz()
    << "; E=" << std::setw(12) << p4NuBeam.E()
    << "}\n";
  s << "      x3={ " << std::setw(12) << x4NuUser.X()
    << "," << std::setw(12) << x4NuUser.Y()
    << "," << std::setw(12) << x4NuUser.Z()
    << "; t=" << std::setw(12) << x4NuUser.T()
    << "}";  // no \n on last line

  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::NuChoice& nuchoice)
{
  os << nuchoice.AsString();
  return os;
}
