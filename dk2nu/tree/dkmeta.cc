/**
 * \class dkmeta
 * \file  dkmeta.cc
 *
 * \brief A class that defines the "dkmeta" object used as the primary
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
 * $Id: dkmeta.cc,v 1.1 2012-11-07 01:35:47 rhatcher Exp $
 */

#include "dkmeta.h"
ClassImp(dkmeta)

#include <iostream>

//-----------------------------------------------------------------------------
dkmeta::dkmeta() { Clear(); }
dkmeta::~dkmeta() { ; }
//-----------------------------------------------------------------------------
void dkmeta::Clear(const std::string &)
{
  //const int    kUnsetInt    = -1;
  const double kUnsetDouble = -1.e4;
  const std::string kUnsetString = "<<unset-string>>";

  beamsim     = kUnsetString;
  physics     = kUnsetString;
  physcuts    = kUnsetString;
  tgtcfg      = kUnsetString;
  horncfg     = kUnsetString;
  dkvolcfg    = kUnsetString;

  beam0x      = kUnsetDouble;
  beam0y      = kUnsetDouble;
  beam0z      = kUnsetDouble;
  beamhwidth  = kUnsetDouble;
  beamvwidth  = kUnsetDouble;
  beamdxdz    = kUnsetDouble;
  beamdydz    = kUnsetDouble;
  
  xloc.clear();
  yloc.clear();
  zloc.clear();

  nameloc.clear();

  vintnames.clear();
  vdblnames.clear();
}

std::string dkmeta::AsString(const std::string& opt) const
{
  std::string s("dkmeta: ");
  s += opt;
  return s;
}

std::ostream& operator<<(std::ostream& os, const dkmeta& entry)
{
  os << entry.AsString() << std::endl;
  return os;
}
