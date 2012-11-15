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
 * \version $Revision: 1.2 $
 *
 * \date $Date: 2012-11-15 09:09:27 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: dkmeta.cc,v 1.2 2012-11-15 09:09:27 rhatcher Exp $
 */

#include "dkmeta.h"
#include "dflt.h"

#include <iostream>
#include <iomanip>
#include <sstream>

//-----------------------------------------------------------------------------
ClassImp(bsim::Location)
bsim::Location::Location() { clear(); }
bsim::Location::~Location() { ; }
bsim::Location::Location(double xi, double yi, double zi, std::string namei)
  : x(xi), y(yi), z(zi), name(namei) { ; }
void bsim::Location::clear(const std::string &)
{
  x = y = z = 0;
  name = "<<no-location>>";
}
std::string bsim::Location::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << " {" 
    << std::setw(12) << x << ", "
    << std::setw(12) << y << ", "
    << std::setw(12) << z << "} \"" << name << "\"";
  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::Location& location)
{
  os << location.AsString();
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::DkMeta)
bsim::DkMeta::DkMeta() { clear(); }
bsim::DkMeta::~DkMeta() { ; }
void bsim::DkMeta::clear(const std::string &)
{
  const std::string kUnsetString = "<<unset-string>>";

  job  = bsim::kDfltInt;
  pots = 0;  // initialize to sensible value

  beamsim     = kUnsetString;
  physics     = kUnsetString;
  physcuts    = kUnsetString;
  tgtcfg      = kUnsetString;
  horncfg     = kUnsetString;
  dkvolcfg    = kUnsetString;

  beam0x      = bsim::kDfltDouble;
  beam0y      = bsim::kDfltDouble;
  beam0z      = bsim::kDfltDouble;
  beamhwidth  = bsim::kDfltDouble;
  beamvwidth  = bsim::kDfltDouble;
  beamdxdz    = bsim::kDfltDouble;
  beamdydz    = bsim::kDfltDouble;
  
  location.clear();

  vintnames.clear();
  vdblnames.clear();
}
std::string bsim::DkMeta::AsString(const std::string& opt) const
{
  std::ostringstream s;
  s << "bsim::DkMeta: \"" << opt << "\" job " << job
    << " pots " << pots << "\n";
  /**
   * add print of job meta data strings here
   */
  size_t nl = location.size();
  s << nl << " locations:\n";
  for ( size_t l = 0; l < nl; ++l ) {
    s << " [" << std::setw(2) << l << "] " << location[l] << "\n";
  }
  s << "beamsim:  \"" << beamsim << "\"\n"
    << "physics:  \"" << physics << "\"\n"
    << "physcuts: \"" << physcuts << "\"\n"
    << "tgtcfg:   \"" << tgtcfg << "\"\n"
    << "horncfg:  \"" << horncfg << "\"\n"
    << "dkvolcfg: \"" << dkvolcfg << "\"\n";
  s << "beam 0={" << beam0x << "," << beam0y << "," << beam0z << "}"
    << " dxdz=" << beamdxdz << " dydz=" << beamdydz << "\n"
    << "width {h=" << beamhwidth << ", v=" << beamvwidth << "}\n";
  size_t ni = vintnames.size();
  if ( ni > 0 ) {
    s << ni << " int names: ";
    for ( size_t i = 0; i < ni; ++i ) { s << " \"" << vintnames[i] << "\""; }
    s << "\n";
  }
  size_t nd = vdblnames.size();
  if ( nd > 0 ) {
    s << nd << " dbl names: ";
    for ( size_t i = 0; i < nd; ++i ) { s << " \"" << vdblnames[i] << "\""; }
    s << "\n";
  }

  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::DkMeta& dkmeta)
{
  os << dkmeta.AsString(); // << std::endl;
  return os;
}
