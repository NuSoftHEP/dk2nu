/**
 * \class dk2nu
 * \file  dk2nu.cc
 *
 * \brief A class that defines the "dk2nu" object used as the primary
 *        branch for a TTree for the output of neutrino flux simulations
 *        such as g4numi, g4numi_flugg, etc.
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.3 $
 *
 * \date $Date: 2012-11-07 21:12:01 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: dk2nu.cc,v 1.3 2012-11-07 21:12:01 rhatcher Exp $
 */

#include "dk2nu.h"
ClassImp(dk2nu)

#include <iostream>
#include <iomanip>
#include <sstream>

//-----------------------------------------------------------------------------
dk2nu::dk2nu() { Clear(); }
dk2nu::~dk2nu() { ; }
//-----------------------------------------------------------------------------
void dk2nu::Clear(const std::string &)
{
  const int    kUnsetInt    = -1;
  const double kUnsetDouble = -1.e4;

  job         = kUnsetInt;
  potnum      = kUnsetInt;

  norig       = kUnsetInt;
  ndecay      = kUnsetInt;
  ntype       = kUnsetInt;
  vx          = kUnsetDouble;
  vy          = kUnsetDouble;
  vz          = kUnsetDouble;
  pdpx        = kUnsetDouble;
  pdpy        = kUnsetDouble;
  pdpz        = kUnsetDouble;
  ppdxdz      = kUnsetDouble;
  ppdydz      = kUnsetDouble;
  pppz        = kUnsetDouble;
  ppenergy    = kUnsetDouble;
  ppmedium    = kUnsetDouble;
  ptype       = kUnsetInt;
  ppvx        = kUnsetDouble;
  ppvy        = kUnsetDouble;
  ppvz        = kUnsetDouble;
  muparpx     = kUnsetDouble;
  muparpy     = kUnsetDouble;
  muparpz     = kUnsetDouble;
  mupare      = kUnsetDouble;
  necm        = kUnsetDouble;
  nimpwt      = kUnsetDouble;

  xpoint      = kUnsetDouble;
  ypoint      = kUnsetDouble;
  zpoint      = kUnsetDouble;
  tvx         = kUnsetDouble;
  tvy         = kUnsetDouble;
  tvz         = kUnsetDouble;
  tpx         = kUnsetDouble;
  tpy         = kUnsetDouble;
  tpz         = kUnsetDouble;
  tptype      = kUnsetInt;
  tgen        = kUnsetInt;
  tgptype     = kUnsetInt;

  apdg.clear();
  trackid.clear();
  parentid.clear();
  startx.clear();
  starty.clear();
  startz.clear();
  startpx.clear();
  startpy.clear();
  startpz.clear();

  stopx.clear();
  stopy.clear();
  stopz.clear();
  stoppx.clear();
  stoppy.clear();
  stoppz.clear();

  pprodpx.clear();
  pprodpy.clear();
  pprodpz.clear();

  ivol.clear();
  fvol.clear();
  proc.clear();

  nupx.clear();
  nupy.clear();
  nupz.clear();
  nuenergy.clear();
  nuwgt.clear();

  trkx.clear();
  trky.clear();
  trkz.clear();
  trkpx.clear();
  trkpy.clear();
  trkpz.clear();

  vint.clear();
  vdbl.clear();

}

std::string dk2nu::AsString(const std::string& opt) const
{

  std::ostringstream s;
  s << "dk2nu: \"" << opt << "\" job " << job << " pot# " << potnum << "\n";
  size_t nloc = nuwgt.size();
  bool printloc = true; // ( opt.find("l") != std::string::npos);
  s << nloc << " locations " << (printloc?"":"(suppressed)") << "\n";
  if ( ! printloc ) nloc = 0;
  for ( size_t iloc = 0; iloc < nloc; ++iloc ) {
    s << "[" << std::setw(2) << iloc << "]"
      << " p4=(" << std::setw(12) << nupx[iloc] << "," 
      << std::setw(12) << nupy[iloc] << "," 
      << std::setw(12) << nupz[iloc] << ")"
      << " E=" << std::setw(12) << nuenergy[iloc]
      << " wgt=" << std::setw(12) << nuwgt[iloc]
      << "\n";
  }
  s << "norig " << norig << " ndecay " << ndecay << " ntype " << ntype << "\n";

  s << "production vtx=(" << vx << "," << vy << "," << vz << ")\n";
  return s.str();
}

std::ostream& operator<<(std::ostream& os, const dk2nu& entry)
{
  os << entry.AsString() << std::endl;
  return os;
}
