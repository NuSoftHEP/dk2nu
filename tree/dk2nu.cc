/**
 * \class Dk2Nu
 * \file  dk2nu.cc
 *
 * \brief A class that defines the "dk2nu" object used as the primary
 *        branch for a TTree for the output of neutrino flux simulations
 *        such as g4numi, g4numi_flugg, etc.
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.4 $
 *
 * \date $Date: 2012-11-15 09:09:27 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: dk2nu.cc,v 1.4 2012-11-15 09:09:27 rhatcher Exp $
 */

#include "dk2nu.h"
#include "dflt.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "TMath.h"

//-----------------------------------------------------------------------------
ClassImp(bsim::NuRay)
bsim::NuRay::NuRay() { clear(); }
bsim::NuRay::~NuRay() { ; }
bsim::NuRay::NuRay(double pxi, double pyi, double pzi, double Ei, double wgti)
  : px(pxi), py(pyi), pz(pzi), E(Ei), wgt(wgti) { ; }
void bsim::NuRay::clear(const std::string &)
{ 
  px = py = pz = E = wgt = bsim::kDfltDouble;
}
std::string bsim::NuRay::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "p3={ " << std::setw(12) << px
    << "," << std::setw(12) << py 
    << "," << std::setw(12) << pz 
    << "} E=" << std::setw(10) << E
    << " wgt=" << std::setw(12) << wgt;
  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::NuRay& nuray)
{
  os << nuray.AsString();
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::Decay)
bsim::Decay::Decay() { clear(); }
bsim::Decay::~Decay() { ; }
void bsim::Decay::clear(const std::string &)
{ 
  norig = ndecay = ntype = bsim::kDfltInt;
  vx = vy = vz = pdpx = pdpy = pdpz = bsim::kDfltDouble;
  ppdxdz = ppdydz = pppz = ppenergy = bsim::kDfltDouble;
  ppmedium = bsim::kDfltDouble; 
  ptype = bsim::kDfltInt;
  muparpx = muparpy = muparpz = mupare = bsim::kDfltDouble;
  necm = nimpwt = bsim::kDfltDouble;

}
std::string bsim::Decay::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "norig " << norig << " ndecay " << ndecay 
    << " ntype " << ntype << " ptype " << ptype << "\n";
  s << "necm " << necm << " nimpwt " << nimpwt
    << " ppmedium " << ppmedium << "\n";
  s << "v={" << vx << "," << vy << "," << vz << "} ";
  s << "pdp={" << pdpx << "," << pdpy << "," << pdpz << "}";
  // ppdxdz, ppdydz,pppz,ppenergy
  // muparpx,muparpy,muparpz
  //last line shouldn't have endl << "\n";
  
  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::Decay& decay)
{
  os << decay.AsString();
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::Ancestor)
bsim::Ancestor::Ancestor() { clear(); }
bsim::Ancestor::~Ancestor() { ; }
void bsim::Ancestor::clear(const std::string &)
{
  pdg = 0;  // 0 is not a legal PDG code
  startx = starty = startz = bsim::kDfltDouble;
  startpx = startpy = startpz = bsim::kDfltDouble;
  stoppx = stoppy = stoppz = bsim::kDfltDouble;
  polx = poly = polz = bsim::kDfltDouble;
  pprodpx = pprodpy = pprodpz = bsim::kDfltDouble;
  nucleus = 0;  // not a legal PDG code
  proc = "<<no-process>>";
  ivol = "<<no-volume>>";
}
std::string bsim::Ancestor::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "pdg=" << std::setw(5) << pdg << " \"" << proc << "\""
    << " \"" << ivol << "\" " << nucleus << "\n";
  s << " startx {" << startx << "," << starty << "," << startz << "}";
  //last line shouldn't have endl << "\n";
  return s.str();
}
void bsim::Ancestor::SetStartXYZ(Double_t x, Double_t y, Double_t z)
{ startx = x; starty = y; startz = z; }
void bsim::Ancestor::SetStartP(Double_t px, Double_t py, Double_t pz)
{ startpx = px; startpy = py; startpz = pz; }
void bsim::Ancestor::SetStopP(Double_t px, Double_t py, Double_t pz)
{ stoppx = px; stoppy = py; stoppz = pz; }
void bsim::Ancestor::SetPProdP(Double_t px, Double_t py, Double_t pz)
{ pprodpx = px; pprodpy = py; pprodpz = pz; }

Double_t bsim::Ancestor::r() const 
{ return TMath::Sqrt(startx*startx+starty*starty); }
Double_t bsim::Ancestor::startpt() const
{ return TMath::Sqrt(startpx*startpx+startpy*startpy); }
Double_t bsim::Ancestor::startp() const
{ return TMath::Sqrt(startpx*startpx+startpy*startpy+startpz*startpz); }
Double_t bsim::Ancestor::stoppt() const
{ return TMath::Sqrt(stoppx*stoppx+stoppy*stoppy); }
Double_t bsim::Ancestor::stopp() const
{ return TMath::Sqrt(stoppx*stoppx+stoppy*stoppy+stoppz*stoppz); }
Double_t bsim::Ancestor::pprodpt() const
{ return TMath::Sqrt(pprodpx*pprodpx+pprodpy*pprodpy); }
Double_t bsim::Ancestor::pprodp() const
{ return TMath::Sqrt(pprodpx*pprodpx+pprodpy*pprodpy+pprodpz*pprodpz); }

std::ostream& operator<<(std::ostream& os, const bsim::Ancestor& ancestor)
{
  os << ancestor.AsString();
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::TgtExit)
bsim::TgtExit::TgtExit() { clear(); }
bsim::TgtExit::~TgtExit() { ; }
void bsim::TgtExit::clear(const std::string &)
{
  tvx = tvy = tvz = tpx = tpy = tpz = bsim::kDfltDouble;
  tptype = tgen = bsim::kDfltInt;
}
std::string bsim::TgtExit::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "tgtexit: tptype=" << tptype << " tgen=" << tgen;
  s << " v={" << tvx << "," << tvy << "," << tvz << "}\n";
  s << " p={" << tpx << "," << tpy << "," << tpz << "}";
  //last line shouldn't have endl; << "\n";
  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::TgtExit& tgtexit)
{
  os << tgtexit.AsString();
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::Traj)
bsim::Traj::Traj() { clear(); }
bsim::Traj::~Traj() { ; }
void bsim::Traj::clear(const std::string &)
{
  trkx = trky = trkz = trkpx = trkpy = trkpz = bsim::kDfltDouble;
}
std::string bsim::Traj::AsString(const std::string& /* opt */) const
{
  std::ostringstream s;
  s << "bsim::Traj: ";
  s << " v={" << trkx << "," << trky << "," << trkz << "}";
  s << " p={" << trkpx << "," << trkpy << "," << trkpz << "}";
  //last line shouldn't have endl; << "\n";
  return s.str();
}
std::ostream& operator<<(std::ostream& os, const bsim::Traj& traj)
{
  os << traj.AsString() << std::endl;
  return os;
}

//-----------------------------------------------------------------------------
ClassImp(bsim::Dk2Nu)
bsim::Dk2Nu::Dk2Nu() { clear(); }
bsim::Dk2Nu::~Dk2Nu() { ; }
void bsim::Dk2Nu::clear(const std::string &)
{ 
  job    = bsim::kDfltInt;
  potnum = 0;
  nuray.clear();     /// clear the vector
  decay.clear();     /// clear the object
  ancestor.clear();  /// clear the vector

  ppvx  = bsim::kDfltDouble;
  ppvy  = bsim::kDfltDouble;
  ppvz  = bsim::kDfltDouble;

  tgtexit.clear();  /// clear the object
  traj.clear();     /// clear the vector

  flagbits = 0;
  vint.clear();     /// clear the vector
  vdbl.clear();     /// clear the vector

}
std::string bsim::Dk2Nu::AsString(const std::string& opt) const
{
  std::ostringstream s;
  s << "bsim::Dk2Nu: \"" << opt << "\" job " << job 
    << " pot# " << potnum << "\n";

  size_t nloc = nuray.size();
  bool printloc = true; // ( opt.find("l") != std::string::npos);
  s << nloc << " locations " << (printloc?"":"(suppressed)") << "\n";
  if ( ! printloc ) nloc = 0;
  for ( size_t iloc = 0; iloc < nloc; ++iloc ) {
    s << "[" << std::setw(2) << iloc << "] " << nuray[iloc] << "\n";
  }

  s << decay << "\n";

  size_t nanc = ancestor.size();
  bool printanc = true; // ( opt.find("a") != std::string::npos);
  s << nanc << " ancestors " << (printanc?"":"(suppressed)") << "\n";
  if ( ! printanc ) nanc = 0;
  for ( size_t ianc = 0; ianc < nanc; ++ianc ) {
    s << "[" << std::setw(2) << ianc << "] " << ancestor[ianc] << "\n";
  }

  s << "ppv{xyz}={" << ppvx << "," << ppvy << "," << ppvz << "}\n";
  s << tgtexit << "\n";

  size_t ntraj = traj.size();
  bool printtraj = true; // ( opt.find("t") != std::string::npos);
  s << ntraj << " trajectories " << (printtraj?"":"(suppressed)") << "\n";
  if ( ! printtraj ) ntraj = 0;
  for ( size_t itraj = 0; itraj < ntraj; ++itraj ) {
    s << "[" << std::setw(2) << itraj << "] " << traj[itraj] << "\n";
  }

  size_t ni = vint.size();
  if ( ni > 0 ) {
    s << ni << " int : ";
    for ( size_t i = 0; i < ni; ++i ) { s << " " << vint[i]; }
    s << "\n";
  }
  size_t nd = vdbl.size();
  if ( nd > 0 ) {
    s << nd << " dbl : ";
    for ( size_t i = 0; i < nd; ++i ) { s << " " << vdbl[i]; }
    s << "\n";
  }
  s << "flagbits: 0x" << std::hex << std::setfill('0') 
    << std::setw(8) << flagbits << std::dec << std::setfill(' ');

  return s.str();
}
void bsim::Dk2Nu::Print(Option_t* option) const
{
  std::cout << AsString(option) << std::endl;
}
size_t bsim::Dk2Nu::indxnu() const { return ancestor.size()-1; }
size_t bsim::Dk2Nu::indxp() const { return ancestor.size()-2; }
size_t bsim::Dk2Nu::indxgp() const { return ancestor.size()-3; }
bool   bsim::Dk2Nu::overflow() const { return (flagbits&bsim::kFlgOverflow); }

std::ostream& operator<<(std::ostream& os, const bsim::Dk2Nu& dk2nu)
{
  os << dk2nu.AsString(); // << std::endl;
  return os;
}
