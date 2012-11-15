#include <iostream>
#include <iomanip>

#include "dk2nu/tree/calcLocationWeights.h"

#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"

/// user interface
void bsim::calcLocationWeights(const bsim::DkMeta* dkmeta, bsim::Dk2Nu* dk2nu)
{
  size_t nloc = dkmeta->location.size();
  for (size_t iloc = 0; iloc < nloc; ++iloc ) {
    // skip calculation for random location ... should already be filled
    const std::string rkey = "random decay";
    if ( dkmeta->location[iloc].name == rkey ) {
      if ( iloc != 0 ) {
        std::cerr << "calcLocationWeights \"" << rkey << "\""
                  << " isn't the 0-th entry" << std::endl;
        assert(0);
      }
      if ( dk2nu->nuray.size() != 1 ) {
        std::cerr << "calcLocationWeights \"" << rkey << "\""
                  << " nuenergy[" << iloc << "] not filled" << std::endl;
        assert(0);
      }
      continue;
    }
    TVector3 xyzDet(dkmeta->location[iloc].x,
                    dkmeta->location[iloc].y,
                    dkmeta->location[iloc].z);  // position to evaluate
    double enu_xy = 0;  // give a default value
    double wgt_xy = 0;  // give a default value
    int status = bsim::calcEnuWgt(dk2nu,xyzDet,enu_xy,wgt_xy);
    if ( status != 0 ) {
      std::cerr << "bsim::calcEnuWgt returned " << status << " for " 
                << dkmeta->location[iloc].name << std::endl;
    }
    // with the recalculated energy compute the momentum components
    TVector3 xyzDk(dk2nu->decay.vx,dk2nu->decay.vy,dk2nu->decay.vz);  // origin of decay
    TVector3 p3 = enu_xy * (xyzDet - xyzDk).Unit();
    bsim::NuRay anuray(p3.x(), p3.y(), p3.z(), enu_xy, wgt_xy);
    dk2nu->nuray.push_back(anuray);
  }
}

//___________________________________________________________________________
int bsim::calcEnuWgt(const bsim::Decay& decay, const TVector3& xyz,
                     double& enu, double& wgt_xy)
{
  // Neutrino Energy and Weight at arbitrary point
  // Based on:
  //   NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
  //   Title:   Neutrino Beam Simulation using PAW with Weighted Monte Carlos
  //   Author:  Rick Milburn
  //   Date:    1995-10-01

  // History:
  // jzh  3/21/96 grab R.H.Milburn's weighing routine
  // jzh  5/ 9/96 substantially modify the weighting function use dot product 
  //              instead of rotation vecs to get theta get all info except 
  //              det from ADAMO banks neutrino parent is in Particle.inc
  //              Add weighting factor for polarized muon decay
  // jzh  4/17/97 convert more code to double precision because of problems 
  //              with Enu>30 GeV
  // rwh 10/ 9/08 transliterate function from f77 to C++

  // Original function description:
  //   Real function for use with PAW Ntuple To transform from destination
  //   detector geometry to the unit sphere moving with decaying hadron with
  //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decays will
  //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
  //   fraction of decays that hit the target.  For a given target point and 
  //   area, and given x-y components of decay transverse location and slope,
  //   and given decay distance from target ans given decay GAMMA and 
  //   rest-frame neutrino energy, the lab energy at the target and the 
  //   fractional solid angle in the rest-frame are determined.
  //   For muon decays, correction for non-isotropic nature of decay is done.

  // Arguments:
  //    dk2nu    :: contains current decay information
  //    xyz      :: 3-vector of position to evaluate
  //                in *beam* frame coordinates  (cm units)
  //    enu      :: resulting energy
  //    wgt_xy   :: resulting weight
  // Return:
  //    (int)    :: error code
  // Assumptions:
  //    Energies given in GeV
  //    Particle codes have been translated from GEANT into PDG codes

  // for now ... these masses _should_ come from TDatabasePDG 
  // but use these hard-coded values to "exactly" reproduce old code
  //
  const double kPIMASS = 0.13957;
  const double kKMASS  = 0.49368;
  const double kK0MASS = 0.49767;
  const double kMUMASS = 0.105658389;
  const double kOMEGAMASS = 1.67245;

  const int kpdg_nue       =   12;  // extended Geant 53
  const int kpdg_nuebar    =  -12;  // extended Geant 52
  const int kpdg_numu      =   14;  // extended Geant 56
  const int kpdg_numubar   =  -14;  // extended Geant 55

  const int kpdg_muplus     =   -13;  // Geant  5
  const int kpdg_muminus    =    13;  // Geant  6
  const int kpdg_pionplus   =   211;  // Geant  8
  const int kpdg_pionminus  =  -211;  // Geant  9
  const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short    =   310;  // Geant 16
  const int kpdg_k0mix      =   311;  
  const int kpdg_kaonplus   =   321;  // Geant 11
  const int kpdg_kaonminus  =  -321;  // Geant 12
  const int kpdg_omegaminus =  3334;  // Geant 24
  const int kpdg_omegaplus  = -3334;  // Geant 32

  const double kRDET = 100.0;   // set to flux per 100 cm radius

  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();

  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error


  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = kPIMASS;
  switch ( decay.ptype ) {
  case kpdg_pionplus:
  case kpdg_pionminus:
    parent_mass = kPIMASS;
    break;
  case kpdg_kaonplus:
  case kpdg_kaonminus:
    parent_mass = kKMASS;
    break;
  case kpdg_k0long:
  case kpdg_k0short:
  case kpdg_k0mix:
    parent_mass = kK0MASS;
    break;
  case kpdg_muplus:
  case kpdg_muminus:
    parent_mass = kMUMASS;
    break;
  case kpdg_omegaminus:
  case kpdg_omegaplus:
    parent_mass = kOMEGAMASS;
    break;
  default:
    std::cerr << "bsim::calcEnuWgt unknown particle type " << decay.ptype
              << std::endl << std::flush;
    assert(0);
    return 1;
  }

  double parentp2 = ( decay.pdpx*decay.pdpx +
                      decay.pdpy*decay.pdpy +
                      decay.pdpz*decay.pdpz );
  double parent_energy = TMath::Sqrt( parentp2 +
                                     parent_mass*parent_mass);
  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;
  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = decay.necm;
  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos-decay.vx)*(xpos-decay.vx) +
                            (ypos-decay.vy)*(ypos-decay.vy) +
                            (zpos-decay.vz)*(zpos-decay.vz) );

  double emrat = 1.0;
  double costh_pardet = -999., theta_pardet = -999.;

  // boost correction, but only if parent hasn't stopped
  if ( parentp > 0. ) {
    costh_pardet = ( decay.pdpx*(xpos-decay.vx) +
                     decay.pdpy*(ypos-decay.vy) +
                     decay.pdpz*(zpos-decay.vz) ) 
                     / ( parentp * rad);
    if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
    if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
    theta_pardet = TMath::ACos(costh_pardet);

    // Weighted neutrino energy in beam, approx, good for small theta
    emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  }

  enu = emrat * enuzr;  // the energy ... normally

  // Get solid angle/4pi for detector element
  double sangdet = ( kRDET*kRDET / 
                     ( (zpos-decay.vz)*(zpos-decay.vz) ) ) / 4.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

  // Done for all except polarized muon decay
  // in which case need to modify weight 
  // (must be done in double precision)
  if ( decay.ptype == kpdg_muplus || decay.ptype == kpdg_muminus) {
    double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

    // Boost neu neutrino to mu decay CM
    beta[0] = decay.pdpx / parent_energy;
    beta[1] = decay.pdpy / parent_energy;
    beta[2] = decay.pdpz / parent_energy;
    p_nu[0] = (xpos-decay.vx)*enu/rad;
    p_nu[1] = (ypos-decay.vy)*enu/rad;
    p_nu[2] = (zpos-decay.vz)*enu/rad;
    partial = gamma * 
      (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
    partial = enu - partial/(gamma+1.0);
    // the following calculation is numerically imprecise
    // especially p_dcm_nu[2] leads to taking the difference of numbers 
    //  of order ~10's and getting results of order ~0.02's
    // for g3numi we're starting with floats (ie. good to ~1 part in 10^7)
    p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
    p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
    p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
    p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
                               p_dcm_nu[1]*p_dcm_nu[1] +
                               p_dcm_nu[2]*p_dcm_nu[2] );

    // Boost parent of mu to mu production CM
    double particle_energy = decay.ppenergy;
    gamma = particle_energy/parent_mass;
    beta[0] = decay.ppdxdz * decay.pppz / particle_energy;
    beta[1] = decay.ppdydz * decay.pppz / particle_energy;
    beta[2] =                    decay.pppz / particle_energy;
    partial = gamma * ( beta[0]*decay.muparpx + 
                        beta[1]*decay.muparpy + 
                        beta[2]*decay.muparpz );
    partial = decay.mupare - partial/(gamma+1.0);
    p_pcm_mp[0] = decay.muparpx - beta[0]*gamma*partial;
    p_pcm_mp[1] = decay.muparpy - beta[1]*gamma*partial;
    p_pcm_mp[2] = decay.muparpz - beta[2]*gamma*partial;
    double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
                                 p_pcm_mp[1]*p_pcm_mp[1] +
                                 p_pcm_mp[2]*p_pcm_mp[2] );

    const double eps = 1.0e-30;  // ? what value to use
    if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
      return 3; // mu missing parent info?
    }
    // Calc new decay angle w.r.t. (anti)spin direction
    double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
                     p_dcm_nu[1]*p_pcm_mp[1] +
                     p_dcm_nu[2]*p_pcm_mp[2] ) /
                   ( p_dcm_nu[3]*p_pcm );
    if ( costh >  1.0 ) costh =  1.0;
    if ( costh < -1.0 ) costh = -1.0;
    // Calc relative weight due to angle difference
    double wgt_ratio = 0.0;
    switch ( decay.ntype ) {
    case kpdg_nue:
    case kpdg_nuebar:
      wgt_ratio = 1.0 - costh;
      break;
    case kpdg_numu:
    case kpdg_numubar:
    {
      double xnu = 2.0 * enuzr / kMUMASS;
      wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
      break;
    }
    default:
      return 2; // bad neutrino type
    }
    wgt_xy = wgt_xy * wgt_ratio;

  } // ptype is muon

  return 0;
}
//___________________________________________________________________________

int bsim::calcEnuWgt(const bsim::Dk2Nu* dk2nu, const TVector3& xyz,
                     double& enu, double& wgt_xy)
{
  return bsim::calcEnuWgt(dk2nu->decay,xyz,enu,wgt_xy);
}
//___________________________________________________________________________
