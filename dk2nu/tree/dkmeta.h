/**
 * \class dkmeta
 * \file  dkmeta.h
 *
 * \brief A class that defines the "dkmeta" object used as the 
 *        branch for a TTree for the output of meta-data from 
 *        neutrino flux simulations such as g4numi, g4numi_flugg, etc.
 *        This tree has one entry of this type for the file.  Kept
 *        as a tree so files can be chained.
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2012-11-07 01:35:47 $
 *
 * Contact: rhatcher@fnal.gov
 *
 * $Id: dkmeta.h,v 1.1 2012-11-07 01:35:47 rhatcher Exp $
 *
 * Notes tagged with "DKMETA" are questions that should be answered
 */

#ifndef DKMETA_H
#define DKMETA_H

#include "TROOT.h"
#include "TObject.h"

#include <vector>
#include <string>

class dkmeta;
std::ostream& operator<<(std::ostream& os, const dkmeta& entry);

class dkmeta
{
private:
  ClassDef(dkmeta,3) // KEEP THIS UP-TO-DATE!  increment for each change

public:
  /**
   *   Public methods for constructing/destruction and resetting the data
   */
  dkmeta();
  virtual     ~dkmeta();
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
   *  General Info:
   */
   Int_t    job;           ///< identifying job # (keep files distinct)
   Double_t pots;          ///< protons-on-target

   /**
    * DKMETA:  
    * formatted strings are most flexible ...
    * but not necessarily convenient to use
    * ??? Should parts of these be standardized ??? 
    */
   std::string beamsim;    ///< e.g. "flugg" or "g4numi/<tag>"
   std::string physics;    ///< e.g. "fluka08", "g4.9.3p01"
   std::string physcuts;   ///< tracking cuts    e.g. "threshold=0.1GeV"
   std::string tgtcfg;     ///< target config    e.g. "minos/epoch3/-10cm"
   std::string horncfg;    ///< horn config      e.g. "FHC/185A/LE/h1xoff=1mm"
   std::string dkvolcfg;   ///< decay vol config e.g. "helium" or "vacuum"

  /**
   *=======================================================================
   *  Beam Info:
   */
   Double_t beam0x;       ///< x of beam center at start
   Double_t beam0y;       ///< y of beam center at start
   Double_t beam0z;       ///< z of beam start
   Double_t beamhwidth;   ///< horizontal width of beam
   Double_t beamvwidth;   ///< vertical width of beam
   Double_t beamdxdz;     ///< beam slope dx/dz
   Double_t beamdydz;     ///< beam slope dy/dz

  /**
   *=======================================================================
   *  Detector Position Info:
   *  Values are in beam coordinate system w/ units of "cm"
   */
   std::vector<Double_t> xloc;   ///< x positions of detectors
   std::vector<Double_t> yloc;   ///< y positions of detectors
   std::vector<Double_t> zloc;   ///< z positions of detectors

   std::vector<std::string> nameloc; ///< names of detector locations (e.g. "NOvA-ND-3x3")

  /**
   *=======================================================================
   *  Special Info:
   *  Document extensibility enhancements 
   */
   std::vector<std::string> vintnames;    ///< names of elements for user defined vector of integers
   std::vector<std::string> vdblnames;    ///< names of elements for user defined vector of doubles
   
 };

#endif
