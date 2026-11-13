//
// test reading locations text file
// this script can be run using:
//       root -b -q load_dk2nu.C test_read_locations.C+
// 
// rhatcher@fnal.gov  2012-10-02
//====================================================================

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

// dk2nu includes (run load_dk2nu.C first to set the include path)
#include "tree/dkmeta.h"
#include "tree/readWeightLocations.h"

void test_read_locations(std::string locfilename = "${DK2NU}/etc/locations.txt") 
{
  std::vector<bsim::Location> location;
  bsim::readWeightLocations(locfilename, location);

  size_t nl = location.size();
  std::cout << nl << " locations:\n";
  for ( size_t l = 0; l < nl; ++l ) {
    std::cout << " [" << std::setw(2) << l << "] " << location[l] << "\n";
  }
  std::cout << endl;

}
