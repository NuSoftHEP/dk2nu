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

void test_read_locations(std::string locfilename = "../etc/locations.txt") 
{
  std::vector<std::string> nameloc;
  std::vector<double> xloc, yloc, zloc;

  readWeightLocations(locfilename, nameloc, xloc, yloc, zloc);

  std::cout << "Read " << nameloc.size() << " locations" << std::endl;
  for (size_t iloc = 0; iloc < nameloc.size(); ++iloc ) {
    std::cout << "{ " << setw(10) << xloc[iloc]
              << ", " << setw(10) << yloc[iloc]
              << ", " << setw(10) << zloc[iloc]
              << " } \"" << nameloc[iloc] << "\""
              << std::endl;
  }

}
