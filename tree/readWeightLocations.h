#include <string>
#include <vector>
#include <iomanip>

/// bsim namespace for beam simulation classes and functions
namespace bsim {
  class Location;
  class DkMeta;
}

namespace bsim {
  /// Read a text file that contains a header line followed by
  /// quartets of "<xpos> <ypos> <zpos> <text string>" on separate 
  /// lines.  Fill the supplied vectors.  Trim off leading/trailing 
  /// blanks and quotes (single/double) from the string.
  /// Convention has it that positions are given in (cm).
  void readWeightLocations(std::string locfilename,
                           std::vector<bsim::Location>& locations);

  /// a variant that will fill the dkmeta object
  void readWeightLocations(std::string locfilename, bsim::DkMeta* dkmeta);

  /// print the locations
  void printWeightLocations(std::vector<bsim::Location>& locations);

  /// a variant that prints the locations in the dkmeta object
  void printWeightLocations(bsim::DkMeta* dkmeta);

} // end-of-namespace "bsim"
