#include <string>
#include <vector>
#include <iomanip>

class dkmeta;

/// Read a text file that contains a header line followed by
/// quartets of "<xpos> <ypos> <zpos> <text string>" on separate 
/// lines.  Fill the supplied vectors.  Trim off leading/trailing 
/// blanks and quotes (single/double) from the string.
/// Convention has it that positions are given in (cm).
void readWeightLocations(std::string locfilename,
                         std::vector<std::string>& nameloc,
                         std::vector<double>& xloc,
                         std::vector<double>& yloc,
                         std::vector<double>& zloc);

/// a variant that will fill the dkmeta object
void readWeightLocations(std::string locfilename, dkmeta* dkmetaObj);
