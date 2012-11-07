#include <iostream>
#include <fstream>
#include <iomanip>

#include "dk2nu/tree/readWeightLocations.h"

#include "dk2nu/tree/dkmeta.h"

/// Read a text file that contains a header line followed by
/// quartets of "<xpos> <ypos> <zpos> <text string>" on separate 
/// lines.  Fill the supplied vectors.  Trim off leading/trailing 
/// blanks and quotes (single/double) from the string.
/// Convention has it that positions are given in (cm).
void readWeightLocations(std::string locfilename,
                         std::vector<std::string>& nameloc,
                         std::vector<double>& xloc,
                         std::vector<double>& yloc,
                         std::vector<double>& zloc)
{

  std::ifstream locfile(locfilename.c_str());

  int iline=0;

  // read/skip header line in text file
  char header[1000];
  locfile.getline(header,sizeof(header));
  ++iline;

  // read lines
  char tmp[1001];
  size_t tmplen = sizeof(tmp);
  while ( ! locfile.eof() ) {
    double x, y, z;
    locfile >> x >> y >> z;
    locfile.getline(tmp,tmplen-1);
    size_t i = locfile.gcount();
    // make sure the c-string is null terminated
    size_t inull = i;
    //if ( inull < 0 ) inull = 0;
    if ( inull > tmplen-1 ) inull = tmplen-1;
    tmp[inull] = '\0';
    std::string name(tmp);
    // ignore leading & trailing blanks (and any single/double quotes)
    size_t ilast  = name.find_last_not_of(" \t\n'\"");
    name.erase(ilast+1,std::string::npos);  // trim tail
    size_t ifirst = name.find_first_not_of(" \t\n'\"");
    name.erase(0,ifirst);  // trim head

    ++iline;
    if ( ! locfile.good() ) {
      //if ( verbose)
      //  std::cout << "stopped reading on line " << iline << std::endl;
      break;
    }
    nameloc.push_back(name);
    xloc.push_back(x);
    yloc.push_back(y);
    zloc.push_back(z);
  }

}

/// a variant that will fill the dkmeta object
void readWeightLocations(std::string locfilename, dkmeta* dkmetaObj)
{
  ///  read & print the locations where weights are to be calculated
  std::vector<std::string>& nameloc = dkmetaObj->nameloc;
  std::vector<double>& xloc         = dkmetaObj->xloc;
  std::vector<double>& yloc         = dkmetaObj->yloc;
  std::vector<double>& zloc         = dkmetaObj->zloc;

  /// make an entry for the random decay
  nameloc.push_back("random decay");
  xloc.push_back(0);  // positions for random case are bogus
  yloc.push_back(0);
  zloc.push_back(0);

  /// read and parse the text file for additional positions
  /// use the vector version
  readWeightLocations(locfilename, nameloc, xloc, yloc, zloc);
}
