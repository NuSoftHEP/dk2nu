#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "tree/readWeightLocations.h"

#include "tree/dkmeta.h"

#include "TSystem.h"

/// Read a text file that contains a header line followed by
/// quartets of "<xpos> <ypos> <zpos> <text string>" on separate 
/// lines.  Fill the supplied vectors.  Trim off leading/trailing 
/// blanks and quotes (single/double) from the string.
/// Convention has it that positions are given in (cm).
void bsim::readWeightLocations(std::string locfilename,
                               std::vector<bsim::Location>& locations)
{

  const char* path = gSystem->ExpandPathName(locfilename.c_str());
  std::ifstream locfile(path);

  int iline=0;

  char comment_buffer[1000];

  // read lines
  char tmp[1001];
  size_t tmplen = sizeof(tmp);
  while ( ! locfile.eof() ) {
    char c;
    while ( ( c = locfile.get() ) == ' ' ) {}; // eat leading spaces
    if ( c == '#' ) {
      // comment_buffer line
      locfile.getline(comment_buffer,sizeof(comment_buffer));
      continue;
    }
    locfile.putback(c);
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
    bsim::Location alocation(x,y,z,name);
    locations.push_back(alocation);
  }

}

/// a variant that will fill the dkmeta object
void bsim::readWeightLocations(std::string locfilename, bsim::DkMeta* dkmeta)
{
  ///  read & print the locations where weights are to be calculated

  std::vector<bsim::Location>& locations = dkmeta->location;
  /// make an entry for the random decay
  bsim::Location rndmloc(0,0,0,"random decay");
  locations.push_back(rndmloc);

  /// read and parse the text file for additional positions
  /// use the vector version
  readWeightLocations(locfilename, locations);
}

void bsim::printWeightLocations(std::vector<bsim::Location>& locations)
{
  size_t nl = locations.size();
  std::cout << nl << " locations:\n";
  for ( size_t l = 0; l < nl; ++l ) {
    std::cout << " [" << std::setw(2) << l << "] " << locations[l] << "\n";
  }
}

void bsim::printWeightLocations(bsim::DkMeta* dkmeta)
{
  bsim::printWeightLocations(dkmeta->location);
}
