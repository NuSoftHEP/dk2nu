#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace bsim;

#pragma link C++ class bsim::NuRay+;
#pragma link C++ class bsim::Decay+;
#pragma link C++ class bsim::Ancestor+;
#pragma link C++ class bsim::TgtExit+;
#pragma link C++ class bsim::Traj+;
#pragma link C++ class bsim::Dk2Nu+;

#pragma link C++ class std::vector<bsim::NuRay>+;
#pragma link C++ class std::vector<bsim::Ancestor>+;
#pragma link C++ class std::vector<bsim::Traj>+;

#pragma link C++ class bsim::Location+;
#pragma link C++ class bsim::DkMeta+;

#pragma link C++ class std::vector<bsim::Location>+;

#pragma link C++ class bsim::NuChoice+;

#pragma link C++ function bsim::readWeightLocations;
#pragma link C++ function bsim::printWeightLocations;
#pragma link C++ function bsim::calcLocationWeights;
#pragma link C++ function bsim::calcEnuWgt;

#pragma link C++ function bsim::IsDefault;

#pragma link C++ function operator<<(std::ostream&, const bsim::NuRay&);
#pragma link C++ function operator<<(std::ostream&, const bsim::Decay&);
#pragma link C++ function operator<<(std::ostream&, const bsim::Ancestor&);
#pragma link C++ function operator<<(std::ostream&, const bsim::TgtExit&);
#pragma link C++ function operator<<(std::ostream&, const bsim::Traj&);
#pragma link C++ function operator<<(std::ostream&, const bsim::Dk2Nu&);

#pragma link C++ function operator<<(std::ostream&, const bsim::Location&);
#pragma link C++ function operator<<(std::ostream&, const bsim::DkMeta&);

#pragma link C++ function operator<<(std::ostream&, const bsim::NuChoice&);
#endif

