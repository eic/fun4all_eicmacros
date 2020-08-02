#ifndef MACRO_G4DIRCJLEIC
#define MACRO_G4DIRCJLEIC

#include "GlobalVariables.C"

#include <g4jleic/G4JLeicDIRCSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4jleic.so)

namespace Enable
{
  bool DIRC = false;
}

namespace G4DIRC
{
  double outer_radius = 88.;
  double length = 340.;
}  // namespace G4DIRC

void JLDIRCInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4DIRC::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4DIRC::length / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4DIRC::length / 2.);
}

double JLDIRC(PHG4Reco* g4Reco,
              double radius)
{
  if (radius > 81)
  {
    cout << "radius " << radius << " too large (>81) to fit DIRC" << endl;
    exit(1);
  }
  G4JLeicDIRCSubsystem* jldirc = new G4JLeicDIRCSubsystem();
  jldirc->SetActive();
  jldirc->SuperDetector("JLDIRC");
  g4Reco->registerSubsystem(jldirc);
  radius = G4DIRC::outer_radius;
  return radius;
}
#endif  // MACRO_G4DIRCJLEIC
