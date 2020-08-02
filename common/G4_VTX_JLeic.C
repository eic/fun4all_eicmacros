#ifndef MACRO_G4VTXJLEIC_C
#define MACRO_G4VTXJLEIC_C

#include "GlobalVariables.C"

#include <g4jleic/G4JLeicVTXSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4jleic.so)

namespace Enable
{
  bool VTX = false;
}

namespace G4VTX
{
  double outer_radius = 16.;
  double length = 48;
}  // namespace G4VTX

void VTXInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4VTX::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4VTX::length / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4VTX::length / 2.);
}

double VTX(PHG4Reco* g4Reco,
           double radius)
{
  if (radius > 3.5)
  {
    cout << "radius too large to fit VTX" << endl;
    exit(1);
  }
  G4JLeicVTXSubsystem* jlvtx = new G4JLeicVTXSubsystem();
  jlvtx->SetActive();
  jlvtx->SuperDetector("JLVTX");
  g4Reco->registerSubsystem(jlvtx);
  radius = G4VTX::outer_radius;
  return radius;
}

#endif  // MACRO_G4VTXJLEIC_C
