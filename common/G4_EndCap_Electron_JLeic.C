#ifndef MACRO_G4ENDCAPELECTRON
#define MACRO_G4ENDCAPELECTRON

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool ENDCAP_ELECTRON = false;
}

namespace G4ENDCAPELECTRON
{
  double outer_radius = 244.;  // cm
  double size_z = 60.;
  double place_z = -400. / 2. - size_z / 2.;
}  // namespace G4ENDCAPELECTRON

void EndCap_ElectronInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4ENDCAPELECTRON::outer_radius);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4ENDCAPELECTRON::place_z - G4ENDCAPELECTRON::size_z / 2.);
}

void EndCap_Electron(PHG4Reco* g4Reco)
{
  // here is our silicon:
  double inner_radius = 20.;  // cm
  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("ECELECTRON", 0);
  cyl->set_double_param("radius", inner_radius);
  cyl->set_string_param("material", "G4_Fe");
  cyl->set_double_param("thickness", G4ENDCAPELECTRON::outer_radius - inner_radius);
  cyl->set_double_param("length", G4ENDCAPELECTRON::size_z);
  cyl->set_double_param("place_z", G4ENDCAPELECTRON::place_z);
  cyl->SetActive();
  cyl->SuperDetector("ECELECTRON");
  g4Reco->registerSubsystem(cyl);

  return;
}

#endif  // MACRO_G4ENDCAPELECTRON
