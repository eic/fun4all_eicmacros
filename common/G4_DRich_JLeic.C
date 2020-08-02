#ifndef MACRO_G4DRICHJLEIC
#define MACRO_G4DRICHJLEIC

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool DRICH = false;
}

namespace G4DRICH
{
  double outer_radius = 150.;
  double size_z = 170.;
  double place_z = 400. / 2. + size_z / 2.;
}  // namespace G4DRICH

void DRichInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4DRICH::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4DRICH::place_z + G4DRICH::size_z / 2.);
}

void DRich(PHG4Reco* g4Reco)
{
  double drich_inner_radius = 20.;  // cm
  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("DRICH", 0);
  cyl->set_color(1., 1., 0.2, 0.2);
  cyl->set_double_param("radius", drich_inner_radius);
  cyl->set_string_param("material", "G4_CARBON_DIOXIDE");
  cyl->set_double_param("thickness", G4DRICH::outer_radius - drich_inner_radius);
  cyl->set_double_param("length", G4DRICH::size_z);
  cyl->set_double_param("place_z", G4DRICH::place_z);
  cyl->SetActive();
  cyl->SuperDetector("DRICH");
  g4Reco->registerSubsystem(cyl);
}

#endif  // MACRO_G4DRICHJLEIC
