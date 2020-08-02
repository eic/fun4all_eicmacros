#ifndef MACRO_G4BARRELHCALJLEIC_C
#define MACRO_G4BARRELHCALJLEIC_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool BARREL_HCAL = false;
}

namespace G4BARRELHCAL
{
  double inner_radius = 144.;
  double outer_radius = inner_radius + 100.;
  double length = 460.;
}  // namespace G4BARRELHCAL

void Barrel_HcalInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4BARRELHCAL::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4BARRELHCAL::length / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4BARRELHCAL::length / 2.);
}

double Barrel_Hcal(PHG4Reco* g4Reco,
                   double radius)
{
  // here is our silicon:
  double gap = 2.;
  double thickness = 2.;  // cm
  if (radius > G4BARRELHCAL::inner_radius)
  {
    cout << "inconsistency: radius: " << radius
         << " larger than Barrel Hcal inner radius: " << G4BARRELHCAL::inner_radius << endl;
    gSystem->Exit(-1);
  }

  PHG4CylinderSubsystem* cyl;
  radius = G4BARRELHCAL::inner_radius;
  for (int ilayer = 0; ilayer < 25; ilayer++)
  {
    if (radius > G4BARRELHCAL::outer_radius)
    {
      continue;
    }
    cyl = new PHG4CylinderSubsystem("BARRELHCAL", ilayer);
    cyl->set_double_param("radius", radius);
    cyl->set_string_param("material", "G4_Fe");
    cyl->set_double_param("thickness", thickness);
    cyl->set_double_param("length", G4BARRELHCAL::length);
    cyl->SetActive();
    cyl->SuperDetector("BARRELHCAL");
    g4Reco->registerSubsystem(cyl);
    radius += gap + thickness;
  }
  return G4BARRELHCAL::outer_radius;
}
#endif  // MACRO_G4BARRELHCALJLEIC_C
