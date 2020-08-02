#ifndef MACRO_G4CTD_JLEIC
#define MACRO_G4CTD_JLEIC

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool CTD = false;
}

namespace G4CTD
{
  double outer_radius = 80.;
  double size_z = 340.;
}  // namespace G4CTD

void CTDInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4CTD::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4CTD::size_z / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4CTD::size_z / 2.);
}

double CTD(PHG4Reco* g4Reco,
           double radius)
{
  // here is our silicon:
  double shift_z = 0;         // shift z from GlobalVariables.C
  double inner_radius = 21.;  // cm
  double si_layer_gap = 5.;
  double si_thickness = 0.001;  // cm
  if (radius > inner_radius)
  {
    cout << "inconsistency: radius: " << radius
         << " larger than CTD inner radius: " << inner_radius << endl;
    gSystem->Exit(-1);
  }

  PHG4CylinderSubsystem* cyl;
  for (int ilayer = 0; ilayer < 15; ilayer++)
  {
    radius = inner_radius + ilayer * si_layer_gap;
    if (radius + si_thickness > G4CTD::outer_radius)
    {
      continue;
    }
    cyl = new PHG4CylinderSubsystem("JLCTD", ilayer);
    cyl->set_color(0.1, 0, 1., 0.1);
    cyl->set_double_param("radius", radius);
    cyl->set_string_param("material", "G4_Si");
    cyl->set_double_param("thickness", si_thickness);
    cyl->set_double_param("place_z", shift_z);
    cyl->SetActive();
    cyl->SuperDetector("JLCTD");
    cyl->set_double_param("length", G4CTD::size_z);
    g4Reco->registerSubsystem(cyl);
  }
  radius += si_thickness + no_overlapp;
  return G4CTD::outer_radius;
}

#endif  // MACRO_G4CTD_JLEIC
