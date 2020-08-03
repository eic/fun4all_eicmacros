#ifndef MACRO_G4MAGNETJLEIC_C
#define MACRO_G4MAGNETJLEIC_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool MAGNET = false;
  bool MAGNET_ABSORBER = false;
  bool MAGNET_OVERLAPCHECK = false;
  int MAGNET_VERBOSITY = 0;
}  // namespace Enable

namespace G4MAGNET
{
  double magnet_inner_radius = 130.;
  double magnet_outer_radius = 143.;
  double magnet_length = 400.;
  double magfield_rescale = 1;
  string magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/SolenoidMag3D.TABLE");

}  // namespace G4MAGNET

void MagnetInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4MAGNET::magnet_outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4MAGNET::magnet_length / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4MAGNET::magnet_length / 2.);
}

double Magnet(PHG4Reco* g4Reco,
              double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::MAGNET_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::MAGNET_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::MAGNET_VERBOSITY);

  if (radius > G4MAGNET::magnet_inner_radius)
  {
    cout << "inconsistency: radius: " << radius
         << " larger than Magnet inner radius: " << G4MAGNET::magnet_inner_radius << endl;
    gSystem->Exit(-1);
  }

// right now Cleo is a massive Al cylinder which is likely not correct
// need to know the radiation length and adjust the material for that
// or get an idea about the actual cleo magnet design
  radius = G4MAGNET::magnet_inner_radius;
  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("MAGNET", 0);
  cyl->set_double_param("radius", G4MAGNET::magnet_inner_radius);
  cyl->set_double_param("length", G4MAGNET::magnet_length);
  cyl->set_double_param("thickness", G4MAGNET::magnet_outer_radius - G4MAGNET::magnet_inner_radius);
  cyl->set_string_param("material", "Al5083");
  cyl->SuperDetector("MAGNET");
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  radius = G4MAGNET::magnet_outer_radius;  // outside of magnet

  if (verbosity > 0)
  {
    cout << "========================= G4_Magnet.C::Magnet() ===========================" << endl;
    cout << " MAGNET Material Description:" << endl;
    cout << "  inner radius = " << G4MAGNET::magnet_inner_radius << " cm" << endl;
    cout << "  outer radius = " << G4MAGNET::magnet_outer_radius << " cm" << endl;
    cout << "  length = " << G4MAGNET::magnet_length << " cm" << endl;
    cout << "===========================================================================" << endl;
  }

  radius += no_overlapp;

  return radius;
}

#endif  // MACRO_G4MAGNETJLEIC_C
