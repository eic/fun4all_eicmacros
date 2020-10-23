#ifndef MACRO_G4GEMJLEIC
#define MACRO_G4GEMJLEIC

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool GEM = false;
}

namespace G4GEM
{
  double outer_radius = 115.;
}

void GemInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4GEM::outer_radius);
}

void Gem(PHG4Reco* g4Reco)
{
  // here is our silicon:
  double gem_inner_radius = 0.;    // cm
  double gem_outer_radius = 115.;  // cm
  double size_z = 1.;
  PHG4CylinderSubsystem* cyl;
  for (int ilayer = 0; ilayer < 8; ilayer++)
  {
    double irad = gem_inner_radius + 1. + 0.5 * ilayer;
    double orad = G4GEM::outer_radius - 25. + 2. * ilayer;
    cyl = new PHG4CylinderSubsystem("GemHadron", ilayer);
    cyl->set_double_param("radius", irad);
    cyl->set_string_param("material", "G4_CARBON_DIOXIDE");
    cyl->set_double_param("thickness", orad - irad);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_double_param("length", size_z);
    double place_z = 340. / 2. + 5. + 3. * ilayer;
    cyl->set_double_param("place_z", place_z);
    cyl->SetActive();
    cyl->SuperDetector("GEMHADRON");
    g4Reco->registerSubsystem(cyl);
    BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, place_z + size_z);
  }
  for (int ilayer = 0; ilayer < 8; ilayer++)
  {
    double irad = gem_inner_radius + 1. + 0.5 * ilayer;
    double orad = G4GEM::outer_radius - 25. + 2. * ilayer;
    cyl = new PHG4CylinderSubsystem("GEMELECTRON", ilayer);
    cyl->set_double_param("radius", gem_inner_radius);
    cyl->set_string_param("material", "G4_CARBON_DIOXIDE");
    cyl->set_double_param("thickness", orad - irad);
    cyl->set_int_param("lengthviarapidity", 0);
    cyl->set_double_param("length", size_z);
    double place_z = -340. / 2. - 5. - 3. * ilayer;
    cyl->set_double_param("place_z", place_z);
    cyl->SetActive();
    cyl->SuperDetector("GEMELECTRON");
    g4Reco->registerSubsystem(cyl);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, place_z - size_z);
  }

  return;
}

#endif  // MACRO_G4GEMJLEIC
