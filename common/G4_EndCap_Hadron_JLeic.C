#ifndef MACRO_G4ENDCAPHADRON
#define MACRO_G4ENDCAPHADRON

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libg4detectors.so)

namespace Enable
{
  bool ENDCAP_HADRON = false;
}

namespace G4ENDCAPHADRON
{
  double outer_radius = 243.;  // cm
  double length = 250.;
  double place_z = 400. / 2. + length / 2. + 170.;
}  // namespace G4ENDCAPHADRON

void EndCap_HadronInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4ENDCAPHADRON::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4ENDCAPHADRON::place_z + G4ENDCAPHADRON::length / 2.);
}

void EndCap_Hadron(PHG4Reco *g4Reco)
{
  double hadron_inner_radius = 80.;  // cm
  PHG4CylinderSubsystem *mothervol = new PHG4CylinderSubsystem("EndCapHadronContainer", 0);
  mothervol->set_color(0.3, 0, 3., 0.1);
  mothervol->set_double_param("radius", 20.);
  mothervol->set_string_param("material", "G4_AIR");
  mothervol->set_double_param("thickness", G4ENDCAPHADRON::outer_radius - 20.);
  mothervol->set_double_param("length", G4ENDCAPHADRON::length);
  mothervol->set_double_param("place_z", G4ENDCAPHADRON::place_z);
  g4Reco->registerSubsystem(mothervol);
  double size_z = 2.;
  double gap = 2.;
  double z_start = size_z / 2. - 250. / 2.;
  int nlayer = 25;
  for (int i = 0; i < nlayer; i++)
  {
    PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("ECHADRON", i);
    cyl->SetMotherSubsystem(mothervol);
    cyl->set_color(0.6, 0, 0.6, 1);
    cyl->set_double_param("radius", hadron_inner_radius);
    cyl->set_string_param("material", "G4_Fe");
    cyl->set_double_param("thickness", G4ENDCAPHADRON::outer_radius - hadron_inner_radius);
    cyl->set_double_param("length", size_z);
    cyl->set_double_param("place_z", z_start);
    cyl->SetActive();
    cyl->SuperDetector("ECHADRON");
    z_start += gap + size_z;
    g4Reco->registerSubsystem(cyl);
  }
  return;
}
#endif  // MACRO_G4ENDCAPHADRON
