#ifndef MACRO_G4MAGNETBEAST_C
#define MACRO_G4MAGNETBEAST_C

#include "GlobalVariables.C"

#include <eicdetectors/BeastMagnetSubsystem.h>

#include <g4main/PHG4Reco.h>

R__LOAD_LIBRARY(libeicdetectors.so)

namespace Enable
{
  bool MAGNET = false;
  bool MAGNET_ABSORBER = false;
}  // namespace Enable

namespace G4MAGNET
{
  double magnet_outer_radius = 300.;
  double magnet_length = 500.;
  double magfield_rescale = 1;
  string magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat");

}  // namespace G4MAGNET

void MagnetInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4MAGNET::magnet_outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4MAGNET::magnet_length / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4MAGNET::magnet_length / 2.);
}

void Magnet(PHG4Reco* g4Reco)

{
  bool AbsorberActive = Enable::ABSORBER || Enable::MAGNET_ABSORBER;

    BeastMagnetSubsystem *beast = new BeastMagnetSubsystem();
    beast->set_string_param("GDMPath",(string(getenv("CALIBRATIONROOT")) + string("/Magnet/BeastSolenoid.gdml")));
    beast->set_string_param("TopVolName","SOLENOID");
    beast->SetActive(AbsorberActive);
    beast->SuperDetector("MAGNET");
    g4Reco->registerSubsystem(beast);
  return;
}

#endif  // MACRO_G4MAGNETBEAST_C
