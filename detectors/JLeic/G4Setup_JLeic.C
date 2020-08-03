#ifndef MACRO_G4SETUPJLEIC_C
#define MACRO_G4SETUPJLEIC_C

#include "GlobalVariables.C"

#include "G4_Barrel_Hcal_JLeic.C"
#include "G4_BeamLine_JLeic.C"
#include "G4_BlackHole.C"
#include "G4_CTD_JLeic.C"
#include "G4_DIRC_JLeic.C"
#include "G4_DRich_JLeic.C"
#include "G4_EndCap_Electron_JLeic.C"
#include "G4_EndCap_Hadron_JLeic.C"
#include "G4_Gem_JLeic.C"
#include "G4_Magnet_Cleo.C"
#include "G4_Pipe_EIC.C"
#include "G4_User.C"
#include "G4_VTX_JLeic.C"
#include "G4_World.C"

#include <g4decayer/EDecayType.hh>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <phfield/PHFieldConfig.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libg4decayer.so)

void G4Init()
{
  // load detector/material macros and execute Init() function

  if (Enable::PIPE) PipeInit();

  if (Enable::VTX) VTXInit();

  if (Enable::CTD) CTDInit();

  //----------------------------------------
  // MAGNET

  if (Enable::MAGNET) MagnetInit();

  if (Enable::GEM) GemInit();

  if (Enable::DIRC) JLDIRCInit();

  if (Enable::BARREL_HCAL) Barrel_HcalInit();

  if (Enable::DRICH) DRichInit();

  if (Enable::ENDCAP_ELECTRON) EndCap_ElectronInit();

  if (Enable::ENDCAP_HADRON) EndCap_HadronInit();

  if (Enable::BEAMLINE) BeamLineInit();
}

void G4Setup()
{
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4Reco *g4Reco = new PHG4Reco();

  WorldInit(g4Reco);

  // global coverage used for length of cylinders if lengthviarapidity is set
  // probably needs to be adjusted for JLeic
  g4Reco->set_rapidity_coverage(1.1);

  if (G4P6DECAYER::decayType != EDecayType::kAll)
  {
    g4Reco->set_force_decay(G4P6DECAYER::decayType);
  }

  double fieldstrength;
  istringstream stringline(G4MAGNET::magfield);
  stringline >> fieldstrength;
  if (stringline.fail())
  {  // conversion to double fails -> we have a string

    g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::kFieldCleo);
  }
  else
  {
    g4Reco->set_field(fieldstrength);  // use const soleniodal field
  }
  g4Reco->set_field_rescale(G4MAGNET::magfield_rescale);

  double radius = 0.;

  //----------------------------------------
  // PIPE
  if (Enable::PIPE) radius = Pipe(g4Reco, radius);

  //----------------------------------------
  // VTX
  if (Enable::VTX) radius = VTX(g4Reco, radius);

  //----------------------------------------
  // CTD
  if (Enable::CTD) radius = CTD(g4Reco, radius);

  //----------------------------------------
  // DIRC
  if (Enable::DIRC) radius = JLDIRC(g4Reco, radius);

  //----------------------------------------
  // MAGNET

  if (Enable::MAGNET) radius = Magnet(g4Reco, radius);

  //----------------------------------------
  // Barrel Hcal

  if (Enable::BARREL_HCAL) radius = Barrel_Hcal(g4Reco, radius);

  //----------------------------------------
  // Gem (hadron and electron going side

  if (Enable::GEM) Gem(g4Reco);

  if (Enable::DRICH) DRich(g4Reco);

  if (Enable::ENDCAP_ELECTRON) EndCap_Electron(g4Reco);

  if (Enable::ENDCAP_HADRON) EndCap_Hadron(g4Reco);

  if (Enable::BEAMLINE) BeamLine(g4Reco);

  if (Enable::USER) UserDetector(g4Reco);

  //----------------------------------------
  // BLACKHOLE

  if (Enable::BLACKHOLE) BlackHole(g4Reco, radius);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);
  // finally adjust the world size in case the default is too small
  WorldSize(g4Reco, radius);

  se->registerSubsystem(g4Reco);
  return;
}

void DstCompress(Fun4AllDstOutputManager *out)
{
  if (out)
  {
    out->StripNode("G4HIT_PIPE");
  }
}

#endif  // MACRO_G4SETUPJLEIC_C
