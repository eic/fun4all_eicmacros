#ifndef MACRO_G4SETUPBEAST_C
#define MACRO_G4SETUPBEAST_C

#include <GlobalVariables.C>

#include <G4_BlackHole.C>
#include <G4_Magnet_Beast.C>
#include <G4_Pipe_EIC.C>
#include <G4_User.C>
#include <G4_World.C>

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

  //----------------------------------------
  // MAGNET

  if (Enable::MAGNET) MagnetInit();
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
    g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::kFieldBeast);
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
  // MAGNET

  if (Enable::MAGNET) radius = Magnet(g4Reco, radius);

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

#endif  // MACRO_G4SETUPBEAST_C
