#ifndef MACRO_DISPLAYON_C
#define MACRO_DISPLAYON_C

#include <g4main/PHG4Reco.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)

namespace Enable
{
  bool DISPLAY = false;
}

// This starts the QT based G4 gui which takes control
// when x'ed out it will return a pointer to PHG4Reco so
// the gui can be startrd again
PHG4Reco *QTGui()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco("PHG4RECO");
  g4->InitRun(se->topNode());
  g4->ApplyDisplayAction();
  g4->StartGui();
  return g4;
}

// stupid macro to turn on the geant4 display
// we ask Fun4All for a pointer to PHG4Reco
// using the ApplyCommand will start up the
// G4 cmd interpreter and graphics system
// the vis.mac contains the necessary commands to
// start up the visualization, the next event will
// be displayed. Do not execute this macro
// before PHG4Reco was registered with Fun4All
PHG4Reco * DisplayOn(const char *mac = "vis.mac")
{
  char cmd[100];
  Fun4AllServer *se = Fun4AllServer::instance();
  PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco("PHG4RECO");
  g4->InitRun(se->topNode());
  g4->ApplyDisplayAction();
  sprintf(cmd, "/control/execute %s", mac);
  g4->ApplyCommand(cmd);
  // draw by particle type and set nice color
  g4->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e+ steelblue");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e- steelblue");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set pi+ red");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set pi- red");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set kaon- mediumvioletred");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set kaon+ mediumvioletred");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set kaon0 mediumvioletred");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set proton+ orange");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set proton- orange");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set neutron lightgrey");
  // g4->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set gamma wheat");
//     g4->ApplyCommand("/vis/scene/add/trajectories smooth rich");
  // remove neutrons and neutrinos
  g4->ApplyCommand("/vis/filtering/trajectories/create/particleFilter");
  // g4->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/add neutron");
  g4->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/add neutrino");
  g4->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/invert true");

  // set background white for presentations
  g4->ApplyCommand("/vis/viewer/set/background white");

  
  return g4;
}
// print out the commands I always forget
void displaycmd()
{
  cout << "draw axis: " << endl;
  cout << " g4->ApplyCommand(\"/vis/scene/add/axes 0 0 0 50 cm\")" << endl;
  cout << "zoom" << endl;
  cout << " g4->ApplyCommand(\"/vis/viewer/zoom 1\")" << endl;
  cout << "viewpoint:" << endl;
  cout << " g4->ApplyCommand(\"/vis/viewer/set/viewpointThetaPhi 0 0\")" << endl;
  cout << "panTo:" << endl;
  cout << " g4->ApplyCommand(\"/vis/viewer/panTo 0 0 cm\")" << endl;
  cout << "print to eps:" << endl;
  cout << " g4->ApplyCommand(\"/vis/ogl/printEPS\")" << endl;
  cout << "set background color:" << endl;
  cout << " g4->ApplyCommand(\"/vis/viewer/set/background white\")" << endl;
}



#endif

