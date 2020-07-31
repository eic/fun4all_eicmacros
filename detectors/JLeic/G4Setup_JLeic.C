#ifndef MACRO_G4SETUPJLEIC_C
#define MACRO_G4SETUPJLEIC_C

#include "GlobalVariables.C"

#include "G4_BlackHole.C"
#include "G4_Magnet_JLeic.C"
#include "G4_Pipe_EIC.C"
#include "G4_User.C"
#include "G4_World.C"
#include "G4_VTX_JLeic.C"
#include "G4_CTD_JLeic.C"

#include "G4_Gem.C"
#include "G4_JLDIRC.C"
#include "G4_Barrel_Hcal.C"
#include "G4_BeamLine.C"
#include "G4_DRich.C"
#include "G4_EndCap_Electron.C"
#include "G4_EndCap_Hadron.C"

#include <g4decayer/EDecayType.hh>

#include <g4eval/PHG4DstCompressReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <phfield/PHFieldConfig.h>
#include <g4main/HepMCNodeReader.h>
class SubsysReco;
R__LOAD_LIBRARY(libg4decayer.so)
R__LOAD_LIBRARY(libg4detectors.so)

void G4Init(const bool do_gem = true,
            const bool do_jldirc = true,
            const bool do_barrel_hcal = true,
            const bool do_drich = true,
            const bool do_endcap_electron = true,
            const bool do_endcap_hadron = true,
            const bool do_beamline = true
	    )
  {

  // load detector/material macros and execute Init() function

  if (Enable::PIPE) PipeInit();

  if (Enable::VTX) VTXInit();

  if (Enable::CTD) CTDInit();

  //----------------------------------------
  // MAGNET

  if (Enable::MAGNET) MagnetInit();

  if (do_gem)
    {
      gROOT->LoadMacro("G4_Gem.C");
      GemInit();
    }
  if (do_jldirc)
    {
      gROOT->LoadMacro("G4_JLDIRC.C");
      JLDIRCInit();
    }

  if (do_barrel_hcal)
    {
      gROOT->LoadMacro("G4_Barrel_Hcal.C");
      Barrel_HcalInit();
    }

  if (do_drich)
    {
      gROOT->LoadMacro("G4_DRich.C");
      DRichInit();
    }

  if (do_endcap_electron)
    {
      gROOT->LoadMacro("G4_EndCap_Electron.C");
      EndCap_ElectronInit();
    }

  if (do_endcap_hadron)
    {
      gROOT->LoadMacro("G4_EndCap_Hadron.C");
      EndCap_HadronInit();
    }
  if (do_beamline)
  {
      gROOT->LoadMacro("G4_BeamLine.C");
      BeamLineInit();
  }

}


int G4Setup(const int absorberactive = 0,
            const bool do_gem = true,
            const bool do_jldirc = true,
            const bool do_barrel_hcal = true,
            const bool do_drich = true,
            const bool do_endcap_electron = true,
            const bool do_endcap_hadron = true,
            const bool do_beamline = true)
 {
  
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();


  PHG4Reco* g4Reco = new PHG4Reco();

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
  if (stringline.fail()) { // conversion to double fails -> we have a string

      g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::kFieldCleo);
  } else {
    g4Reco->set_field(fieldstrength); // use const soleniodal field
  }
  g4Reco->set_field_rescale(G4MAGNET::magfield_rescale);
  
  double radius = 0.;

  //----------------------------------------
  // PIPE
  if (Enable::PIPE)
  {
    radius = Pipe(g4Reco, radius);
  }

  //----------------------------------------
  // VTX
  if (Enable::VTX) radius = VTX(g4Reco, radius);

  //----------------------------------------
  // CTD
  if (Enable::CTD) radius = CTD(g4Reco, radius);
  
  //----------------------------------------
  // DIRC
  if (do_jldirc) radius = JLDIRC(g4Reco, radius, absorberactive);
  
  //----------------------------------------
  // MAGNET
  
  if (Enable::MAGNET) radius = Magnet(g4Reco, radius);

  //----------------------------------------
  // Barrel Hcal
 
  if (do_barrel_hcal) radius = Barrel_Hcal(g4Reco, radius, 0, absorberactive);


  //----------------------------------------
  // Gem (hadron and electron going side
  
  if (do_gem) double tmp = Gem(g4Reco, radius, 0, absorberactive);

  if (do_drich) double tmp =  DRich(g4Reco, radius, 0, absorberactive);

  if (do_endcap_electron) double tmp =  EndCap_Electron(g4Reco, radius, 0, absorberactive);

  if (do_endcap_hadron) double tmp =  EndCap_Hadron(g4Reco, radius, 0, absorberactive);

  if (do_beamline) double tmp = BeamLine(g4Reco, radius, 0, absorberactive);

  if (Enable::USER)
  {
    UserDetector(g4Reco);
  }

  //----------------------------------------
  // BLACKHOLE
  
  if (Enable::BLACKHOLE)
  {
    BlackHole(g4Reco, radius);
  }

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);
  // finally adjust the world size in case the default is too small
  WorldSize(g4Reco, radius);

  se->registerSubsystem( g4Reco );
  return 0;
}
void DstCompress(Fun4AllDstOutputManager *out)
{
  if (out)
  {
    out->StripNode("G4HIT_PIPE");
  }
}

#endif // MACRO_G4SETUPJLEIC_C
