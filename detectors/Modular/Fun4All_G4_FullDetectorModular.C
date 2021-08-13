#ifndef MACRO_FUN4ALLG4EICDETECTOR_C
#define MACRO_FUN4ALLG4EICDETECTOR_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_ModularDetector.C>
#include <G4_Bbc.C>
#include <G4_CaloTrigger.C>
#include <G4_DSTReader_ModularDetector.C>
#include <G4_FwdJets.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_Production.C>
#include <G4_User.C>
#include <QA.C>

#include <TROOT.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <PHPy6GenTrigger.h>
#include <PHPy6ParticleTrigger.h>
#include <PHPy6JetTrigger.h>
#include <phool/recoConsts.h>

#include <eiceval/EventEvaluatorEIC.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libeiceval.so)

void ParseTString(TString &specialSetting);

int Fun4All_G4_FullDetectorModular(
    const int nEvents                 = 1,
    const double particlemomMin       = -1,
    const double particlemomMax       = -1,
    TString specialSetting            = "ALLSILICON-TTLEM",
    TString generatorSettings         = "e10p250MB",
    const string &inputFile           = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const string &outputFile          = "G4EICDetector.root",
    const string &embed_input_file    = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip                    = 0,
    const string &outdir              = ".")
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  //PHRandomSeed::Verbosity(1);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as initial seed
  // which will produce identical results so you can debug your code
  //rc->set_IntFlag("RANDOMSEED", 12345);
  
  // switching IPs by comment/uncommenting the following lines
  // used for both beamline setting and for the event generator crossing boost
  Enable::IP6 = true;
  // Enable::IP8 = true;
  
  //===============
  // Input options
  //===============

  // pythia6
  // Use Pythia 6
  if(particlemomMin==-1 && particlemomMax==-1){
    Input::PYTHIA6 = true;
  }
  // Simple multi particle generator in eta/phi/pt ranges
  Input::SIMPLE = false;
  if (particlemomMin>-1 && particlemomMax>-1){
    Input::SIMPLE = true;
    Input::SIMPLE_VERBOSITY = 0;
 }
  // Input::SIMPLE_NUMBER = 2; // if you need 2 of them
 
  Input::VERBOSITY = 0;
  INPUTHEPMC::filename = inputFile;


  Enable::QA = false;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  InputInit();
  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::SIMPLE){
    if (generatorSettings.Contains("SimplePion"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 1);
    else if (generatorSettings.Contains("SimpleKaon"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("kaon-", 1);
    else if (generatorSettings.Contains("SimpleProton"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("proton", 1);
    else if (generatorSettings.Contains("SimplePhoton"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("gamma", 1);
    else if (generatorSettings.Contains("SimpleNeutron"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("neutron", 1);
    else if (generatorSettings.Contains("SimpleElectron"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("e-", 1);    
    else if (generatorSettings.Contains("SimplePiZero"))
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi0", 1);
    else {
      std::cout << "You didn't specify which particle you wanted to generate, exiting" << std::endl;
      return 0;
    }
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                              PHG4SimpleEventGenerator::Uniform,
                                                                              PHG4SimpleEventGenerator::Uniform);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 5.);
    if (generatorSettings.Contains("central"))
      INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.8, 1.2);
    else if (generatorSettings.Contains("bck"))
      INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-4, -1.7);
    else 
      INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(1.1, 4.0);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_p_range(particlemomMin, particlemomMax);
  }
  if(particlemomMin>-1 && particlemomMax == -1){
    PHG4ParticleGenerator *gen = new PHG4ParticleGenerator("PGENERATOR");
    gen->set_name("pi-");
    // gen->set_name("pi0");
    gen->set_vtx(0, 0, 0);
    gen->set_eta_range(-3.5, 3.5);            // around midrapidity
    if(particlemomMin > -1)
      gen->set_mom_range(particlemomMin, particlemomMin);                   // fixed 4 GeV/c
    else
      gen->set_mom_range(1, 60);                   // fixed 4 GeV/c
    gen->set_phi_range(0., 2* M_PI);  // 0-90 deg
    // gen->Verbosity(1);  // 0-90 deg
    se->registerSubsystem(gen);
  }
  // pythia6
  if (Input::PYTHIA6){
    if (generatorSettings.Contains("e10p250MB") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep.cfg");
    else if (generatorSettings.Contains("e10p250pTHard5") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep_MinPartonP5GeV.cfg");
    else if (generatorSettings.Contains("e10p250pTQ210"))
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep_QSquare10GeV.cfg");
    else if (generatorSettings.Contains("e10p250pTHard10"))
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep_MinPartonP10GeV.cfg");
    else if (generatorSettings.Contains("e10p250pTHard20"))
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep_MinPartonP20GeV.cfg");
    else if (generatorSettings.Contains("e5p100MB") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e5p100.cfg");
    else if (generatorSettings.Contains("e5p100pTHard5") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e5p100_MinPartonP5GeV.cfg");
    else if (generatorSettings.Contains("e10p275MB") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e10p275.cfg");
    else if (generatorSettings.Contains("e10p275pTHard5") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e10p275_MinPartonP5GeV.cfg");
    else if (generatorSettings.Contains("e10p275pTHard10") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e10p275_MinPartonP10GeV.cfg");
    else if (generatorSettings.Contains("e18p275MB") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e18p275.cfg");
    else if (generatorSettings.Contains("e18p275pTHard5") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e18p275_MinPartonP5GeV.cfg");
    else if (generatorSettings.Contains("e18p275pTHard10") )
      INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_e18p275_MinPartonP10GeV.cfg");
    else 
      INPUTGENERATOR::Pythia6->set_config_file(generatorSettings.Data());
    
    if (generatorSettings.Contains("FPartTrigg")){
      PHPy6ParticleTrigger *ptrig = new PHPy6ParticleTrigger();
      ptrig->SetPtLow(1);
      ptrig->SetEtaHighLow(1,5);
      INPUTGENERATOR::Pythia6->register_trigger(ptrig);
    }
    if (generatorSettings.Contains("FJetTrigg")){
      PHPy6JetTrigger *trig = new PHPy6JetTrigger();
      trig->SetEtaHighLow(1,5);
      if (generatorSettings.Contains("pTHard5"))
        trig->SetMinJetPt(5);
      else if (generatorSettings.Contains("pTHard10"))
        trig->SetMinJetPt(10);
      else if (generatorSettings.Contains("pTHard20"))
        trig->SetMinJetPt(20);
      else 
        trig->SetMinJetPt(1);
      trig->SetJetR(0.7);
      INPUTGENERATOR::Pythia6->register_trigger(trig);
    }
    //! apply EIC beam parameter following EIC CDR
    Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia6);

    
  }

  // register all input generators with Fun4All
  InputRegister();


  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;
  Enable::DSTOUT_COMPRESS = false;  // Compress DST files

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  Enable::DSTREADER = false;

  // turn the display on (default off)
  if(specialSetting.Contains("display"))
    Enable::DISPLAY = true;

  bool enableCentral = true;
  bool enableForwardTracking = true;
  bool enableElectronSide = true;

  //======================
  // What to run
  //======================
  // Global options (enabled for all subsystems - if implemented)
  //  Enable::ABSORBER = true;
//    Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  //  Enable::BBC = true;
  Enable::BBCFAKE = true; // Smeared vtx and t0, use if you don't want real BBC in simulation

  // whether to simulate the Be section of the beam pipe
  Enable::PIPE = true;
  // EIC beam pipe extension beyond the Be-section:
  G4PIPE::use_forward_pipes = true;
  //EIC hadron far forward magnets and detectors. IP6 and IP8 are incompatible (pick either or);
  Enable::HFARFWD_MAGNETS = false;
  Enable::HFARFWD_VIRTUAL_DETECTORS = false;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // geometry - tracking
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // backward GEM
  if (specialSetting.Contains("EGEM")){
    Enable::EGEM = true;
    if (specialSetting.Contains("EGEMOO")) // only last 2 EGEM layers
      Enable::EGEM_FULL = false;
    else 
      Enable::EGEM_FULL = true;
  }
  //Forward GEM
  if (specialSetting.Contains("FGEM")){
    Enable::FGEM = true;
     // FGEM settings
    if (specialSetting.Contains("FGEMOrg")){
      Enable::FGEM_ORIG = true;
    } else {
      Enable::FGEM_ORIG = false;
    }
  }
  
  // barrel tracker (LANL)
  if (specialSetting.Contains("BARREL"))
    Enable::BARREL = true;
  if(specialSetting.Contains("FST"))
    Enable::FST = true;
  
  // all silicon tracker version (LBL)
  if(specialSetting.Contains("ALLSILICON")){
    Enable::ALLSILICON = true;
    Enable::ALLSILICON_ABSORBER = true;
  }
  
  // LGAD layers
  if(specialSetting.Contains("TTL")){
    Enable::FTTL = true;
    Enable::ETTL = true;
    Enable::CTTL = true;
    G4TTL::SETTING::optionCEMC    = true;
  }
  // mvtx/tpc tracker
  if(specialSetting.Contains("MVTX")){
    Enable::MVTX = true;
    Enable::TPC = true;
  }
  
  if(specialSetting.Contains("ENDCAPTPC"))
   Enable::TPC_ENDCAP = true;

  Enable::TRACKING = true;
  Enable::TRACKING_EVAL = Enable::TRACKING && true;
  if (specialSetting.Contains("INNERTRACKING")) {
    Enable::TRACKING_INNER = true;
  }
  if (specialSetting.Contains("TREXTOUT"))
    Enable::TRACKING_EVAL_DETAILED = Enable::TRACKING_EVAL && true;

  G4TRACKING::DISPLACED_VERTEX = true;  // this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2ndary vertexes
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // geometry - barrel
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // PID detectors
  Enable::DIRC = true;
  G4DIRC::SETTING::USECEMCGeo   = true;
  
  // sPHENIX SPACAL reuse
  Enable::CEMC = true;
  //  Enable::CEMC_ABSORBER = true;
  
  // sPHENIX HCal inner reuse
  Enable::HCALIN = true;
  G4HCALIN::SETTING::USECEMCGeo = true;
  //  Enable::HCALIN_ABSORBER = true;
  
  if (specialSetting.Contains("BECAL") ){
    Enable::BECAL = true;
    // need to switch of CEMC & HCALin
    Enable::CEMC    = false;
    Enable::HCALIN  = true; // for now deactivated due to crash
  }
  Enable::MAGNET = true;

  // sPHENIX HCal outer reuse
  Enable::HCALOUT = true;
  //  Enable::HCALOUT_ABSORBER = true;


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // geometry - 'hadron' direction
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // PID detectors - RICH's
  Enable::RICH = true;
//   Enable::AEROGEL = true;

  // PHENIX EMCal shashlik reuse
  Enable::FEMC = true;
  //  Enable::FEMC_ABSORBER = true;

  // STAR forward HCal
  Enable::FHCAL = true;
  if(specialSetting.Contains("FEMCSTANDALONE") || specialSetting.Contains("LFHCAL"))
    Enable::FHCAL = false;
  Enable::FHCAL_VERBOSITY = 0;
  //  Enable::FHCAL_ABSORBER = true;
  //  Enable::FHCAL_SUPPORT = true; // make support active volume


  Enable::DRCALO = false;
  if(specialSetting.Contains("DRCALO")){
    Enable::DRCALO = true;
    G4TTL::SETTING::optionDR = 1;
    if(!specialSetting.Contains("FwdConfig") && !specialSetting.Contains("FwdSquare")){
      Enable::FEMC = false;
      Enable::FHCAL = false;
    }
  }
  Enable::DRCALO_VERBOSITY = 0;
  //  Enable::DRCALO_ABSORBER = true;

  // PSD like HCal
  if ( specialSetting.Contains("LFHCAL")){
    Enable::LFHCAL = true;
    Enable::LFHCAL_ABSORBER = false;
  }

  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // EICDetector geometry - 'electron' direction
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // PID detectors - RICH's
  Enable::mRICH = true;

  Enable::EEMC  = true;
  Enable::EEMCH = false;
  if (specialSetting.Contains("EEMCH")){
    Enable::EEMCH = true;
    Enable::EEMC  = false;
    G4EEMCH::SETTING::USECEMCGeo  = true;
  }
  Enable::EHCAL = true;
  if(specialSetting.Contains("noEHCAL"))
    Enable::EHCAL = false;
  Enable::EHCAL_VERBOSITY = 0;
  //  Enable::EHCAL_ABSORBER = true;

  Enable::PLUGDOOR = false;

  
  // projections to calorimeters
  if (specialSetting.Contains("TRACKEVALHITS")){
    G4TRACKING::PROJECTION_CEMC   = false;
    G4TRACKING::PROJECTION_FEMC   = false;
    G4TRACKING::PROJECTION_FHCAL  = false;
    G4TRACKING::PROJECTION_EEMC   = false;
    G4TRACKING::PROJECTION_EHCAL  = false;
    G4TRACKING::PROJECTION_DRCALO = false;
  } else {
    G4TRACKING::PROJECTION_CEMC   = Enable::CEMC && true;
    G4TRACKING::PROJECTION_FEMC   = Enable::FEMC && true;
    G4TRACKING::PROJECTION_FHCAL  = Enable::FHCAL && true;
    G4TRACKING::PROJECTION_EEMC   = ( Enable::EEMC || Enable::EEMCH ) && true;  
    G4TRACKING::PROJECTION_EHCAL  = Enable::EHCAL && true;
    G4TRACKING::PROJECTION_DRCALO = Enable::DRCALO && true;
  }


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // special settings for Calo standalone studies
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // deactivate all respective detector systems for standalone studies
  if(specialSetting.Contains("STANDALONE") ){
    Enable::PIPE = false;
    G4PIPE::use_forward_pipes = false;
    Enable::HFARFWD_MAGNETS = false;
    Enable::HFARFWD_VIRTUAL_DETECTORS = false;
    Enable::TPC_ENDCAP = false;
    G4TRACKING::PROJECTION_CEMC   = false;
    G4TRACKING::PROJECTION_FEMC   = false;
    G4TRACKING::PROJECTION_FHCAL  = false;
    G4TRACKING::PROJECTION_EHCAL  = false;
    G4TRACKING::PROJECTION_DRCALO = false;
    G4TRACKING::PROJECTION_EEMC   = false;
    Enable::MAGNET = false;
    Enable::DIRC = false;
    Enable::RICH = false;
    Enable::mRICH = false;
    Enable::AEROGEL = false;
    Enable::CEMC = false;
    Enable::HCALOUT = false;
    Enable::HCALIN = false;
    Enable::EHCAL = false;
    Enable::EEMC = false;
    Enable::EEMCH = false;
    Enable::FEMC = false;
    Enable::FHCAL = false;
    Enable::LFHCAL = false;
    Enable::BECAL = false;
    Enable::FTTL = false;
    Enable::CTTL = false;
    Enable::ETTL = false;
    if(specialSetting.Contains("PIPE")){
      Enable::PIPE = true;
      G4PIPE::use_forward_pipes = true;    
    }
    if(specialSetting.Contains("Magnet"))
      Enable::MAGNET = true;
    if(specialSetting.Contains("ALLSILICON"))
      Enable::ALLSILICON = true;
    if(specialSetting.Contains("CEMC"))
      Enable::CEMC = true;
    if(specialSetting.Contains("DR"))
      Enable::DRCALO = true;
    if(specialSetting.Contains("FEMC"))
      Enable::FEMC = true;
    if(specialSetting.Contains("FHCAL") && !specialSetting.Contains("LFHCAL"))
      Enable::FHCAL = true;
    if(specialSetting.Contains("LFHCAL"))
      Enable::LFHCAL = true;
    if(specialSetting.Contains("BECAL"))
      Enable::BECAL = true;
    if(specialSetting.Contains("EHCAL"))
      Enable::EHCAL = true;
    if(specialSetting.Contains("EEMCH"))
      Enable::EEMCH = true;
    if(specialSetting.Contains("CHCAL")){
      Enable::HCALIN   = true;
      Enable::HCALOUT  = true;
    }
    if(specialSetting.Contains("DIRC"))
      Enable::DIRC = true;
    
    if(specialSetting.Contains("FWDCALO")){
      Enable::FEMC    = true;
      Enable::FHCAL   = true;
    }
    if(specialSetting.Contains("FWDLCALO")){
      Enable::FEMC    = true;
      Enable::LFHCAL  = true;
    }
    if(specialSetting.Contains("BARCALO")){
      Enable::BECAL    = true;
      Enable::HCALIN   = true;
      Enable::HCALOUT  = true;
      Enable::MAGNET = true;
    }
    if(specialSetting.Contains("BCKCALO")){
      Enable::EHCAL    = true;
      Enable::EEMCH    = true;
    }
    if(specialSetting.Contains("TTL")){
      // Enable::PIPE = true;
      // G4PIPE::use_forward_pipes = true;
      // LGAD layers
      if(specialSetting.Contains("FTTL"))
        Enable::FTTL = true;
      if(specialSetting.Contains("ETTL"))
        Enable::ETTL = true;
      if(specialSetting.Contains("CTTL")){
        Enable::CTTL = true;
        // Enable::DIRC = true;
        // Enable::CEMC = true;
        // Enable::BECAL = true;
        // Enable::ALLSILICON = true;
        G4DIRC::SETTING::USECEMCGeo   = false;
        G4TTL::SETTING::optionCEMC    = false;
      }
    }
  }

  // Automatic settings based on previous selections:
  Enable::CEMC_CELL       = Enable::CEMC && true;
  Enable::CEMC_TOWER      = Enable::CEMC_CELL && true;
  Enable::CEMC_CLUSTER    = Enable::CEMC_TOWER && true;
  Enable::CEMC_EVAL       = Enable::CEMC_CLUSTER && false;

  Enable::HCALIN_CELL     = Enable::HCALIN && true;
  Enable::HCALIN_TOWER    = Enable::HCALIN_CELL && true;
  Enable::HCALIN_CLUSTER  = Enable::HCALIN_TOWER && true;
  Enable::HCALIN_EVAL     = Enable::HCALIN_CLUSTER && false;

  Enable::HCALOUT_CELL    = Enable::HCALOUT && true;
  Enable::HCALOUT_TOWER   = Enable::HCALOUT_CELL && true;
  Enable::HCALOUT_CLUSTER = Enable::HCALOUT_TOWER && true;
  Enable::HCALOUT_EVAL    = Enable::HCALOUT_CLUSTER && false;

  Enable::DRCALO_CELL     = Enable::DRCALO && true;
  Enable::DRCALO_TOWER    = Enable::DRCALO_CELL && true;
  Enable::DRCALO_CLUSTER  = Enable::DRCALO_TOWER && true;
  Enable::DRCALO_EVAL     = Enable::DRCALO_CLUSTER && true;

  Enable::FEMC_CELL     = Enable::FEMC && true;
  Enable::FEMC_TOWER    = Enable::FEMC_CELL && true;
  Enable::FEMC_CLUSTER  = Enable::FEMC_TOWER && true;
  Enable::FEMC_EVAL     = Enable::FEMC_CLUSTER && false;

  Enable::FHCAL_CELL    = Enable::FHCAL && true;
  Enable::FHCAL_TOWER   = Enable::FHCAL_CELL && true;
  Enable::FHCAL_CLUSTER = Enable::FHCAL_TOWER && true;
  Enable::FHCAL_EVAL    = Enable::FHCAL_CLUSTER && false;

  Enable::EEMC_CELL     = Enable::EEMC && true;
  Enable::EEMC_TOWER    = Enable::EEMC_CELL && true;
  Enable::EEMC_CLUSTER  = Enable::EEMC_TOWER && true;
  Enable::EEMC_EVAL     = Enable::EEMC_CLUSTER && false;

  Enable::EEMCH_CELL     = Enable::EEMCH && true;
  Enable::EEMCH_TOWER    = Enable::EEMCH_CELL && true;
  Enable::EEMCH_CLUSTER  = Enable::EEMCH_TOWER && true;
  Enable::EEMCH_EVAL     = Enable::EEMCH_CLUSTER && false;

  Enable::EHCAL_CELL    = Enable::EHCAL && true;
  Enable::EHCAL_TOWER   = Enable::EHCAL_CELL && true;
  Enable::EHCAL_CLUSTER = Enable::EHCAL_TOWER && true;
  Enable::EHCAL_EVAL    = Enable::EHCAL_CLUSTER && false;

  Enable::LFHCAL_CELL    = Enable::LFHCAL && true;
  Enable::LFHCAL_TOWER   = Enable::LFHCAL_CELL && true;
  Enable::LFHCAL_CLUSTER = Enable::LFHCAL_TOWER && true;
  Enable::LFHCAL_EVAL    = Enable::LFHCAL_CLUSTER && false;

  Enable::BECAL_CELL    = Enable::BECAL && true;
  Enable::BECAL_TOWER   = Enable::BECAL_CELL && true;
  Enable::BECAL_CLUSTER = Enable::BECAL_TOWER && true;
  Enable::BECAL_EVAL    = Enable::BECAL_CLUSTER && false;
  
  // Other options
  Enable::GLOBAL_RECO = true;
  Enable::GLOBAL_FASTSIM = true;

  Enable::CALOTRIGGER = true && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;

  // Select only one jet reconstruction- they currently use the same
  // output collections on the node tree!
  Enable::JETS = false;
  Enable::JETS_EVAL = Enable::JETS && true;

  Enable::FWDJETS = false;
  Enable::FWDJETS_EVAL = Enable::FWDJETS && true;

  // HI Jet Reco for jet simulations in Au+Au (default is false for
  // single particle / p+p simulations, or for Au+Au simulations which
  // don't care about jets)
  Enable::HIJETS = false && Enable::JETS && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;
  //Enable::BLACKHOLE_SAVEHITS = false; // turn off saving of bh hits
  //BlackHoleGeometry::visible = true;

  //Enable::USER = true;

  //---------------
  // World Settings
  //---------------
  //  G4WORLD::PhysicsList = "QGSP_BERT"; //FTFP_BERT_HP best for calo
  //  G4WORLD::WorldMaterial = "G4_AIR"; // set to G4_GALACTIC for material scans

  //---------------
  // Magnet Settings
  //---------------
  if(specialSetting.Contains("3T")){
    const string magfield = "3.0"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
    G4MAGNET::magfield = magfield;
    //  G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
  } else {
    G4MAGNET::magfield_rescale = -1.4 / 1.5;  // make consistent with expected Babar field strength of 1.4T
  }
  //---------------
  // Pythia Decayer
  //---------------
  // list of decay types in
  // $OFFLINE_MAIN/include/g4decayer/EDecayType.hh
  // default is All:
  // G4P6DECAYER::decayType = EDecayType::kAll;

  // translate the option TString into subsystem namespace options
  ParseTString(specialSetting); 

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------

  // If "readhepMC" is also set, the Upsilons will be embedded in Hijing events, if 'particles" is set, the Upsilons will be embedded in whatever particles are thrown
  if (!Input::READHITS) G4Setup(specialSetting);

  //------------------
  // Detector Division
  //------------------
  if (Enable::BBC || Enable::BBCFAKE) Bbc_Reco();
  if (Enable::CEMC_CELL) CEMC_Cells();
  if (Enable::HCALIN_CELL) HCALInner_Cells();
  if (Enable::HCALOUT_CELL) HCALOuter_Cells();
  if (Enable::EEMC_CELL) EEMC_Cells();
  if (Enable::EEMCH_CELL) EEMCH_Cells();

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------
  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------
  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();

  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();

  //-----------------------------
  // e, h direction Calorimeter  towering and clustering
  //-----------------------------
  if (Enable::FEMC_TOWER) FEMC_Towers();
  if (Enable::FEMC_CLUSTER) FEMC_Clusters();

  if (Enable::FHCAL_TOWER) FHCAL_Towers();
  if (Enable::FHCAL_CLUSTER) FHCAL_Clusters();

  if (Enable::DRCALO_TOWER) DRCALO_Towers();
  if (Enable::DRCALO_CLUSTER) DRCALO_Clusters();

  if (Enable::LFHCAL_TOWER) LFHCAL_Towers();
  if (Enable::LFHCAL_CLUSTER) LFHCAL_Clusters();

  if (Enable::EEMC_TOWER) EEMC_Towers();
  if (Enable::EEMC_CLUSTER) EEMC_Clusters();

  if (Enable::EEMCH_TOWER) EEMCH_Towers();
  if (Enable::EEMCH_CLUSTER) EEMCH_Clusters();

  if (Enable::EHCAL_TOWER) EHCAL_Towers();
  if (Enable::EHCAL_CLUSTER) EHCAL_Clusters();

  if (Enable::BECAL_TOWER) BECAL_Towers();
  if (Enable::BECAL_CLUSTER) BECAL_Clusters();
  
  if (Enable::DSTOUT_COMPRESS) ShowerCompress();

  //--------------
  // SVTX tracking
  //--------------
  if (Enable::TRACKING) Tracking_Reco(specialSetting);

  //-----------------
  // Global Vertexing
  //-----------------
  if (Enable::GLOBAL_RECO) Global_Reco();
  else if (Enable::GLOBAL_FASTSIM) Global_FastSim();

  //-----------------
  // Calo Trigger Simulation
  //-----------------
  if (Enable::CALOTRIGGER) CaloTrigger_Sim();

  //---------
  // Jet reco
  //---------
  if (Enable::JETS) Jet_Reco();
  if (Enable::HIJETS) HIJetReco();
  if (Enable::FWDJETS) Jet_FwdReco();

  string outputroot = outdir + "/" + outputFile;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos){
    outputroot.erase(pos, remove_this.length());
  }

  if (Enable::DSTREADER) G4DSTreader_EICDetector(outputroot + "_DSTReader.root");

  //----------------------
  // Simulation evaluation
  //----------------------
  bool doFullEventTree = true;
  if(doFullEventTree){
    EventEvaluatorEIC *eval = new EventEvaluatorEIC("EVENTEVALUATOR",  outputroot + "_eventtree.root");
    eval->Verbosity(0);
    if(specialSetting.Contains("GEOMETRYTREE"))
      eval->set_do_GEOMETRY(true);
    if (Enable::FHCAL)
      eval->set_do_FHCAL(true);
    if (Enable::FEMC)
      eval->set_do_FEMC(true);
    if (Enable::DRCALO)
      eval->set_do_DRCALO(true);
    if (Enable::LFHCAL)
      eval->set_do_LFHCAL(true);
    if (Enable::EHCAL)
      eval->set_do_EHCAL(true);
    if (Enable::EEMC || Enable::EEMCH)
      eval->set_do_EEMC(true);
    if (Enable::EEMCH && G4EEMCH::SETTING::USEHYBRID)
      eval->set_do_EEMCG(true);
    if (Enable::CEMC)
      eval->set_do_CEMC(true);
    if (Enable::HCALIN)
      eval->set_do_HCALIN(true);
    if (Enable::BECAL)
      eval->set_do_BECAL(true);
    if (Enable::HCALOUT)
      eval->set_do_HCALOUT(true);
    if (Enable::FHCAL || Enable::FEMC || Enable::EHCAL || Enable::EEMC ||  Enable::EEMCH || Enable::CEMC || Enable::HCALIN || Enable::HCALOUT )
      eval->set_do_CLUSTERS(true);
    if (Enable::TRACKING){
      eval->set_do_TRACKS(true);
      eval->set_do_HITS(true);
      eval->set_do_PROJECTIONS(true);
      if (G4TRACKING::DISPLACED_VERTEX) eval->set_do_VERTEX(true);
    }
    if (Input::PYTHIA6){
      eval->set_do_HEPMC(true);
      eval->set_do_store_event_level_info(true);
    }
    eval->set_do_MCPARTICLES(true);
    se->registerSubsystem(eval);
  }

  if (specialSetting.Contains("TRACKEVALHITS")) Tracking_Eval(outputroot + "_g4tracking_eval.root", specialSetting);
  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  //--------------
  // Set up Output Manager
  //--------------
  if (Enable::PRODUCTION){
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT){
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }

  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY){
    DisplayOn();
    // gROOT->ProcessLine("PHG4Reco *g4 = QTGui();"); // alternative to DisplayOn
    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }
  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0){
    return 0;
  }
  // if we run any of the particle generators and use 0 it'll run forever
  if (nEvents == 0 && !Input::READHITS && !Input::HEPMC && !Input::READEIC){
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->skip(skip);
  se->run(nEvents);

  if (Enable::QA) QA_Output(outputroot + "_qa.root");

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION){
    Production_MoveOutput();
  }
  gSystem->Exit(0);
  return 0;
}

void ParseTString(TString &specialSetting)
{
// Barrel settings
  if (specialSetting.Contains("BARRELV1"))
  {
    G4BARREL::SETTING::BARRELV1 = true;
  }
  else if (specialSetting.Contains("BARRELV2"))
  {
    G4BARREL::SETTING::BARRELV2 = true;
  }
  else if (specialSetting.Contains("BARRELV3"))
  {
    G4BARREL::SETTING::BARRELV3 = true;
  }
  else if (specialSetting.Contains("BARRELV4"))
  {
    G4BARREL::SETTING::BARRELV4 = true;
  }
  else if (specialSetting.Contains("BARREL"))
  {
    G4BARREL::SETTING::BARRELV0 = true;
  }

// FST settings
  if (specialSetting.Contains("FSTV1"))
  {
    G4FST::SETTING::FSTV1 = true;
  }
  else if (specialSetting.Contains("FSTV2"))
  {
    G4FST::SETTING::FSTV2 = true;
  }
  else if (specialSetting.Contains("FSTV3"))
  {
    G4FST::SETTING::FSTV3 = true;
  }
  else if (specialSetting.Contains("FSTV4"))
  {
    G4FST::SETTING::FSTV4 = true;
  }
  else if (specialSetting.Contains("FSTV41"))
  {
    G4FST::SETTING::FSTV41 = true;
  }
  else if (specialSetting.Contains("FSTV42"))
  {
    G4FST::SETTING::FSTV42 = true;
  }
  else if (specialSetting.Contains("FSTVTPC"))
  {
    G4FST::SETTING::FST_TPC = true;
  }
  else if (specialSetting.Contains("FST"))
  {
    G4FST::SETTING::FSTV0 = true;
  }

  // FHCAL/FEMC settings
  if (specialSetting.Contains("fsPHENIX"))
  {
    G4FEMC::SETTING::fsPHENIX = Enable::FEMC && true;
  }
  else if (specialSetting.Contains("EC2x"))
  {
    G4FEMC::SETTING::EC2x = Enable::FEMC && true;
  }
  else if (specialSetting.Contains("ROS"))
  {
    G4FEMC::SETTING::readoutsplit = Enable::FEMC && true;
  }
  
  if (specialSetting.Contains("FullEtaAcc")) // common for FHCAL and FEMC
  {
    G4FHCAL::SETTING::FullEtaAcc = Enable::FHCAL && true;
    G4FEMC::SETTING::FullEtaAcc = true;
  }
  if (specialSetting.Contains("ASYM")) // common for FHCAL and FEMC
  {
    G4FHCAL::SETTING::asymmetric = Enable::FHCAL && true;
    G4LFHCAL::SETTING::asymmetric = Enable::LFHCAL &&true;
    G4FEMC::SETTING::asymmetric = Enable::FEMC && true;
  }
  if (specialSetting.Contains("XDEPTH")) // common for FHCAL and FEMC
  {
    G4FHCAL::SETTING::extradepth  = Enable::FHCAL && true;
    G4LFHCAL::SETTING::longer     = Enable::LFHCAL && true;
  }
  if (specialSetting.Contains("wDR")) // common for FHCAL and FEMC
  {
    G4FHCAL::SETTING::wDR = Enable::FHCAL && true;
    G4LFHCAL::SETTING::wDR = Enable::LFHCAL && true;
    G4FEMC::SETTING::wDR = Enable::FEMC && true;
  }

  if (specialSetting.Contains("FwdConfig")) // common for FHCAL and FEMC
  {
    G4DRCALO::SETTING::FwdConfig = Enable::DRCALO && true;
    G4FHCAL::SETTING::wDR = Enable::FHCAL && true;
    G4LFHCAL::SETTING::wDR = Enable::LFHCAL && true;
    G4FEMC::SETTING::wDR = Enable::FEMC && true;
    G4TTL::SETTING::optionDR = 1;
  }
  if (specialSetting.Contains("FwdSquare")) // common for FHCAL and FEMC
  {
    G4DRCALO::SETTING::FwdSquare = Enable::DRCALO &&true;
    G4FHCAL::SETTING::FwdSquare = Enable::DRCALO && true;
    G4LFHCAL::SETTING::FwdSquare = Enable::LFHCAL && true;
    G4FEMC::SETTING::FwdSquare = Enable::FEMC &&true;
    G4TTL::SETTING::optionDR = 1;
  }
  
  if (specialSetting.Contains("HC2x"))
  {
    G4FHCAL::SETTING::HC2x = true;
    G4LFHCAL::SETTING::HC2x = true;
  } 
  if (specialSetting.Contains("FHCFeTungsten"))
  {
    G4FHCAL::SETTING::Absorber_FeTungsten = 1;
  } 
  if (specialSetting.Contains("FHCFeTungsten"))
  {
    G4FHCAL::SETTING::Absorber_FeTungsten = 1;
  } 
  else if (specialSetting.Contains("HC4x"))
  {
    G4FHCAL::SETTING::HC4x = Enable::FHCAL && true;
  }

  if (specialSetting.Contains("towercalib1"))
  {
    G4FHCAL::SETTING::towercalib1 = Enable::FHCAL &&true;
  }
  else if (specialSetting.Contains("towercalibSiPM"))
  {
    G4FHCAL::SETTING::towercalibSiPM = Enable::FHCAL &&true;
  }
  else if (specialSetting.Contains("towercalibHCALIN"))
  {
    G4FHCAL::SETTING::towercalibHCALIN = Enable::FHCAL &&true;
  }
  else if (specialSetting.Contains("towercalib3"))
  {
    G4FHCAL::SETTING::towercalib3 = Enable::FHCAL &&true;
  }
  // DRCALO settings
  if (specialSetting.Contains("DRTungsten"))
  {
    G4DRCALO::SETTING::Tungsten = Enable::DRCALO && true;
  }
  if (specialSetting.Contains("DRQuartz"))
  {
    G4DRCALO::SETTING::Quartz = Enable::DRCALO && true;
  }
  if (specialSetting.Contains("DRPMMA"))
  {
    G4DRCALO::SETTING::PMMA = Enable::DRCALO &&true;
  }  
  if (specialSetting.Contains("DRTUBES"))
  {
    G4DRCALO::SETTING::Tubes = Enable::DRCALO &&true;
  }  
  
  // EEMCH setting
  if (specialSetting.Contains("purePbWO4"))
    G4EEMCH::SETTING::USEHYBRID = false;
  else 
    G4EEMCH::SETTING::USEHYBRID = true;
  
  if (specialSetting.Contains("BECAL")){
    G4EEMCH::SETTING::USECEMCGeo  = false;
    G4DIRC::SETTING::USECEMCGeo   = false;
    G4TTL::SETTING::optionCEMC    = false;
    G4HCALIN::SETTING::USECEMCGeo = false;
  }
  if (specialSetting.Contains("EEMCH"))
    G4TTL::SETTING::optionEEMCH   = true;
  else 
    G4TTL::SETTING::optionEEMCH   = false;
  
  if (specialSetting.Contains("TTLEMd"))
    G4TTL::SETTING::optionGeo    = 1;
  else if (specialSetting.Contains("TTLEMl"))
    G4TTL::SETTING::optionGeo    = 2;
  else if (specialSetting.Contains("TTLEMs"))
    G4TTL::SETTING::optionGeo    = 3;
  else if (specialSetting.Contains("TTLF"))
    G4TTL::SETTING::optionGeo    = 4;

  if (specialSetting.Contains("TTLBasicGeo")){
    G4TTL::SETTING::optionBasicGeo    = true;
  } else {
    // deactivate DIRC basic supports in case the updated TTL is used -> already contains supports
    G4DIRC::SETTING::USEskinSupports = false;
  }

  if (specialSetting.Contains("ACLGAD"))
    G4TTL::SETTING::optionGran    = 2;
  else if (specialSetting.Contains("LGLGAD"))
    G4TTL::SETTING::optionGran    = 3;
  
  if (specialSetting.Contains("ALLSILICONV3"))
    G4ALLSILICON::SETTING::geomVersion = 3;
  
}

#endif
