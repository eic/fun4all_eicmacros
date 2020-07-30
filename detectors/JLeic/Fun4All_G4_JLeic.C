#ifndef MACRO_FUN4ALLG4JLEIC_C
#define MACRO_FUN4ALLG4JLEIC_C

#include "GlobalVariables.C"
#include "GlobalVariables_JLEIC.C"

#include "G4Setup_JLeic.C"
#include "DisplayOn.C"
#include "G4_Tracking_EIC.C"
#include "G4_Input.C"

#include <phool/PHRandomSeed.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <g4histos/G4HitNtuple.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4ParticleGeneratorVectorMeson.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/HepMCNodeReader.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <phool/recoConsts.h>
#include <phpythia6/PHPythia6.h>
#include <phpythia8/PHPythia8.h>
#include <phhepmc/Fun4AllHepMCPileupInputManager.h>
#include <phhepmc/Fun4AllHepMCInputManager.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libphhepmc.so)
R__LOAD_LIBRARY(libPHPythia6.so)
R__LOAD_LIBRARY(libPHPythia8.so)
R__LOAD_LIBRARY(libg4histos.so)

using namespace std;


int Fun4All_G4_JLeic(
    const int nEvents = 1,
    const string &inputFile = "/sphenix/data/data02/review_2017-08-02/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_e-_eta0_8GeV-0002.root",
    const string &outputFile = "G4JLeic.root",
    const string &embed_input_file = "/sphenix/data/data02/review_2017-08-02/sHijing/fm_0-4.list",
    const int skip = 0,
    const string &outdir = ".")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You can either set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);

  //===============
  // Input options
  //===============
  // First enable the input generators
  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  //  Input::READHITS = true;
  INPUTREADHITS::filename = inputFile;

  // Or:
  // Use particle generator
  // And
  // Further choose to embed newly simulated events to a previous simulation. Not compatible with `readhits = true`
  // In case embedding into a production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
  // E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
  //  Input::EMBED = true;
  INPUTEMBED::filename = embed_input_file;

  Input::SIMPLE = true;
  //Input::SIMPLE_VERBOSITY = 1;

  //  Input::PYTHIA6 = true;

  // Input::PYTHIA8 = true;

  //  Input::GUN = true;
  //Input::GUN_VERBOSITY = 1;

  // Upsilon generator
  //  Input::UPSILON = true;
  Input::UPSILON_VERBOSITY = 0;

//  Input::HEPMC = true;
  Input::VERBOSITY = 0;
  INPUTHEPMC::filename = inputFile;

  // Event pile up simulation with collision rate in Hz MB collisions.
  //Input::PILEUPRATE = 100e3;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  // This creates the input generator(s)
  InputInit();

  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  if (Input::SIMPLE)
  {
    INPUTGENERATOR::SimpleEventGenerator->add_particles("pi-", 5);
    if (Input::HEPMC || Input::EMBED)
    {
      INPUTGENERATOR::SimpleEventGenerator->set_reuse_existing_vertex(true);
      INPUTGENERATOR::SimpleEventGenerator->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
    }
    else
    {
      INPUTGENERATOR::SimpleEventGenerator->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                             PHG4SimpleEventGenerator::Uniform,
                                                                             PHG4SimpleEventGenerator::Uniform);
      INPUTGENERATOR::SimpleEventGenerator->set_vertex_distribution_mean(0., 0., 0.);
      INPUTGENERATOR::SimpleEventGenerator->set_vertex_distribution_width(0., 0., 5.);
    }
    INPUTGENERATOR::SimpleEventGenerator->set_eta_range(-1, 1);
    INPUTGENERATOR::SimpleEventGenerator->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator->set_pt_range(0.1, 20.);
  }
  // Upsilons
  if (Input::UPSILON)
  {
    INPUTGENERATOR::VectorMesonGenerator->add_decay_particles("e", 0);
    INPUTGENERATOR::VectorMesonGenerator->set_rapidity_range(-1, 1);
    INPUTGENERATOR::VectorMesonGenerator->set_pt_range(0., 10.);
    // Y species - select only one, last one wins
    INPUTGENERATOR::VectorMesonGenerator->set_upsilon_1s();
  }
  // particle gun
  if (Input::GUN)
  {
    INPUTGENERATOR::Gun->AddParticle("pi-", 0, 1, 0);
    INPUTGENERATOR::Gun->set_vtx(0, 0, 0);
  }

  //--------------
  // Set Input Manager specific options
  //--------------
  // can only be set after InputInit() is called

  if (Input::HEPMC)
  {
    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_width(100e-4,100e-4,30,0);//optional collision smear in space, time
//    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_mean(0,0,0,0);//optional collision central position shift in space, time
    // //optional choice of vertex distribution function in space, time
    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_function(PHHepMCGenHelper::Gaus,PHHepMCGenHelper::Gaus,PHHepMCGenHelper::Gaus,PHHepMCGenHelper::Gaus);
    //! embedding ID for the event
    //! positive ID is the embedded event of interest, e.g. jetty event from pythia
    //! negative IDs are backgrounds, .e.g out of time pile up collisions
    //! Usually, ID = 0 means the primary Au+Au collision background
    //INPUTMANAGER::HepMCInputManager->set_embedding_id(2);
  }
  // register all input generators with Fun4All
  InputRegister();

// set up production relatedstuff
//   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

//  Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  //  Enable::DSTREADER = true;

  // turn the display on (default off)
  Enable::DISPLAY = false;

  //======================
  // What to run
  //======================

  // whether to simulate the Be section of the beam pipe
  Enable::PIPE = true;
  // EIC beam pipe extension beyond the Be-section:
  //G4PIPE::use_forward_pipes = true;

  bool do_ctd = true;

  bool do_vtx = true;

  bool do_jldirc = true;

  Enable::MAGNET = true;

  bool do_barrel_hcal = true;

  bool do_gem = true;

  bool do_drich = true;

  bool do_endcap_electron = true;

  bool do_endcap_hadron = false;

  bool do_beamline = true;

  bool do_tracking = true;
  //---------------
  // Load libraries
  //---------------

  // establish the geometry and reconstruction setup
  G4Init(do_ctd, do_vtx, do_gem, do_jldirc, do_barrel_hcal, do_drich, do_endcap_electron, do_endcap_hadron, do_beamline);

  int absorberactive = 1;  // set to 1 to make all absorbers active volumes

//  const string magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"); // default map from the calibration database
  const string magfield = "2.0"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
  const float magfield_rescale = 1; // scale map if needed

  if (!Input::READHITS)
  {
    //---------------------
    // Detector description
    //---------------------

    G4Setup(absorberactive, magfield,
            do_ctd, do_vtx, do_gem, do_jldirc, do_barrel_hcal, do_drich, do_endcap_electron, do_endcap_hadron, do_beamline, magfield_rescale);
  }


  if (do_tracking)
  {
    Tracking_Reco();
    Tracking_Eval("trkeval.root");
  }
   G4HitNtuple *g4h = new G4HitNtuple("G4HitNtuple","/phenix/scratch/pinkenbu/g4hitntuple.root");
   g4h->AddNode("PIPE", 0);
   g4h->AddNode("JLVTX",1);
   g4h->AddNode("JLCTD",2);
   g4h->AddNode("JLDIRC",3);
   g4h->AddNode("MAGNET",4);
   g4h->AddNode("BARRELHCAL",5);
   se->registerSubsystem(g4h);

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }

  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  if (Enable::DISPLAY)
    {
      DisplayOn();
      // prevent macro from finishing so can see display
      int i;
      cout << "***** Enter any integer to proceed" << endl;
      cin >> i;
    }

  se->run(nEvents);

  //-----
  // Exit
  //-----


  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}

#endif // MACRO_FUN4ALLG4JLEIC_C
