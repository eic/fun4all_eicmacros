#ifndef MACRO_G4EEMC_C
#define MACRO_G4EEMC_C

#include <GlobalVariables.C>

#include <g4calo/RawTowerBuilderByHitIndex.h>
#include <g4calo/RawTowerDigitizer.h>
#include <g4eiccalos/PHG4CrystalCalorimeterSubsystem.h>
#include <g4eiccalos/PHG4ForwardCalCellReco.h>

#include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>

#include <caloreco/RawClusterBuilderFwd.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawTowerCalibration.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4eiccalos.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)

namespace Enable
{
  bool EEMC = false;
  bool EEMC_ABSORBER = false;
  bool EEMC_CELL = false;
  bool EEMC_TOWER = false;
  bool EEMC_CLUSTER = false;
  bool EEMC_EVAL = false;
  bool EEMC_OVERLAPCHECK = false;
  int EEMC_VERBOSITY = 0;
}  // namespace Enable


namespace G4EEMC
{
  int use_projective_geometry = 0;

  //  double Gdz = 18. + 0.0001;  // These 2 paras are only served as the dimension of the black hole
  //  double Gz0 = -170.;
  double Gdz = 20. + 0.1;
  double Gz0 = -211.;
  
  // Digitization (default photon digi):
  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kSimple_photon_digitization;
  // directly pass the energy of sim tower to digitized tower
  // kNo_digitization
  // simple digitization with photon statistics, single amplitude ADC conversion and pedestal
  // kSimple_photon_digitization
  // digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  // kSiPM_photon_digitization

  enum enu_Eemc_clusterizer
  {
    kEemcGraphClusterizer,
    kEemcTemplateClusterizer
  };
  //default template clusterizer, as developed by Sasha Bazilevsky
  enu_Eemc_clusterizer Eemc_clusterizer = kEemcTemplateClusterizer;
  // graph clusterizer
  //enu_Eemc_clusterizer Eemc_clusterizer = kEemcGraphClusterizer;

}  // namespace G4EEMC


void EEMCInit()
{
  if (G4EEMC::use_projective_geometry)
  {
    BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 81.);
  }
  else
  {
    BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 65.6);
  }
  // from towerMap_EEMC_v006.txt
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4EEMC::Gz0 - G4EEMC::Gdz / 2.);
}


void EEMCSetup(PHG4Reco *g4Reco)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::EEMC_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::EEMC_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::EEMC_VERBOSITY);

  /** Use dedicated EEMC module */
  ostringstream mapping_eemc_1, mapping_eemc_2;
    

  
  PHG4CrystalCalorimeterSubsystem *eemc_crystal = new PHG4CrystalCalorimeterSubsystem("EEMC");
  eemc_crystal->SuperDetector("EEMC");
  eemc_crystal->SetActive();
  if (AbsorberActive)
    eemc_crystal->SetAbsorberActive();
  
  if (!G4EEMC::use_projective_geometry)
  {
    mapping_eemc_1 << getenv("CALIBRATIONROOT") << "/CrystalCalorimeter/mapping/crystal_mapping/tower_map_crystal.txt";
    eemc_crystal->set_string_param("mappingtower", mapping_eemc_1.str());
  }
  eemc_crystal->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(eemc_crystal);
  

  
  
  PHG4CrystalCalorimeterSubsystem *eemc_glass = new PHG4CrystalCalorimeterSubsystem("EEMC_glass");
  eemc_glass->SuperDetector("EEMC_glass");
  eemc_glass->SetActive();
  if (AbsorberActive)
    eemc_glass->SetAbsorberActive();

  if (!G4EEMC::use_projective_geometry)
    {
      mapping_eemc_2 << getenv("CALIBRATIONROOT") << "/CrystalCalorimeter/mapping/crystal_mapping/tower_map_glass.txt";
      eemc_glass->set_string_param("mappingtower", mapping_eemc_2.str());
    }
  
  eemc_glass->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(eemc_glass);
  
  
}

void EEMC_Cells()
{}

void EEMC_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  ostringstream mapping_eemc_1, mapping_eemc_2;
  mapping_eemc_1 << getenv("CALIBRATIONROOT") << "/CrystalCalorimeter/mapping/crystal_mapping/tower_map_crystal.txt";
  mapping_eemc_2 << getenv("CALIBRATIONROOT") << "/CrystalCalorimeter/mapping/crystal_mapping/tower_map_glass.txt";

  // CMS lead tungstate barrel ECAL at 18 degree centrigrade: 4.5 photoelectrons per MeV
  // lead tungsten test in Orsay is 15~20 p.e. per MeV, sci-glass is 5 p.e. per MeV
  const double EEMC_photoelectron_per_GeV_crystal = 15000;
  const double EEMC_photoelectron_per_GeV_glass = 5000;

  //the original values are [8, 16], no noise case[0, 0], really high case[80, 160]
  const double crystal_pedestal_ADC = 0, crystal_zero_suppression_ADC = 0;
  const double glass_pedestal_ADC = 0, glass_zero_suppression_ADC = 0;
  


  
  RawTowerBuilderByHitIndex *tower_EEMC_crystal = new RawTowerBuilderByHitIndex("TowerBuilder_EEMC_crystal");
  tower_EEMC_crystal->Detector("EEMC");
  tower_EEMC_crystal->set_sim_tower_node_prefix("SIM");
  tower_EEMC_crystal->GeometryTableFile(mapping_eemc_1.str());
  se->registerSubsystem(tower_EEMC_crystal);
  
  // Calorimeter digitization 
  RawTowerDigitizer *TowerDigitizer_EEMC_crystal = new RawTowerDigitizer("EEMCRawTowerDigitizer_crystal");
  TowerDigitizer_EEMC_crystal->Detector("EEMC");
  TowerDigitizer_EEMC_crystal->Verbosity(verbosity);
  TowerDigitizer_EEMC_crystal->set_raw_tower_node_prefix("RAW");
  TowerDigitizer_EEMC_crystal->set_digi_algorithm(G4EEMC::TowerDigi);
  TowerDigitizer_EEMC_crystal->set_pedstal_central_ADC(0);
  TowerDigitizer_EEMC_crystal->set_pedstal_width_ADC(crystal_pedestal_ADC);  // eRD1 test beam setting
  TowerDigitizer_EEMC_crystal->set_photonelec_ADC(1);     //not simulating ADC discretization error
  TowerDigitizer_EEMC_crystal->set_photonelec_yield_visible_GeV(EEMC_photoelectron_per_GeV_crystal);
  TowerDigitizer_EEMC_crystal->set_zero_suppression_ADC(crystal_zero_suppression_ADC);  // eRD1 test beam setting
  se->registerSubsystem(TowerDigitizer_EEMC_crystal);

  // Calorimeter calibration 
  RawTowerCalibration *TowerCalibration_EEMC_crystal = new RawTowerCalibration("EEMCRawTowerCalibration_crystal");
  TowerCalibration_EEMC_crystal->Detector("EEMC");
  TowerCalibration_EEMC_crystal->Verbosity(verbosity);
  TowerCalibration_EEMC_crystal->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  if (G4EEMC::TowerDigi == RawTowerDigitizer::kNo_digitization)
    TowerCalibration_EEMC_crystal->set_calib_const_GeV_ADC(1.);
  else
    TowerCalibration_EEMC_crystal->set_calib_const_GeV_ADC(1. / EEMC_photoelectron_per_GeV_crystal);
  TowerCalibration_EEMC_crystal->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration_EEMC_crystal);
  

    
  
  RawTowerBuilderByHitIndex *tower_EEMC_glass = new RawTowerBuilderByHitIndex("TowerBuilder_EEMC_glass");
  tower_EEMC_glass->Detector("EEMC_glass");
  tower_EEMC_glass->set_sim_tower_node_prefix("SIM");
  tower_EEMC_glass->GeometryTableFile(mapping_eemc_2.str());
  se->registerSubsystem(tower_EEMC_glass);

  RawTowerDigitizer *TowerDigitizer_EEMC_glass = new RawTowerDigitizer("EEMCRawTowerDigitizer_glass");
  TowerDigitizer_EEMC_glass->Detector("EEMC_glass");
  TowerDigitizer_EEMC_glass->Verbosity(verbosity);
  TowerDigitizer_EEMC_glass->set_raw_tower_node_prefix("RAW");
  TowerDigitizer_EEMC_glass->set_digi_algorithm(G4EEMC::TowerDigi);
  TowerDigitizer_EEMC_glass->set_pedstal_central_ADC(0);
  TowerDigitizer_EEMC_glass->set_pedstal_width_ADC(glass_pedestal_ADC);  // eRD1 test beam setting
  TowerDigitizer_EEMC_glass->set_photonelec_ADC(1);     //not simulating ADC discretization error
  TowerDigitizer_EEMC_glass->set_photonelec_yield_visible_GeV(EEMC_photoelectron_per_GeV_glass);
  TowerDigitizer_EEMC_glass->set_zero_suppression_ADC(glass_zero_suppression_ADC);  // eRD1 test beam setting
  se->registerSubsystem(TowerDigitizer_EEMC_glass);  
  
  RawTowerCalibration *TowerCalibration_EEMC_glass = new RawTowerCalibration("EEMCRawTowerCalibration_glass");
  TowerCalibration_EEMC_glass->Detector("EEMC_glass");
  TowerCalibration_EEMC_glass->Verbosity(verbosity);
  TowerCalibration_EEMC_glass->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
  if (G4EEMC::TowerDigi == RawTowerDigitizer::kNo_digitization)
    TowerCalibration_EEMC_glass->set_calib_const_GeV_ADC(1.);
  else
    TowerCalibration_EEMC_glass->set_calib_const_GeV_ADC(1. / EEMC_photoelectron_per_GeV_glass);
  TowerCalibration_EEMC_glass->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration_EEMC_glass);
  

  
}



void EEMC_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EEMC_VERBOSITY);
  Fun4AllServer *se = Fun4AllServer::instance();

  
  if (G4EEMC::Eemc_clusterizer == G4EEMC::kEemcTemplateClusterizer)
  {
    
    RawClusterBuilderTemplate *ClusterBuilder_crystal = new RawClusterBuilderTemplate("EEMCRawClusterBuilderTemplate_crystal");
    ClusterBuilder_crystal->Detector("EEMC");
    ClusterBuilder_crystal->Verbosity(2);
    se->registerSubsystem(ClusterBuilder_crystal);
    
    
    RawClusterBuilderTemplate *ClusterBuilder_glass = new RawClusterBuilderTemplate("EEMCRawClusterBuilderTemplate_glass");
    ClusterBuilder_glass->Detector("EEMC_glass");
    ClusterBuilder_glass->Verbosity(verbosity);
    se->registerSubsystem(ClusterBuilder_glass);
    
  }
  else if (G4EEMC::Eemc_clusterizer == G4EEMC::kEemcGraphClusterizer)
  {
        
    RawClusterBuilderFwd *ClusterBuilder_crystal = new RawClusterBuilderFwd("EEMCRawClusterBuilderFwd_crystal");
    ClusterBuilder_crystal->Detector("EEMC");
    ClusterBuilder_crystal->Verbosity(verbosity);
    ClusterBuilder_crystal->Verbosity(2);
    se->registerSubsystem(ClusterBuilder_crystal);
    
    
    RawClusterBuilderFwd *ClusterBuilder_glass = new RawClusterBuilderFwd("EEMCRawClusterBuilderFwd_glass");
    ClusterBuilder_glass->Detector("EEMC_glass");
    ClusterBuilder_glass->Verbosity(verbosity);
    se->registerSubsystem(ClusterBuilder_glass);
    
  }
  else
  {
    cout << "EEMC_Clusters - unknown clusterizer setting " << G4EEMC::Eemc_clusterizer << endl;
    gSystem->Exit(1);
  }
  return;
}


void EEMC_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::EEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

     
  CaloEvaluator *eval_crystal = new CaloEvaluator("EEMCEVALUATOR", "EEMC", outputfile.c_str());
  eval_crystal->Verbosity(verbosity);
  eval_crystal->set_do_cluster_eval(true);
  se->registerSubsystem(eval_crystal);
  
  
  CaloEvaluator *eval_glass = new CaloEvaluator("EEMCEVALUATOR", "EEMC_glass", outputfile.c_str());
  eval_glass->Verbosity(verbosity);
  eval_glass->set_do_cluster_eval(true);
  se->registerSubsystem(eval_glass);  
  
  
  /*
  CaloEvaluator *eval_crystal = new CaloEvaluator("EEMCEVALUATOR", "EEMC_crystal", outputfile.c_str());
  eval_crystal->Verbosity(verbosity);
  eval_crystal->set_do_cluster_eval(true);
  se->registerSubsystem(eval_crystal);
  
   
  CaloEvaluator *eval_glass = new CaloEvaluator("EEMCEVALUATOR", "EEMC_glass", outputfile.c_str());
  eval_glass->Verbosity(verbosity);
  eval_glass->set_do_cluster_eval(true);
  se->registerSubsystem(eval_glass);  
  */
  return;
}
#endif
