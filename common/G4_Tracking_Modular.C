#ifndef MACRO_G4TRACKINGMODULAR_C
#define MACRO_G4TRACKINGMODULAR_C

#include <GlobalVariables.C>

#include <G4_CEmc_EIC.C>
#include <G4_FEMC_EIC.C>
#include <G4_FHCAL.C>
#include <G4_GEM_EIC.C>
#include <G4_Mvtx_EIC.C>
#include <G4_TPC_EIC.C>
#include <G4_AllSilicon.C>
#include <G4_TTL_EIC.C>

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <trackreco/PHRaveVertexing.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <fun4all/Fun4AllServer.h>

#include <vector>

R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

namespace Enable
{
  bool TRACKING = false;
  bool TRACKING_EVAL = false;
  bool TRACKING_EVAL_DETAILED = false;
  int TRACKING_VERBOSITY = 0;
}  // namespace Enable

namespace G4TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_CEMC = false;
  bool PROJECTION_FEMC = false;
  bool PROJECTION_FHCAL = false;
}  // namespace G4TRACKING

//-----------------------------------------------------------------------------//
void TrackingInit()
{
  TRACKING::TrackNodeName = "TrackMap";
}
//-----------------------------------------------------------------------------//
void Tracking_Reco(TString specialSetting = "")
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
  //  kalman->Verbosity();
  //  kalman->Smearing(false);
  if (G4TRACKING::DISPLACED_VERTEX){
    // do not use truth vertex in the track fitting,
    // which would lead to worse momentum resolution for prompt tracks
    // but this allows displaced track analysis including DCA and vertex finding
    kalman->set_use_vertex_in_fitting(false);
    kalman->set_vertex_xy_resolution(0);  // do not smear the vertex used in the built-in DCA calculation
    kalman->set_vertex_z_resolution(0);   // do not smear the vertex used in the built-in DCA calculation
    kalman->enable_vertexing(true);       // enable vertex finding and fitting
  } else {
    // constraint to a primary vertex and use it as part of the fitting level arm
    kalman->set_use_vertex_in_fitting(true);
    kalman->set_vertex_xy_resolution(50e-4);
    kalman->set_vertex_z_resolution(50e-4);
  }

  kalman->set_sub_top_node_name("TRACKS");
  kalman->set_trackmap_out_name(TRACKING::TrackNodeName);

  //-------------------------
  // Barrel upgrade (LANL)
  //-------------------------
  if (Enable::BARREL){
    kalman->add_phg4hits("G4HIT_BARREL",              //      const std::string& phg4hitsNames,
                         PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                         5e-4,                        //      const float radres,
                         5e-4,                        //      const float phires,
                         5e-4,                        //      const float lonres,
                         1,                           //      const float eff,
                         0);                          //      const float noise
  }
  
  //-------------------------
  // FST (LANL)
  //-------------------------
  if (Enable::FST) {
    for (int i = 0; i < 5; i++) {
      kalman->add_phg4hits(Form("G4HIT_FST_%d", i),           //      const std::string& phg4hitsNames,
                           PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                           5e-4,                              //      const float radres,
                           5e-4,                              //      const float phires,
                           50e-4 / sqrt(12.),                 //      const float lonres,
                           1,                                 //      const float eff,
                           0);                                //      const float noise
    }
  }

  
  //-------------------------
  // Barrel ALLSILICON version (LBL)
  //-------------------------
  char nodename[100];
  if (Enable::ALLSILICON){
    // CENTRAL BARREL
    for (int i = 10; i < 16; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_CENTRAL_%d", i);
      kalman->add_phg4hits(
          nodename,                    // const std::string& phg4hitsNames
          PHG4TrackFastSim::Cylinder,  // const DETECTOR_TYPE phg4dettype
          999.,                        // radial-resolution [cm] (this number is not used in cylindrical geometry)
          5.8e-4,                      // azimuthal (arc-length) resolution [cm]
          5.8e-4,                      // longitudinal (z) resolution [cm]
          1,                           // efficiency (fraction)
          0                            // hit noise
      );
    }
    // FORWARD DISKS
    for (int i = 20; i < 25; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_FORWARD_%d", i);
      kalman->add_phg4hits(
          nodename,                          // const std::string& phg4hitsNames
          PHG4TrackFastSim::Vertical_Plane,  // const DETECTOR_TYPE phg4dettype
          5.8e-4,                            // radial-resolution [cm]
          5.8e-4,                            // azimuthal (arc-length) resolution [cm]
          999.,                              // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
          1,                                 // efficiency (fraction)
          0                                  // hit noise
      );
    }
    // BACKWARD DISKS
    for (int i = 30; i < 35; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_BACKWARD_%d", i);
      kalman->add_phg4hits(
          nodename,                          // const std::string& phg4hitsNames
          PHG4TrackFastSim::Vertical_Plane,  // const DETECTOR_TYPE phg4dettype
          5.8e-4,                            // radial-resolution [cm]
          5.8e-4,                            // azimuthal (arc-length) resolution [cm]
          999.,                              // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
          1,                                 // efficiency (fraction)
          0                                  // hit noise
      );
    }
  }
  
  //-------------------------
  // Timing tracking Layer (TTL - ORNL/Rice)
  //-------------------------
  // central barrel 
  if (Enable::CTTL){
    Float_t pitch=500e-4;
    for (int i = 0; i < 2; i++){
      kalman->add_phg4hits(Form("G4HIT_CTTL_%d", i),           //      const std::string& phg4hitsNames,
                          PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                          999,                               //      const float radres,
                          pitch/sqrt(12),                    //      const float phires,
                          pitch/sqrt(12),                    //      const float lonres, *ignored in plane detector*
                          0.95,                              //      const float eff,
                          0);                                //      const float noise
    }
    kalman -> add_cylinder_state("CTTL_0", 92);
    kalman -> add_cylinder_state("CTTL_1", 114.7);
  }
  
  // electron going direction
  if (Enable::ETTL){
    Float_t pitch=500e-4;
    for (int i = 0; i < 2; i++){
      kalman->add_phg4hits(Form("G4HIT_ETTL_%d", i),           //      const std::string& phg4hitsNames,
                          PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                          pitch/sqrt(12),                    //      const float radres,
                          pitch/sqrt(12),                    //      const float phires,
                          999.,                              //      const float lonres, *ignored in plane detector*
                          0.95,                              //      const float eff,
                          0);                                //      const float noise
    }
    kalman -> add_zplane_state("ETTL_0", -155.5);
    kalman -> add_zplane_state("ETTL_1", -158.5);
  }

  // forward hadron going direction
  if (Enable::FTTL){
    float pitch=200e-4;
    if (specialSetting.Contains("FTTLS3LC")){
      pitch=500e-4;
      for (int i = 0; i < 3; i++){
        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 287);
      kalman -> add_zplane_state("FTTL_1", 289);
      kalman -> add_zplane_state("FTTL_2", 340);
    } else if (specialSetting.Contains("FTTLS2LF")){
      for (int i = 0; i < 4; i++){
        if (i==0 || i==2 )      pitch=200e-4;   // inner rings with higher granualrity
        else           pitch=500e-4;

        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 289);
      kalman -> add_zplane_state("FTTL_1", 289);
      kalman -> add_zplane_state("FTTL_2", 340);
      kalman -> add_zplane_state("FTTL_3", 340);
     
    
    } else if (specialSetting.Contains("FTTLS2LC")){
      pitch=500e-4;
      for (int i = 0; i < 2; i++){
        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 289);
      kalman -> add_zplane_state("FTTL_1", 340);
    } else if (specialSetting.Contains("FTTLSE2LF")){
      for (int i = 0; i < 4; i++){
        if (i==0 || i==2 )      pitch=200e-4;   // inner rings with higher granualrity
        else           pitch=500e-4;

        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 287);
      kalman -> add_zplane_state("FTTL_1", 287);
      kalman -> add_zplane_state("FTTL_2", 289);
      kalman -> add_zplane_state("FTTL_3", 289);

    } else if (specialSetting.Contains("FTTLSE2LC")){
      pitch=500e-4;
      for (int i = 0; i < 2; i++){
        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 287);
      kalman -> add_zplane_state("FTTL_1", 289);
   } else {
      for (int i = 0; i < 6; i++){
        if (i==0 || i==2 || i ==4)      pitch=200e-4;   // inner rings with higher granualrity
        else           pitch=500e-4;

        kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                            PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                            pitch/sqrt(12),                    //      const float radres,
                            pitch/sqrt(12),                    //      const float phires,
                            999.,                              //      const float lonres, *ignored in plane detector*
                            0.95,                              //      const float eff,
                            0);                                //      const float noise
      }
    
      // add tof-like projection for FST layer 5 at z = 280 cm
      kalman -> add_zplane_state("FTTL_0", 287);
      kalman -> add_zplane_state("FTTL_1", 287);
      kalman -> add_zplane_state("FTTL_2", 289);
      kalman -> add_zplane_state("FTTL_3", 289);
      kalman -> add_zplane_state("FTTL_4", 340);
      kalman -> add_zplane_state("FTTL_5", 340);
    }  
  }

  //-------------------------
  // MVTX (sPHENIX) with better resol
  //-------------------------
  if (Enable::MVTX){
    //   MAPS
    kalman->add_phg4hits(
        "G4HIT_MVTX",                //      const std::string& phg4hitsNames,
        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
        5e-4,                        //      const float radres,
        5e-4,                        //      const float phires,
        5e-4,                        //      const float lonres,
        1,                           //      const float eff,
        0                            //      const float noise
    );
  }
  //-------------------------
  // TPC
  //-------------------------
  if (Enable::TPC){
    kalman->add_phg4hits(
        "G4HIT_TPC",                 //      const std::string& phg4hitsNames,
        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
        1,                           //      const float radres,
        200e-4,                      //      const float phires,
        500e-4,                      //      const float lonres,
        1,                           //      const float eff,
        0                            //      const float noise
    );
  }
  
  //-------------------------
  // EGEM
  //-------------------------
  if (Enable::EGEM){
    // GEM, 70um azimuthal resolution, 1cm radial strips
    int minEGEM = 0;
    if (!Enable::EGEM_FULL) // all to disable inner most discs
      minEGEM = 2;
    for (int i = minEGEM; i < 4; i++){
      kalman->add_phg4hits(
          Form("G4HIT_EGEM_%d", i),          //      const std::string& phg4hitsNames,
          PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
          1. / sqrt(12.),                    //      const float radres,
          70e-4,                             //      const float phires,
          100e-4,                            //      const float lonres,
          1,                                 //      const float eff,
          0                                  //      const float noise
      );
    }
  }
  //-------------------------
  // FGEM
  //-------------------------
  if (Enable::FGEM) {
    // GEM2, 70um azimuthal resolution, 1cm radial strips
    for (int i = 2; i < 5; i++) {
      kalman->add_phg4hits(Form("G4HIT_FGEM_%d", i),          //      const std::string& phg4hitsNames,
                           PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                           1. / sqrt(12.),                    //      const float radres,
                           70e-4,                             //      const float phires,
                           100e-4,                            //      const float lonres,
                           1,                                 //      const float eff,
                           0);                                //      const float noise
    }
  }
  //-------------------------
  // FEMC
  //-------------------------
  // Saved track states (projections)
  if (Enable::FEMC && G4TRACKING::PROJECTION_FEMC){
    kalman->add_state_name("FEMC");
  }

  //-------------------------
  // FHCAL
  //-------------------------
  if (Enable::FHCAL && G4TRACKING::PROJECTION_FHCAL) {
    kalman->add_state_name("FHCAL");
  }
  //-------------------------
  // CEMC
  //-------------------------

  if (Enable::CEMC && G4TRACKING::PROJECTION_CEMC){
    kalman->add_state_name("CEMC");
  }
  se->registerSubsystem(kalman);

  return;
}

//-----------------------------------------------------------------------------//

void Tracking_Eval(const std::string &outputfile, TString specialSetting = "")
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();

  //----------------
  // Fast Tracking evaluation
  //----------------
  PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
  fast_sim_eval->set_trackmapname(TRACKING::TrackNodeName);
  
  // save tracking planes in forward direction for single layer studies
  if (Enable::TRACKING_EVAL_DETAILED){
    if (Enable::FST){
      for (Int_t i = 0; i < 5; i++){
        fast_sim_eval->AddProjection(Form("FST_%d",i));
      }
    } 
    if (Enable::ALLSILICON){
      for (int i = 20; i < 25; i++) {
        fast_sim_eval->AddProjection(Form("LBLVTX_FORWARD_%d", i));
      }
    }
  }
  
  // create projections on timing layers to read t and spacial coordinates
  if (Enable::FTTL){
    if (specialSetting.Contains("FTTLS3LC")){
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
      fast_sim_eval->AddProjection("FTTL_2");
    } else if (specialSetting.Contains("FTTLS2LF")){
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
      fast_sim_eval->AddProjection("FTTL_2");
      fast_sim_eval->AddProjection("FTTL_3");
    } else if (specialSetting.Contains("FTTLS2LC")){
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
    } else if (specialSetting.Contains("FTTLSE2LF")){
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
      fast_sim_eval->AddProjection("FTTL_2");
      fast_sim_eval->AddProjection("FTTL_3");
    } else if (specialSetting.Contains("FTTLSE2LC")){
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
   } else {
      fast_sim_eval->AddProjection("FTTL_0");
      fast_sim_eval->AddProjection("FTTL_1");
      fast_sim_eval->AddProjection("FTTL_2");
      fast_sim_eval->AddProjection("FTTL_3");
      fast_sim_eval->AddProjection("FTTL_4");
      fast_sim_eval->AddProjection("FTTL_5");
    }  
  }
  if (Enable::ETTL){
    fast_sim_eval->AddProjection("ETTL_0");
    fast_sim_eval->AddProjection("ETTL_1");
  }
  if (Enable::CTTL){
    fast_sim_eval->AddProjection("CTTL_0");
    fast_sim_eval->AddProjection("CTTL_1");
  }
  
  // write to output file
  fast_sim_eval->set_filename(outputfile);
  se->registerSubsystem(fast_sim_eval);
}
#endif
