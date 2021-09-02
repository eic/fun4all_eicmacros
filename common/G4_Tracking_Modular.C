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
  bool TRACKING_INNER = false;
  int TRACKING_VERBOSITY = 0;
}  // namespace Enable

namespace G4TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_EEMC = false;
  bool PROJECTION_CEMC = false;
  bool PROJECTION_FEMC = false;
  bool PROJECTION_FHCAL = false;
  bool PROJECTION_EHCAL = false;
  bool PROJECTION_DRCALO = false;
}  // namespace G4TRACKING

namespace TRACKING
{
  std::string TrackNodeNameInner = "TrackMapInner";
} // namespace TRACKING

//-----------------------------------------------------------------------------//
void TrackingInit()
{
  TRACKING::TrackNodeName = "TrackMap";
  TRACKING::TrackNodeNameInner = "TrackMapInner";
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
  kalman->Verbosity(verbosity);
  //  kalman->Smearing(false);

  // Store a list of FastSim objects. This way, we can uniformly access multiple objects,
  // setting the same objets on each when appropriate.
  std::vector<PHG4TrackFastSim *> kalmanTrackers = {kalman};

  // Option to perform an additional, separate tracking using only the inner detectors.
  PHG4TrackFastSim *kalmanInnerTracking = nullptr;
  if (Enable::TRACKING_INNER) {
    kalmanInnerTracking = new PHG4TrackFastSim("PHG4TrackFastSimInnerTracking");
    kalmanInnerTracking->Verbosity(verbosity);
    //  kalmanInnerTracking->Smearing(false);
    // Add to the list of trackers so we can uniformly access it as needed.
    kalmanTrackers.push_back(kalmanInnerTracking);
  }
  for (auto k : kalmanTrackers) {
    if (G4TRACKING::DISPLACED_VERTEX){
      // do not use truth vertex in the track fitting,
      // which would lead to worse momentum resolution for prompt tracks
      // but this allows displaced track analysis including DCA and vertex finding
      k->set_use_vertex_in_fitting(false);
      k->set_vertex_xy_resolution(0);  // do not smear the vertex used in the built-in DCA calculation
      k->set_vertex_z_resolution(0);   // do not smear the vertex used in the built-in DCA calculation
      k->enable_vertexing(true);       // enable vertex finding and fitting
    } else {
      // constraint to a primary vertex and use it as part of the fitting level arm
      k->set_use_vertex_in_fitting(true);
      k->set_vertex_xy_resolution(50e-4);
      k->set_vertex_z_resolution(50e-4);
    }
  }

  kalman->set_sub_top_node_name("TRACKS");
  kalman->set_trackmap_out_name(TRACKING::TrackNodeName);
  if (Enable::TRACKING_INNER) {
    kalmanInnerTracking->set_sub_top_node_name("TRACKS");
    kalmanInnerTracking->set_trackmap_out_name(TRACKING::TrackNodeNameInner);
  }

  //-------------------------
  // Barrel upgrade (LANL)
  //-------------------------
  // Different Barrel versions documented in arXiv:2009.0288
  if (Enable::BARREL){
    int nLayer            = 5;
    if (G4BARREL::SETTING::BARRELV6)
      nLayer            = 4;
    float_t pitch         = 20e-4;                    // default pitch size
    double r[6]           = { 3.64, 4.81, 5.98, 16.0, 22.0, -1};  //cm
    if (specialSetting.Contains("BARRELV4")){
      nLayer  = 6;
      r[3]    = 9.2;
      r[4]    = 17.;
      r[5]    = 27.;
    }
  
    for (Int_t i = 0; i < nLayer; i++){
      for (auto k : kalmanTrackers) {
        k->add_phg4hits(Form("G4HIT_BARREL_%d", i),              //      const std::string& phg4hitsNames,
                       PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                       999.,                        //      const float radres, (not used in cylindrical geom)
                       pitch/sqrt(12),              //      const float phires,
                       pitch/sqrt(12),              //      const float lonres,
                       1,                           //      const float eff,
                       0);                          //      const float noise
        if (Enable::TRACKING_EVAL_DETAILED){
          k->add_cylinder_state(Form("BARREL_%d", i), r[i]);
        }
      }
    }
  }
  
  //-------------------------
  // FST (LANL)
  //-------------------------
  // Different FST versions documented in arXiv:2009.0288
  if (Enable::FST) {
    // FORWARD DISKS
    Int_t nDisks        = 5;  
    float zFWDdisks[6]  = {35, 53, 77, 101, 125, 270};
    if (specialSetting.Contains("FSTV4"))
      nDisks        = 6;
    if (specialSetting.Contains("FSTV2"))
      zFWDdisks[4]  = 270;
    Float_t pitch       = 20e-4;          // default pitch size
    Int_t llargerPitch  = -1;             // layer number above which pitch size is higher, if -1 all initialized with small pitch
    if (specialSetting.Contains("FSTV3") || specialSetting.Contains("FSTV41"))        
      llargerPitch  = 3;
    if (specialSetting.Contains("FSTV42"))        
      llargerPitch  = 4;
    
    // intializing hits for different layers
    for (int i = 0; i < nDisks; i++) {
      if (llargerPitch != -1 && i >= llargerPitch)
        pitch           = 36.4e-4;
      for (auto k : kalmanTrackers) {
        k->add_phg4hits(Form("G4HIT_FST_%d", i),           //      const std::string& phg4hitsNames,
                       PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                       pitch/sqrt(12),                    //      const float radres,
                       pitch/sqrt(12),                    //      const float phires,
                       999.,                              //      const float lonres, number not used in vertical plane geometry
                       1,                                 //      const float eff,
                       0);                                //      const float noise
        if (Enable::TRACKING_EVAL_DETAILED){
          k->add_zplane_state(Form("FST_%d", i), zFWDdisks[i]);
        }
      }
    }
  }
  
  //-------------------------
  // Barrel ALLSILICON version (LBL)
  //-------------------------
  char nodename[100];
  if (Enable::ALLSILICON){
    // CENTRAL BARREL
    float_t pitch = 10e-4;
    float rBarrel[6]  = {3.3, 5.7, 21.0, 22.68, 39.30, 43.23};
    for (int i = 10; i < 16; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_CENTRAL_%d", i);

      for (auto k : kalmanTrackers) {
        k->add_phg4hits(
            nodename,                    // const std::string& phg4hitsNames
            PHG4TrackFastSim::Cylinder,  // const DETECTOR_TYPE phg4dettype
            999.,                        // radial-resolution [cm] (this number is not used in cylindrical geometry)
            pitch/sqrt(12),                      // azimuthal (arc-length) resolution [cm]
            pitch/sqrt(12),                      // longitudinal (z) resolution [cm]
            1,                           // efficiency (fraction)
            0                            // hit noise
        );
        if (Enable::TRACKING_EVAL_DETAILED){
          k->add_cylinder_state(Form("LBLVTX_CENTRAL_%d", i), rBarrel[i-10]);
        }
      }
    }
    // FORWARD DISKS
    float zFWDdisks[5]  = {25, 49, 73, 97, 121};
    for (int i = 20; i < 25; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_FORWARD_%d", i);
      for (auto k : kalmanTrackers) {
        k->add_phg4hits(
            nodename,                          // const std::string& phg4hitsNames
            PHG4TrackFastSim::Vertical_Plane,  // const DETECTOR_TYPE phg4dettype
            pitch/sqrt(12),                            // radial-resolution [cm]
            pitch/sqrt(12),                            // azimuthal (arc-length) resolution [cm]
            999.,                              // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
            1,                                 // efficiency (fraction)
            0                                  // hit noise
        );
        if (Enable::TRACKING_EVAL_DETAILED){
          k->add_zplane_state(Form("LBLVTX_FORWARD_%d", i), zFWDdisks[i-20]);
        }
      }
    }
    
    // BACKWARD DISKS
    float zBWDdisks[5]  = {-25, -49, -73, -97, -121};
    for (int i = 30; i < 35; i++) {
      sprintf(nodename, "G4HIT_LBLVTX_BACKWARD_%d", i);
      for (auto k : kalmanTrackers) {
        k->add_phg4hits(
            nodename,                          // const std::string& phg4hitsNames
            PHG4TrackFastSim::Vertical_Plane,  // const DETECTOR_TYPE phg4dettype
            pitch/sqrt(12),                            // radial-resolution [cm]
            pitch/sqrt(12),                            // azimuthal (arc-length) resolution [cm]
            999.,                              // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
            1,                                 // efficiency (fraction)
            0                                  // hit noise
        );
        if (Enable::TRACKING_EVAL_DETAILED){
          k->add_zplane_state(Form("LBLVTX_BACKWARD_%d", i), zBWDdisks[i-30]);
        }
      }
    }
  }
  
  //-------------------------
  // Timing tracking Layer (TTL - ORNL/Rice)
  //-------------------------
  // position resol improvement 
  float posResImp = sqrt(12);
  float pitchLGAD = 500e-4;
  if (G4TTL::SETTING::optionGran == 2)
    posResImp = sqrt(256);
  else if (G4TTL::SETTING::optionGran == 3)
    pitchLGAD   = 1300e-4;
  float resLGAD   = pitchLGAD/posResImp; 
  
  // central barrel 
  if (Enable::CTTL){
    for (int i = 0; i < G4TTL::layer[1]; i++){
      kalman->add_phg4hits(Form("G4HIT_CTTL_%d", i),           //      const std::string& phg4hitsNames,
                          PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
                          999,                               //      const float radres,
                          resLGAD,                    //      const float phires,
                          resLGAD,                    //      const float lonres, *ignored in plane detector*
                          0.95,                              //      const float eff,
                          0);                                //      const float noise
      kalman -> add_cylinder_state(Form("CTTL_%d",i), G4TTL::positionToVtx[1][i]);
    }    
  }
  
  // electron going direction
  if (Enable::ETTL){
    for (int i = 0; i < G4TTL::layer[0]; i++){
      kalman->add_phg4hits(Form("G4HIT_ETTL_%d", i),           //      const std::string& phg4hitsNames,
                          PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                          resLGAD,                    //      const float radres,
                          resLGAD,                    //      const float phires,
                          999.,                              //      const float lonres, *ignored in plane detector*
                          0.95,                              //      const float eff,
                          0);                                //      const float noise
      kalman -> add_zplane_state(Form("ETTL_%d",i), G4TTL::positionToVtx[0][i]);
    }  
  }

  // forward hadron going direction
  if (Enable::FTTL){
    for (int i = 0; i < G4TTL::layer[2]; i++){
      kalman->add_phg4hits(Form("G4HIT_FTTL_%d", i),           //      const std::string& phg4hitsNames,
                          PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                          resLGAD,                    //      const float radres,
                          resLGAD,                    //      const float phires,
                          999.,                              //      const float lonres, *ignored in plane detector*
                          0.95,                              //      const float eff,
                          0);                                //      const float noise
      kalman -> add_zplane_state(Form("FTTL_%d",i), G4TTL::positionToVtx[2][i]);
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
    int minFGEM = 0;
    if (Enable::FGEM_ORIG){
      minFGEM = 0;
    } else {
      minFGEM = 2;
    }
    // GEM2, 70um azimuthal resolution, 1cm radial strips
    for (int i = minFGEM; i < 5; i++) {
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
  if (Enable::FEMC){
    kalman->add_state_name("FEMC");
  //   kalman->add_zplane_state("FEMC", 310);
  }

  //-------------------------
  // DRCALO
  //-------------------------
  if (Enable::DRCALO) {
    kalman -> add_zplane_state("DRCALO_0", 300);
  }

  //-------------------------
  // FHCAL
  //-------------------------
  if (Enable::FHCAL) {
    kalman->add_state_name("FHCAL");
  //   kalman->add_zplane_state("FHCAL", 350);
  }

  //-------------------------
  // LFHCAL
  //-------------------------
  if (Enable::LFHCAL) {
    // kalman->add_state_name("FHCAL");
    kalman->add_zplane_state("LFHCAL", 350);
  }
  //-------------------------
  // CEMC
  //-------------------------

  if (Enable::CEMC){
    kalman->add_state_name("CEMC");
  }
  
  //-------------------------
  // BECAL
  //-------------------------

  if (Enable::BECAL){
    // kalman->add_state_name("BECAL");
    kalman -> add_cylinder_state("BECAL", 105);
  }
  
  //-------------------------
  // HCALIN
  //-------------------------

  if (Enable::HCALIN){
    kalman->add_state_name("HCALIN");
  }
  
  
  //-------------------------
  // HCALOUT
  //-------------------------

  if (Enable::HCALOUT){
    kalman->add_state_name("HCALOUT");
  }
  
  //-------------------------
  // EEMC
  //-------------------------
  if ((Enable::EEMC || Enable::EEMCH ) )
  {
    kalman->add_state_name("EEMC");
  }
  
  
  //-------------------------
  // EEMC
  //-------------------------
  if ((Enable::EEMCH ) )
  {
    // kalman->add_state_name("EEMCG");
    kalman->add_zplane_state("EEMCG", -210);
  }
  
  //-------------------------
  // EEMC
  //-------------------------
  if ((Enable::EHCAL ) )
  {
    kalman->add_state_name("EHCAL");
    // kalman->add_zplane_state("EHCAL", -350);
  }
  
  //-------------------------
  // DIRC
  //-------------------------
  if ((Enable::DIRC ) )
  {
    // kalman->add_state_name("DIRC");
    kalman->add_zplane_state("DIRC", 185);
  }


  for (auto k : kalmanTrackers) {
    se->registerSubsystem(k);
  }
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
      Int_t nDisks        = 5;  
      if (specialSetting.Contains("FSTV4"))
        nDisks        = 6;
      for (Int_t i = 0; i < nDisks; i++){
        fast_sim_eval->AddProjection(Form("FST_%d",i));
      }
    } 
    if (Enable::BARREL){    
      Int_t nLayer  = 5;
      if (specialSetting.Contains("BARRELV4"))
        nLayer  = 6;

      for (Int_t i = 0; i < nLayer; i++){
        fast_sim_eval->AddProjection(Form("BARREL_%d",i));
      }
    }
    
    if (Enable::ALLSILICON){
      for (int i = 10; i < 16; i++) {
        fast_sim_eval->AddProjection(Form("LBLVTX_CENTRAL_%d", i));
      }
      for (int i = 20; i < 25; i++) {
        fast_sim_eval->AddProjection(Form("LBLVTX_FORWARD_%d", i));
      }
      for (int i = 30; i < 35; i++) {
        fast_sim_eval->AddProjection(Form("LBLVTX_BACKWARD_%d", i));
      }
    }
  }
  
  // create projections on timing layers to read t and spacial coordinates
  if (Enable::FTTL){
    for (int l = 0; l < G4TTL::layer[2]; l++)
      fast_sim_eval->AddProjection(Form("FTTL_%d",l));
  }
  if (Enable::ETTL){
    for (int l = 0; l < G4TTL::layer[0]; l++)
      fast_sim_eval->AddProjection(Form("ETTL_%d",l));
  }
  if (Enable::CTTL){
    for (int l = 0; l < G4TTL::layer[1]; l++)
      fast_sim_eval->AddProjection(Form("CTTL_%d",l));
  }

  if(Enable::FHCAL && G4TRACKING::PROJECTION_FHCAL) fast_sim_eval->AddProjection("FHCAL");
  if(Enable::FEMC && G4TRACKING::PROJECTION_FEMC) fast_sim_eval->AddProjection("FEMC");
  // if(Enable::DRCALO && G4TRACKING::PROJECTION_DRCALO) fast_sim_eval->AddProjection("DRCALO_0");

  if(Enable::DIRC) fast_sim_eval->AddProjection("DIRC");

  // write to output file
  fast_sim_eval->set_filename(outputfile);
  se->registerSubsystem(fast_sim_eval);
}
#endif
