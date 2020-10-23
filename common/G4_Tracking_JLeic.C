#ifndef MACRO_G4TRACKINGJLEIC_C
#define MACRO_G4TRACKINGJLEIC_C

#include <GlobalVariables.C>

#include <G4_VTX_JLeic.C>

#include <fun4all/Fun4AllServer.h>

#include <g4eval/SvtxEvaluator.h>

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

namespace Enable
{
  bool TRACKING = false;
  bool TRACKING_EVAL = false;
  int TRACKING_VERBOSITY = 0;
}  // namespace Enable

namespace G4TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_JLDIRC = false;
}  // namespace G4TRACKING

//-----------------------------------------------------------------------------//
void TrackingInit()
{
  TRACKING::TrackNodeName = "TrackMap";
}

void Tracking_Reco()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
  kalman->Verbosity(verbosity);

  if (G4TRACKING::DISPLACED_VERTEX)
  {
    // do not use truth vertex in the track fitting,
    // which would lead to worse momentum resolution for prompt tracks
    // but this allows displaced track analysis including DCA and vertex finding
    kalman->set_use_vertex_in_fitting(false);
    kalman->set_vertex_xy_resolution(0);  // do not smear the vertex used in the built-in DCA calculation
    kalman->set_vertex_z_resolution(0);   // do not smear the vertex used in the built-in DCA calculation
    kalman->enable_vertexing(true);       // enable vertex finding and fitting
  }
  else
  {
    // constraint to a primary vertex and use it as part of the fitting level arm
    kalman->set_use_vertex_in_fitting(true);
    kalman->set_vertex_xy_resolution(50e-4);
    kalman->set_vertex_z_resolution(50e-4);
  }

  kalman->set_sub_top_node_name("TRACKS");
  kalman->set_trackmap_out_name(TRACKING::TrackNodeName);

  //   VTX
  if (Enable::VTX)
  {
    kalman->add_phg4hits(
        "G4HIT_JLVTX",               //      const std::string& phg4hitsNames,
        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
        20e-4,                       //      const float radres,
        20e-4,                       //      const float phires,
        20e-4,                       //      const float lonres,
        1,                           //      const float eff,
        0                            //      const float noise
    );
  }
  //   CTD
  if (Enable::CTD)
  {
    kalman->add_phg4hits(
        "G4HIT_JLCTD",               //      const std::string& phg4hitsNames,
        PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
        5e-2,                        //      const float radres,
        5e-2,                        //      const float phires,
        5e-2,                        //      const float lonres,
        1,                           //      const float eff,
        0                            //      const float noise
    );
  }
  //
  // GEM0, 70um azimuthal resolution, 1cm radial strips
  if (Enable::GEM)
  {
    kalman->add_phg4hits(
        "G4HIT_GEMHADRON",                 //      const std::string& phg4hitsNames,
        PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
        1. / sqrt(12.),                    //      const float radres,
        70e-4,                             //      const float phires,
        100e-4,                            //      const float lonres,
        1,                                 //      const float eff,
        0                                  //      const float noise
    );
  }

  se->registerSubsystem(kalman);

  return;
}

void Tracking_Eval(std::string const &outputfile)
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
  fast_sim_eval->set_filename(outputfile);
  se->registerSubsystem(fast_sim_eval);
}

#endif  //  MACRO_G4TRACKINGJLEIC_C
