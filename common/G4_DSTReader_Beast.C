#ifndef MACRO_G4DSTREADERBEAST_C
#define MACRO_G4DSTREADERBEAST_C

#include <GlobalVariables.C>

#include <G4_Magnet_Beast.C>
#include <G4_Pipe_EIC.C>

#include <g4eval/PHG4DSTReader.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libg4eval.so)

//////////////////////////////////////////////////////////////////
/*!
  \file G4_DSTReader.C
  \brief Convert DST to human command readable TTree for quick poke around the outputs
  \author  Jin Huang
  \version $Revision:  $
  \date    $Date: $
*/
//////////////////////////////////////////////////////////////////
namespace Enable
{
  bool DSTREADER = false;
  int DSTREADER_VERBOSITY = 0;
}  // namespace Enable

namespace G4DSTREADER
{
  bool save_g4_raw = true;
}  // namespace G4DSTREADER

void G4DSTreader(const string &outputFile = "G4JLeic.root")
{
  //! debug output on screen?
  int verbosity = max(Enable::VERBOSITY, Enable::DSTREADER_VERBOSITY);

  // save a comprehensive  evaluation file
  PHG4DSTReader *ana = new PHG4DSTReader(outputFile);
  ana->set_save_particle(true);
  ana->set_load_all_particle(false);
  ana->set_load_active_particle(true);
  ana->set_save_vertex(true);

  ana->Verbosity(verbosity);

  if (G4DSTREADER::save_g4_raw)
  {
    if (Enable::PIPE)
    {
      if (Enable::ABSORBER || Enable::PIPE_ABSORBER)
        ana->AddNode("PIPE");
    }
    if (Enable::MAGNET)
    {
      if (Enable::ABSORBER || Enable::MAGNET_ABSORBER)
        ana->AddNode("MAGNET");
    }
    if (Enable::BLACKHOLE)
    {
      ana->AddNode("BH_1");
      ana->AddNode("BH_FORWARD_PLUS");
      ana->AddNode("BH_FORWARD_NEG");
    }
  }
  Fun4AllServer *se = Fun4AllServer::instance();
  se->registerSubsystem(ana);
}
#endif  // MACRO_G4DSTREADERJLEIC_C
