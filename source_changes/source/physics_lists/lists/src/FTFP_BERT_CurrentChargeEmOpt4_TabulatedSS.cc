#include "FTFP_BERT_CurrentChargeEmOpt4_TabulatedSS.hh"

#include "G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS.hh"

FTFP_BERT_CurrentChargeEmOpt4_TabulatedSS::FTFP_BERT_CurrentChargeEmOpt4_TabulatedSS(
  G4int ver)
  : FTFP_BERT(ver)
{
  ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(ver));
}

FTFP_BERT_CurrentChargeEmOpt4_TabulatedSS::FTFP_BERT_CurrentChargeEmOpt4_TabulatedSS(
  const G4CCPhysicsOptions& options, G4int ver)
  : FTFP_BERT(ver)
{
  ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge_TabulatedSS(options, ver));
}
