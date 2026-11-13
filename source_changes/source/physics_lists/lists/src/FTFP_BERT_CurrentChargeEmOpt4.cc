#include "FTFP_BERT_CurrentChargeEmOpt4.hh"

#include "G4EmStandardPhysics_option4_CurrentCharge.hh"

FTFP_BERT_CurrentChargeEmOpt4::FTFP_BERT_CurrentChargeEmOpt4(G4int ver)
  : FTFP_BERT(ver)
{
  ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(ver));
}

FTFP_BERT_CurrentChargeEmOpt4::FTFP_BERT_CurrentChargeEmOpt4(
  const G4CCPhysicsOptions& options, G4int ver)
  : FTFP_BERT(ver)
{
  ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(options, ver));
}
