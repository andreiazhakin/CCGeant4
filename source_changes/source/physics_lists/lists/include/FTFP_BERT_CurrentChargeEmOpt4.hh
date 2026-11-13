#ifndef FTFP_BERT_CurrentChargeEmOpt4_hh
#define FTFP_BERT_CurrentChargeEmOpt4_hh 1

#include "FTFP_BERT.hh"
#include "G4CCPhysicsOptions.hh"

class FTFP_BERT_CurrentChargeEmOpt4 : public FTFP_BERT
{
 public:
  explicit FTFP_BERT_CurrentChargeEmOpt4(G4int ver = 1);
  explicit FTFP_BERT_CurrentChargeEmOpt4(const G4CCPhysicsOptions& options, G4int ver = 1);
  ~FTFP_BERT_CurrentChargeEmOpt4() override = default;
};

#endif
