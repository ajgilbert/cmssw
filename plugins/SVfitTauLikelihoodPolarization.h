#ifndef TauAnalysis_CandidateTools_SVfitTauLikelihoodPolarization_h
#define TauAnalysis_CandidateTools_SVfitTauLikelihoodPolarization_h

/** \class SVfitTauLikelihoodPolarization
 *
 * Plugin for computing likelihood for tau lepton decay "leg"
 * to be compatible with decay tau- --> X nu of polarized tau lepton into hadrons,
 * assuming  matrix element of V-A electroweak decay
 * 
 * NOTE: the system of hadrons X can either be pi-, rho- --> pi- pi0, 
 *       a1- --> pi- pi0 pi0 or a1- --> pi- pi+ pi-;
 *       tau decays into pi- pi+ pi- pi0 are **not** supported
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: SVfitTauLikelihoodPolarization.h,v 1.3 2010/09/08 10:20:01 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodPolarizationBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitVMlineShapeIntegral.h"

#include <TFormula.h>

class SVfitTauLikelihoodPolarization : public SVfitLegLikelihoodPolarizationBase<pat::Tau>
{
 public:
  SVfitTauLikelihoodPolarization(const edm::ParameterSet&);
  ~SVfitTauLikelihoodPolarization();

  void beginCandidate(const pat::Tau&);

  bool isFittedParameter(int, int) const;

 private:
  double negLogLikelihoodPolarized(const pat::Tau&, const SVfitLegSolution&, double) const;

  double probOneProngZeroPi0s(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOneProngOnePi0(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOneProngTwoPi0s(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probThreeProngZeroPi0s(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOtherDecayMode(const pat::Tau&, const SVfitLegSolution&, double) const;

//--- auxiliary functions needed for computation of likelihood 
//    for tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decays 
  double compVMa1x(double, double, double, double, double) const;
  double compVMa1DecayProbL(double, double, double, double, double) const;
  double compVMa1DecayProbT(double, double, double, double, double, double) const;

  std::vector<int> supportedTauDecayModes_;
  size_t numSupportedTauDecayModes_;

  mutable TVectorD vRec_;
  TMatrixD mapRecToGenTauDecayModes_;
  mutable TVectorD vGen_;
  mutable TVectorD vProb_;

  struct decayModeEntryType
  {
    decayModeEntryType(const edm::ParameterSet&);
    ~decayModeEntryType();
    TFormula* xSigma_;
    TFormula* xBias_;
    double pMin_;
  };

  decayModeEntryType* rhoParameters_;
  decayModeEntryType* a1NeutralParameters_;
  decayModeEntryType* a1ChargedParameters_;

  SVfitVMlineShapeIntegral* rhoLpolLineShape_;
  SVfitVMlineShapeIntegral* rhoTpolLineShape_;
  SVfitVMlineShapeIntegral* a1LpolLineShape_;
  SVfitVMlineShapeIntegral* a1TpolLineShape_;

  int tauDecayModeOther_index_;

  SVfitLegLikelihoodBase<pat::Tau>* likelihoodPhaseSpace_;

  bool useCollApproxFormulas_;

//--- temporary variables to speed-up computations
//    (computed once in constructor)
  double a1posMassTerm_;
  double a1posMassTerm2_;
  double a1negMassTerm_;
  double a1q_;
  double a1_8piDiv9_;
  double a1_16piDiv9_;
};

#endif
