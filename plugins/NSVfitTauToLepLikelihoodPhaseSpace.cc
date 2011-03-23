#include "TauAnalysis/CandidateTools/plugins/NSVfitTauToLepLikelihoodPhaseSpace.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToLepHypothesis.h"

#include <TMath.h>

#include <limits>

using namespace SVfit_namespace;

template <typename T>
NSVfitTauToLepLikelihoodPhaseSpace<T>::NSVfitTauToLepLikelihoodPhaseSpace(const edm::ParameterSet& cfg)
  : NSVfitSingleParticleLikelihood(cfg)
{
// nothing to be done yet...
}

template <typename T>
NSVfitTauToLepLikelihoodPhaseSpace<T>::~NSVfitTauToLepLikelihoodPhaseSpace()
{
// nothing to be done yet...
}

template <typename T>
void NSVfitTauToLepLikelihoodPhaseSpace<T>::beginJob(NSVfitAlgorithmBase* algorithm)
{
  algorithm->requestFitParameter(prodParticleLabel_, kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_nuInvMass,  pluginName_);
}

template <typename T>
double NSVfitTauToLepLikelihoodPhaseSpace<T>::operator()(const NSVfitSingleParticleHypothesisBase* hypothesis) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with three-body decay,
//    assuming constant matrix element,
//    so that energy and angular distribution of decay products is solely determined by phase-space
//
//    NOTE: the parametrization of the three-body decay phase-space is taken from the PDG:
//          K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010);
//          formula 38.20b
//
  const NSVfitTauToLepHypothesis<T>* hypothesis_T = dynamic_cast<const NSVfitTauToLepHypothesis<T>*>(hypothesis);
  assert(hypothesis_T != 0);

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToLepLikelihoodPhaseSpace::operator()>:" << std::endl;

  double decayAngle = hypothesis_T->decay_angle_rf();
  double sinDecayAngle = TMath::Sin(decayAngle);
  double nuMass = hypothesis_T->p4invis_rf().mass();
  double visMass = hypothesis_T->p4vis_rf().mass();

  if ( this->verbosity_ ) {
    std::cout << " decayAngle = " << decayAngle << std::endl;
    std::cout << " nuMass = " << nuMass << std::endl;
    std::cout << " visMass = " << visMass << std::endl;
  }

  double prob = 0.5*nuMass
               *TMath::Sqrt((tauLeptonMass2 - square(nuMass + visMass))*(tauLeptonMass2 - square(nuMass - visMass)))/(2*tauLeptonMass)
               *sinDecayAngle;

  double nll = 0.;
  if ( prob > 0. ) {
    nll = -TMath::Log(prob);
  } else {
    if ( prob < 0. ) 
      edm::LogWarning ("NSVfitTauToLepLikelihoodPhaseSpace::operator()")
	<< " Unphysical solution: prob = " << prob << " --> returning very large negative number !!";
    nll = std::numeric_limits<float>::max();
  }

  if ( sinDecayAngle > 0. && nuMass > 0. && nuMass < (tauLeptonMass - visMass) ) {
    double logP1 = TMath::Log(nuMass) - TMath::Log(2.);
    double logP3 = 0.5*TMath::Log((tauLeptonMass2 - square(nuMass + visMass))*(tauLeptonMass2 - square(nuMass - visMass)))
                  - TMath::Log(2*tauLeptonMass);
    double nll2 = -(TMath::Log(sinDecayAngle) + logP1 + logP3);
    std::cout << "nll1 = " << nll << ", nll2 = " << nll2 << std::endl;
  }

  if ( this->verbosity_ ) std::cout << "--> nll = " << nll << std::endl;

  return nll;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef NSVfitTauToLepLikelihoodPhaseSpace<pat::Electron> NSVfitTauToElecLikelihoodPhaseSpace;
typedef NSVfitTauToLepLikelihoodPhaseSpace<pat::Muon> NSVfitTauToMuLikelihoodPhaseSpace;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToElecLikelihoodPhaseSpace, "NSVfitTauToElecLikelihoodPhaseSpace");
DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToMuLikelihoodPhaseSpace, "NSVfitTauToMuLikelihoodPhaseSpace");
