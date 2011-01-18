#ifndef TauAnalysis_CandidateTools_SVfitAlgorithm_h
#define TauAnalysis_CandidateTools_SVfitAlgorithm_h

/** \class SVfitAlgorithm
 *
 * Class used to reconstruct di-tau invariant mass via fitting/likelihood techniques.
 *
 * The fit is performed by passing log-likelihood functions to be minimized to Minuit;
 * the actual likelihood functions are implememted as plugins,
 * providing flexibility to reconstruct multiple di-tau invariant mass solutions
 * using different combinations of
 *  kinematic, missing transverse momentum, tau lepton lifetime,...
 * information.
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.24 $
 *
 * $Id: SVfitAlgorithm.h,v 1.24 2010/11/16 19:37:01 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitEventVertexRefitter.h"
#include "TauAnalysis/CandidateTools/interface/SVfitVertexOnTrackFinder.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// DQM Instrumentation includes
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <TFitterMinuit.h>
#include <Minuit2/FCNBase.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom3.h>

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/foreach.hpp>

// forward declaration of SVfitAlgorithm class
template<typename T1, typename T2> class SVfitAlgorithm;

namespace SVfit_namespace
{
  enum fitParameter { kPrimaryVertexShiftX, kPrimaryVertexShiftY, kPrimaryVertexShiftZ,
    kLeg1thetaRest, kLeg1phiLab, kLeg1decayDistanceLab, kLeg1nuInvMass,
    kLeg1thetaVMrho, kLeg1thetaVMa1, kLeg1thetaVMa1r, kLeg1phiVMa1r,
    kLeg2thetaRest, kLeg2phiLab, kLeg2decayDistanceLab, kLeg2nuInvMass,
    kLeg2thetaVMrho, kLeg2thetaVMa1, kLeg2thetaVMa1r, kLeg2phiVMa1r };
  enum tauDecayProducts { kLeg1, kLeg2 };

  const unsigned defaultMinuitMaxInterations = 5000;
}

template<typename T1, typename T2>
class SVfitMinuitFCNadapter : public ROOT::Minuit2::FCNBase
{
  public:
    SVfitMinuitFCNadapter()
        : svFitAlgorithm_(0)
    {}
    ~SVfitMinuitFCNadapter() {}

    void setSVfitAlgorithm(const SVfitAlgorithm<T1,T2>* svFitAlgorithm) { svFitAlgorithm_ = svFitAlgorithm; }

    /// define "objective" function called by Minuit
    double operator()(const std::vector<double>& x) const {
      return svFitAlgorithm_->negLogLikelihood(x);
    }

    /// define increase in "objective" function used by Minuit to determine one-sigma error contours;
    /// in case of (negative) log-likelihood functions, the value needs to be 0.5
    double Up() const { return 0.5; }
  private:
    const SVfitAlgorithm<T1,T2>* svFitAlgorithm_;
};

template<typename T1, typename T2>
class SVfitAlgorithm
{
  public:
    SVfitAlgorithm(const edm::ParameterSet& cfg) :currentDiTau_(0) {
      verbosity_ = cfg.exists("verbosity") ? cfg.getParameter<int>("verbosity") : 0;
      // Check if we want to use the DQM monitors
      dqmStore_ = NULL;
      if (cfg.exists("monitor")) {
        std::cout << "Enabling DQM monitoring of SV fit" << std::endl;
        monitorCfg_ = cfg.getParameter<edm::ParameterSet>("monitor");
        // Load DQM store
        if (edm::Service<DQMStore>().isAvailable()) {
          dqmStore_ = &(*edm::Service<DQMStore>());
          enableLikelihoodMonitoring_ = true;
        } else {
          edm::LogError ("SVfitAlgorithm") <<
            " Failed to access dqmStore --> histograms will NOT be booked !!";
        }
      }
      name_ = cfg.getParameter<std::string>("name");
      if ( verbosity_ ) std::cout << "initializing SVfit algorithm with name = " << name_ << std::endl;

      eventVertexRefitAlgorithm_ = new SVfitEventVertexRefitter(cfg);

      likelihoodsSupportPolarization_ = false;

      typedef std::vector<edm::ParameterSet> vParameterSet;
      vParameterSet cfgLikelihoodFunctions = cfg.getParameter<vParameterSet>("likelihoodFunctions");
      for ( vParameterSet::const_iterator cfgLikelihoodFunction = cfgLikelihoodFunctions.begin();
           cfgLikelihoodFunction != cfgLikelihoodFunctions.end(); ++cfgLikelihoodFunction ) {
        std::string pluginType = cfgLikelihoodFunction->getParameter<std::string>("pluginType");
        typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<T1,T2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
        SVfitDiTauLikelihoodBase<T1,T2>* likelihoodFunction
            = SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgLikelihoodFunction);
        likelihoodsSupportPolarization_ |= likelihoodFunction->supportsPolarization();
        likelihoodFunctions_.push_back(likelihoodFunction);
      }

      minuitMaxInterations_ = cfg.exists("minuitMaxInterations") ?
          cfg.getParameter<unsigned>("minuitMaxInterations") : SVfit_namespace::defaultMinuitMaxInterations;

//--- initialize Minuit
      minuitFCNadapter_.setSVfitAlgorithm(this);

      minuit_ = new TFitterMinuit();
      minuit_->SetMinuitFCN(&minuitFCNadapter_);
//--- set Minuit strategy = 2,
//    in order to enable reliable error estimates
//    ( cf. http://www-cdf.fnal.gov/physics/statistics/recommendations/minuit.html )
      minuit_->SetStrategy(2);
      minuit_->SetMaxIterations(minuitMaxInterations_);

      //std::cout << "<SVfitAlgorithm::SVfitAlgorithm>:" << std::endl;
      //std::cout << " disabling MINUIT output..." << std::endl;
      //minuit_->SetPrintLevel(1);
      minuit_->SetPrintLevel(-1);
      if (verbosity_)
        minuit_->SetPrintLevel(0);
      minuit_->SetErrorDef(0.5);

      minuit_->CreateMinimizer();

      minuitFittedParameterValues_.resize(minuitNumParameters_);

      if ( cfg.exists("parameterizeVertexAlongTrackLeg1") ) {
	if ( verbosity_ ) std::cout << "--> enabling parameterizeVertexAlongTrackLeg1" << std::endl;
        parameterizeVertexAlongTrackLeg1_ = cfg.getParameter<bool>("parameterizeVertexAlongTrackLeg1");
      } else {
	if ( verbosity_ ) std::cout << "--> disabling parameterizeVertexAlongTrackLeg1" << std::endl;
        parameterizeVertexAlongTrackLeg1_ = false;
      }
      if ( cfg.exists("parameterizeVertexAlongTrackLeg2") ) {
	if ( verbosity_ ) std::cout << "--> enabling parameterizeVertexAlongTrackLeg2" << std::endl;
        parameterizeVertexAlongTrackLeg2_ = cfg.getParameter<bool>("parameterizeVertexAlongTrackLeg2");
      } else {
	if ( verbosity_ ) std::cout << "--> disabling parameterizeVertexAlongTrackLeg2" << std::endl;
        parameterizeVertexAlongTrackLeg2_ = false;
      }

      if ( cfg.exists("estUncertainties") ) {
        edm::ParameterSet cfgEstUncertainties = cfg.getParameter<edm::ParameterSet>("estUncertainties");
        numSamplings_ = cfgEstUncertainties.getParameter<int>("numSamplings");

//--- make numSamplings parameter always an odd number,
//    so that the median of numSamplings random values is well defined
        if ( (numSamplings_ % 2) == 0 ) ++numSamplings_;
      }

      //print(std::cout);
    }

    ~SVfitAlgorithm() {
      delete eventVertexRefitAlgorithm_;

      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::iterator it = likelihoodFunctions_.begin();
           it != likelihoodFunctions_.end(); ++it ) {
        delete (*it);
      }
    }

    void beginJob() {
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->beginJob();
      }
    }

    void beginEvent(edm::Event& evt, const edm::EventSetup& es) {
      currentEvent_ = evt.id();
      tauPairIndex_ = 0;
      //std::cout << "<SVfitAlgorithm::beginEvent>:" << std::endl;
      SVfit::track::VertexOnTrackFinder::setEventSetup(es);

      eventVertexRefitAlgorithm_->beginEvent(evt, es);

      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->beginEvent(evt, es);
      }
    }

    void print(std::ostream& stream) const {
      stream << "<SVfitAlgorithm::print>" << std::endl;
      stream << " name = " << name_ << std::endl;
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->print(stream);
      }
      stream << " likelihoodsSupportPolarization = " << likelihoodsSupportPolarization_ << std::endl;
      stream << " minuitNumParameters = " << minuitNumParameters_ << std::endl;
      stream << " numSamplings = " << numSamplings_ << std::endl;
      stream << std::endl;
    }

    std::vector<SVfitDiTauSolution> fit(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate)
    {
      // Increment tau pair counter
      tauPairIndex_++;

      if ( verbosity_ ) {
        std::cout << "<SVfitAlgorithm::fit>:" << std::endl;
        std::cout << " name = " << name_ << std::endl;
      }

      std::vector<SVfitDiTauSolution> solutions;

//--- refit primary event vertex
//    excluding tracks associated to tau decay products
      std::vector<reco::TrackBaseRef> leg1Tracks = leg1TrackExtractor_(*diTauCandidate.leg1());
      std::vector<reco::TrackBaseRef> leg2Tracks = leg2TrackExtractor_(*diTauCandidate.leg2());
      TransientVertex pv = eventVertexRefitAlgorithm_->refit(leg1Tracks, leg2Tracks);
      if ( verbosity_ ) {
        std::cout << " refitted event vertex (#tracks = " << pv.originalTracks().size() << "):"
		  << " x = " << pv.position().x() << " +/- " << TMath::Sqrt(pv.positionError().cxx()) << ","
		  << " y = " << pv.position().y() << " +/- " << TMath::Sqrt(pv.positionError().cyy()) << ","
		  << " z = " << pv.position().z() << " +/- " << TMath::Sqrt(pv.positionError().czz()) << std::endl;
      }

//--- check if at least one likelihood function supports polarization
//    (i.e. returns a polarization depent likelihood value);
//    in case none of the likelihood functions support polarization,
//    run the fit only once (for "unknown" polarization),
//    in order to save computing time
      if ( likelihoodsSupportPolarization_ ) {
        for ( int leg1PolarizationHypothesis = SVfitLegSolution::kLeftHanded;
             leg1PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg1PolarizationHypothesis ) {
          for ( int leg2PolarizationHypothesis = SVfitLegSolution::kLeftHanded;
               leg2PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg2PolarizationHypothesis ) {
            currentDiTauSolution_ = SVfitDiTauSolution((SVfitLegSolution::polarizationHypothesisType)leg1PolarizationHypothesis,
                                                       (SVfitLegSolution::polarizationHypothesisType)leg2PolarizationHypothesis);
            // Associate the tracks
            currentDiTauSolution_.leg1_.tracks_ = leg1Tracks;
            currentDiTauSolution_.leg2_.tracks_ = leg2Tracks;
            // Refit the vertices if possible
            currentDiTauSolution_.leg1_.recoVertex_ =
	      eventVertexRefitAlgorithm_->fitSecondaryVertex(leg1Tracks);
            currentDiTauSolution_.leg2_.recoVertex_ =
	      eventVertexRefitAlgorithm_->fitSecondaryVertex(leg2Tracks);
            fitPolarizationHypothesis(diTauCandidate, currentDiTauSolution_, pv);
            solutions.push_back(currentDiTauSolution_);
          }
        }
      } else {
        currentDiTauSolution_ = SVfitDiTauSolution(SVfitLegSolution::kUnknown, SVfitLegSolution::kUnknown);
        // Associate the tracks
        currentDiTauSolution_.leg1_.tracks_ = leg1Tracks;
        currentDiTauSolution_.leg2_.tracks_ = leg2Tracks;
        // Refit the vertices if possible
        currentDiTauSolution_.leg1_.recoVertex_ =
	  eventVertexRefitAlgorithm_->fitSecondaryVertex(leg1Tracks);
        currentDiTauSolution_.leg2_.recoVertex_ =
	  eventVertexRefitAlgorithm_->fitSecondaryVertex(leg2Tracks);
        fitPolarizationHypothesis(diTauCandidate, currentDiTauSolution_, pv);
        solutions.push_back(currentDiTauSolution_);
      }

      return solutions;
    }

    double negLogLikelihood(const std::vector<double>& x) const {
      ++indexFitFunctionCall_;

      if ( verbosity_ ) {
        std::cout << "<SVfitAlgorithm::negLogLikelihood>:" << std::endl;
        std::cout << " indexFitFunctionCall = " << indexFitFunctionCall_ << std::endl;
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          // only print unfixed parameters
          if ( !minuit_->IsFixed(iParameter) )
	    std::cout << " Parameter #" << iParameter << " "
              << minuit_->GetParName(iParameter) << " = "
              << x[iParameter] << std::endl;
        }
      }

      if ( !currentDiTau_ ) {
        edm::LogError("SVfitAlgorithm::logLikelihood")
	  << " Pointer to currentDiTau has not been initialized --> skipping !!";
        return 0.;
      }


      applyParameters(currentDiTauSolution_, x);

      double negLogLikelihood = 0.;
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        // Only use this likelihood if it is enabled in this iteration.
        if ( (*likelihoodFunction)->isFitted(fitIteration_) ) {
          negLogLikelihood += (**likelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
        }
      }

      if ( verbosity_ ) std::cout << "--> negLogLikelihood = " << negLogLikelihood << ","
				  << " SVfit mass = " << currentDiTauSolution_.p4().mass() << std::endl;

//-- in order to resolve ambiguities and improve convergence of the fit,
//   add "penalty" terms in case leg1phiLab, leg2phiLab parameters
//   are outside of the "physical" interval -pi..+pi
      double penalty = 0.;
      double leg1phiLab = x[SVfit_namespace::kLeg1phiLab];
      if ( TMath::Abs(leg1phiLab) > TMath::Pi() ) penalty += SVfit_namespace::square(TMath::Abs(leg1phiLab) - TMath::Pi());
      double leg2phiLab = x[SVfit_namespace::kLeg2phiLab];
      if ( TMath::Abs(leg2phiLab) > TMath::Pi() ) penalty += SVfit_namespace::square(TMath::Abs(leg2phiLab) - TMath::Pi());

      // Check if we want to monitor the fit progress
      if (dqmStore_ && enableLikelihoodMonitoring_) {
        for (unsigned iParameter = 0; iParameter < minuitNumParameters_;
            ++iParameter ) {
          // Check if we are monitoring this parameter
          MEMap::const_iterator paramMonitorIter = monitors_.find(
              minuit_->GetParName(iParameter));
          if (paramMonitorIter != monitors_.end()) {
            MonitorElement* monitor = paramMonitorIter->second;
            if (!monitor) {
              edm::LogError("SVfitAlgorith::monitor")
                << " parameter " << iParameter << " has a null monitor element!"
                << " monitor will not be filled!";
              continue;
            }
            monitor->getTH1F()->SetBinContent(
                indexFitFunctionCall_ + 1, x[iParameter]);
          }
        }
        // Check if we want to monitor NLL & mass, and current fit iteration.
        MEMap::const_iterator nllMonitorIter = monitors_.find("nll");
        if (nllMonitorIter != monitors_.end()) {
          nllMonitorIter->second->getTH1F()->SetBinContent(
              indexFitFunctionCall_ + 1, negLogLikelihood + penalty);
        }
        MEMap::const_iterator massMonitorIter = monitors_.find("mass");
        if (massMonitorIter != monitors_.end()) {
          massMonitorIter->second->getTH1F()->SetBinContent(
              indexFitFunctionCall_ + 1, currentDiTauSolution_.p4().mass());
        }
        MEMap::const_iterator iterMonitorIter = monitors_.find("iter");
        if (iterMonitorIter != monitors_.end()) {
          iterMonitorIter->second->getTH1F()->SetBinContent(
              indexFitFunctionCall_ + 1, fitIteration_);
        }
      }
      return negLogLikelihood + penalty;
    }

  private:

    unsigned int fitIteration_;

    void fitPolarizationHypothesis(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate,
                                   SVfitDiTauSolution& solution,
                                   const TransientVertex& pv)
    {
      if ( verbosity_ ) std::cout << "<SVfitAlgorithm::fitPolarizationHypothesis>:" << std::endl;

      //--- initialize pointer to current diTau object
      currentDiTau_ = &diTauCandidate;

//--- initialize data-members of diTauSolution object
      if ( pv.isValid() ) {
        currentDiTauSolution_.eventVertexPosition_(0) = pv.position().x();
        currentDiTauSolution_.eventVertexPosition_(1) = pv.position().y();
        currentDiTauSolution_.eventVertexPosition_(2) = pv.position().z();
        currentDiTauSolution_.eventVertexPositionErr_ = pv.positionError().matrix_new();
      }
      // CV: need to add protection against case that primary event vertex is not valid <-- FIXME ?
      currentDiTauSolution_.eventVertexIsValid_ = pv.isValid();
      currentDiTauSolution_.leg1_.p4Vis_ = diTauCandidate.leg1()->p4();
      currentDiTauSolution_.leg2_.p4Vis_ = diTauCandidate.leg2()->p4();

      // Turn off vertex parameterization in case neutrals are present - leads
      // to situations with non-existent solutions.
      if (leg1NeutralActivity_(*diTauCandidate.leg1())) {
        if (verbosity_)
          std::cout << "Disabling vertex parameterization by leg 1, it has neutrals."
            << std::endl;
        parameterizeVertexAlongTrackLeg1_ = false;
      }
      if (leg2NeutralActivity_(*diTauCandidate.leg2())) {
        if (verbosity_)
          std::cout << "Disabling vertex parameterization by leg 2, it has neutrals."
            << std::endl;
        parameterizeVertexAlongTrackLeg2_ = false;
      }

//--- initialize start-values of Minuit fit parameters
//
//    CV: how to deal with measurement errors in the visible momenta of the two tau decay "legs"
//        when setting the maximum mass of the neutrino system ?
//

      double pvPositionXerr, pvPositionYerr, pvPositionZerr;
      if ( pv.isValid() ) {
        pvPositionXerr = TMath::Sqrt(pv.positionError().cxx());
        pvPositionYerr = TMath::Sqrt(pv.positionError().cyy());
        pvPositionZerr = TMath::Sqrt(pv.positionError().czz());
      } else {
        pvPositionXerr = 0.01;
        pvPositionYerr = 0.01;
        pvPositionZerr = 0.01;
      }

      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftX, "pv_dx", 0., pvPositionXerr,  -1.,  +1.);
      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftY, "pv_dy", 0., pvPositionYerr,  -1.,  +1.);
      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftZ, "pv_dz", 0., pvPositionZerr, -50., +50.);
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaRest, "sv1_thetaRest", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      double leg1PhiLabStepSize = ( parameterizeVertexAlongTrackLeg1_ ) ?
	TMath::Pi()/100. : TMath::Pi();
      minuit_->SetParameter(SVfit_namespace::kLeg1phiLab, "sv1_phiLab", 0., leg1PhiLabStepSize, 0., 0.); // do not set limits for phiLab

      double leg1decayDistanceLab0, leg1decayDistanceLabStepSize;
      if ( parameterizeVertexAlongTrackLeg1_ ) {
	leg1decayDistanceLab0 = 0.;
	leg1decayDistanceLabStepSize = 0.0100; // 100 microns
      } else {
	double gamma_cTauLifetime = (diTauCandidate.leg1()->energy()/SVfit_namespace::tauLeptonMass)*SVfit_namespace::cTauLifetime;
	if ( verbosity_ ) {
	  std::cout << "leg1: energy = " << diTauCandidate.leg1()->energy()
		    << " --> gamma = " << (diTauCandidate.leg1()->energy()/SVfit_namespace::tauLeptonMass) << std::endl;
	  std::cout << " gamma_cTauLifetime = " << gamma_cTauLifetime << std::endl;
	}

	// do not set limits for decayDistanceLab; guarantee that decayDistanceLab is positive by taking square-root as fit parameter
	leg1decayDistanceLab0 = TMath::Sqrt(gamma_cTauLifetime);
	leg1decayDistanceLabStepSize = 0.1*leg1decayDistanceLab0;
      }
      minuit_->SetParameter(SVfit_namespace::kLeg1decayDistanceLab, "sv1_decayDistanceLab",
			    leg1decayDistanceLab0, leg1decayDistanceLabStepSize, 0., 0.);

      double leg1NuMass0, leg1NuMassErr, leg1NuMassMax;
      if ( !SVfit_namespace::isMasslessNuSystem<T1>() ) {
        leg1NuMass0 = 0.8;
        leg1NuMassErr = 0.4;
        leg1NuMassMax = SVfit_namespace::tauLeptonMass - diTauCandidate.leg1()->mass();
      } else {
        leg1NuMass0 = 0.;
        leg1NuMassErr = 1.;
        leg1NuMassMax = 0.;
      }
      if ( verbosity_ ) std::cout << " leg1NuMass0 = " << leg1NuMass0 << ", leg1NuMassMax = " << leg1NuMassMax << std::endl;
      minuit_->SetParameter(SVfit_namespace::kLeg1nuInvMass, "sv1_m12", leg1NuMass0, leg1NuMassErr, 0., leg1NuMassMax);
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMrho, "sv1_thetaVMrho", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMa1, "sv1_thetaVMa1", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMa1r, "sv1_thetaVMa1r", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1phiVMa1r, "sv1_phiVMa1r", 0., TMath::Pi(), 0., 0.); // do not set limits for phiVMa1r
      minuit_->SetParameter(SVfit_namespace::kLeg2thetaRest, "sv2_thetaRest", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      double leg2PhiLabStepSize = ( parameterizeVertexAlongTrackLeg2_ ) ? TMath::Pi()/100. : TMath::Pi();
      minuit_->SetParameter(SVfit_namespace::kLeg2phiLab, "sv2_phiLab", 0., leg2PhiLabStepSize, 0., 0.); // do not set limits for phiLab

      double leg2decayDistanceLab0, leg2decayDistanceLabStepSize;
      if ( parameterizeVertexAlongTrackLeg2_ ) {
	leg2decayDistanceLab0 = 0.;
	leg2decayDistanceLabStepSize = 0.0100; // 100 microns
      } else {
	double gamma_cTauLifetime = (diTauCandidate.leg2()->energy()/SVfit_namespace::tauLeptonMass)*SVfit_namespace::cTauLifetime;
	if ( verbosity_ ) {
	  std::cout << "leg2: energy = " << diTauCandidate.leg2()->energy()
		    << " --> gamma = " << (diTauCandidate.leg2()->energy()/SVfit_namespace::tauLeptonMass) << std::endl;
	  std::cout << " gamma_cTauLifetime = " << gamma_cTauLifetime << std::endl;
	}

	// do not set limits for decayDistanceLab; guarantee that decayDistanceLab is positive by taking square-root as fit parameter
	leg2decayDistanceLab0 = TMath::Sqrt(gamma_cTauLifetime);
	leg2decayDistanceLabStepSize = 0.1*leg2decayDistanceLab0;
      }
      minuit_->SetParameter(SVfit_namespace::kLeg2decayDistanceLab, "sv2_decayDistanceLab",
			    leg2decayDistanceLab0, leg2decayDistanceLabStepSize, 0., 0.);

      double leg2NuMass0, leg2NuMassErr, leg2NuMassMax;
      if ( !SVfit_namespace::isMasslessNuSystem<T2>() ) {
        leg2NuMass0 = 0.8;
        leg2NuMassErr = 0.4;
        leg2NuMassMax = SVfit_namespace::tauLeptonMass - diTauCandidate.leg2()->mass();
      } else {
        leg2NuMass0 = 0.;
        leg2NuMassErr = 1.;
        leg2NuMassMax = 0.;
      }
      if ( verbosity_ ) std::cout << " leg2NuMass0 = " << leg2NuMass0 << ", leg2NuMassMax = " << leg2NuMassMax << std::endl;
      minuit_->SetParameter(SVfit_namespace::kLeg2nuInvMass, "sv2_m12", leg2NuMass0, leg2NuMassErr, 0., leg2NuMassMax);
      minuit_->SetParameter(SVfit_namespace::kLeg2thetaVMrho, "sv2_thetaVMrho", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg2thetaVMa1, "sv2_thetaVMa1", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg2thetaVMa1r, "sv2_thetaVMa1r", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg2phiVMa1r, "sv2_phiVMa1r", 0., TMath::Pi(), 0., 0.); // do not set limits for phiVMa1r

      unsigned int fitIterations = 0;

      // Clear out the old DQM monitors
      monitors_.clear();

      // The list of all parameters that will be fitted.
      std::set<std::string> fittedParameterNames;
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	    likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        // Load information about the candidate
        (*likelihoodFunction)->beginCandidate(diTauCandidate);
        // Figure out how many fit iterations our selected likelihoods need
        if ((*likelihoodFunction)->firstFit() > fitIterations) {
          fitIterations = (*likelihoodFunction)->firstFit();
        }
        // Get the list of all fitted parameters names (for all iterations)
        // We use this to book our monitoring histograms.
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          if ((*likelihoodFunction)->isFittedParameter(iParameter)) {
            fittedParameterNames.insert(minuit_->GetParName(iParameter));
          }
        }
      }

      // Build the basic prefix for DQM stuff for this event.
      stringstream meNameBase;
      meNameBase << name_ << "_evt_"
        << currentEvent_.event() << "_"
        << currentEvent_.luminosityBlock() << "_"
        << currentEvent_.run() << "_"
        << tauPairIndex_  << "_";
      // Book monitoring histograms for all of our variables
      if (dqmStore_) {
        dqmStore_->setCurrentFolder(
            monitorCfg_.getParameter<std::string>("dqmDirectory"));
        unsigned int maxFunctionCalls = monitorCfg_.getParameter<unsigned int>(
            "maxFunctionCalls");
        BOOST_FOREACH(const std::string& parName, fittedParameterNames) {
          // Build the name
          std::string meName = meNameBase.str() + parName;
          monitors_[parName] = dqmStore_->book1D(
              meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
          // Check if we want to do scans of each variable at the end of the run
        }
        // Add NLL and mass tracking
        std::string meName = meNameBase.str() + "mass";
        monitors_["mass"] = dqmStore_->book1D(
            meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
        meName = meNameBase.str() + "nll";
        monitors_["nll"] = dqmStore_->book1D(
            meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
        // Track fit iteration & minuit status
        meName = meNameBase.str() + "iter";
        monitors_["iter"] = dqmStore_->book1D(
            meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
        meName = meNameBase.str() + "mstatus";
        monitors_["mstatus"] = dqmStore_->book1D(
            meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
        meName = meNameBase.str() + "edm";
        monitors_["edm"] = dqmStore_->book1D(
            meName, meName, maxFunctionCalls, 0, maxFunctionCalls);
      }

      minuitStatus_ = -20;
      minuitEdm_ = 0;
      // Loop over the desired number of fit iterations.
      indexFitFunctionCall_ = 0;
      for ( fitIteration_ = 0; fitIteration_ < (fitIterations + 1); ++fitIteration_ ) {
//--- lock (i.e. set to fixed values) Minuit parameters
//    which are not constrained by any likelihood function (in this iteration)
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          bool minuitLockParameter = true;
          for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
		likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
            if ( (*likelihoodFunction)->isFittedParameter(iParameter) &&
		 (*likelihoodFunction)->isFitted(fitIteration_) ) {
              minuitLockParameter = false;
              break;
            }
          }

          if (  minuitLockParameter && !minuit_->IsFixed(iParameter) ) minuit_->FixParameter(iParameter);
          if ( !minuitLockParameter &&  minuit_->IsFixed(iParameter) ) minuit_->ReleaseParameter(iParameter);

	  if ( verbosity_ ) {
            std::cout << " Parameter #" << iParameter << " (" << minuit_->GetParName(iParameter) << "): ";
            if ( minuit_->IsFixed(iParameter) ) std::cout << "LOCKED";
            else std::cout << "FITTED";
            std::cout << std::endl;
          }
        }

        minuitNumFreeParameters_ = minuit_->GetNumberFreeParameters();
        minuitNumFixedParameters_ = minuit_->GetNumberTotalParameters() - minuitNumFreeParameters_;

        if ( verbosity_ ) {
          std::cout << " minuitNumParameters = " << minuit_->GetNumberTotalParameters()
		    << " (free = " << minuitNumFreeParameters_ << ", fixed = " << minuitNumFixedParameters_ << ")" << std::endl;
        }
        assert((minuitNumFreeParameters_ + minuitNumFixedParameters_) == minuitNumParameters_);

        if ( verbosity_ ) std::cout << " BEGINNING FIT ITERATION #" << fitIteration_ << std::endl;

        minuitStatus_ = minuit_->Minimize();

        // Get stats about the minimizer;
        Double_t dummyErrDef;
        Int_t dummyNpvar, dummyNparx;
        minuit_->GetStats(minuitNll_, minuitEdm_, dummyErrDef, dummyNpvar, dummyNparx);

        if ( verbosity_ ) {
          std::cout << " DONE WITH FIT ITERATION #" << fitIteration_
            << " FIT RESULT: " << minuitStatus_
            << " NLL: " << minuitNll_ << " EDM: " << minuitEdm_ << std::endl;
          std::cout << " COVARIANCE MATRIX:" << std::endl;
        }
        if (dqmStore_) {
          MEMap::const_iterator statMonitorIter = monitors_.find("mstatus");
          if (statMonitorIter != monitors_.end()) {
            statMonitorIter->second->getTH1F()->SetBinContent(
                indexFitFunctionCall_ + 1, minuitStatus_);
          }
          MEMap::const_iterator edmMonitorIter = monitors_.find("edm");
          if (edmMonitorIter != monitors_.end()) {
            edmMonitorIter->second->getTH1F()->SetBinContent(
                indexFitFunctionCall_ + 1, minuitEdm_);
          }
        }
      }

      edm::LogInfo("SVfitAlgorithm::fit") << " Minuit fit Status = " << minuitStatus_ << std::endl;

      for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
        minuitFittedParameterValues_[iParameter] = minuit_->GetParameter(iParameter);
        if ( verbosity_ ) {
          std::cout << " Parameter #" << iParameter <<  ":"
		    << " " << minuit_->GetParName(iParameter) << " = "
		    << minuitFittedParameterValues_[iParameter] << std::endl;
        }
      }

      char dummyToGetName[100];
      // Do scans of each parameter that was fitted.
      if (dqmStore_ && monitorCfg_.getParameter<bool>("doScan")) {
        // We dont' want to monitor calls to the likelihood while we are messing
        // around with the scans.
        enableLikelihoodMonitoring_ = false;
        for (unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter) {
          // Check if we fitted this parameter
          std::string parName = minuit_->GetParName(iParameter);
          if (fittedParameterNames.count(parName)) {
            // book our scan of this parameter
            Double_t parValue, parErr, vLow, vHigh;
            minuit_->GetParameter(iParameter, dummyToGetName,
                parValue, parErr, vLow, vHigh);
            double xlow, xhigh;
            if (vLow == 0. && vHigh == 0.) {
              // Unbounded, define using errors
              xlow = parValue - 5*parErr;
              xhigh = parValue + 5*parErr;
            } else {
              xlow = vLow;
              xhigh = vHigh;
            }
            std::string meName = meNameBase.str() + parName + "_scan" ;
            unsigned int steps = 500;
            MonitorElement* monitor = dqmStore_->book1D(
              meName, meName, steps, xlow, xhigh);
            double stepsize = (xhigh - xlow)/steps;
            // Make a copy of the paramters
            vector<double> myxs = minuitFittedParameterValues_;
            // Scan over the parameter
            for(size_t ix = 0; ix < steps; ++ix) {
              double newx = (ix + 0.5)*stepsize;
              myxs[iParameter] = newx;
              double nll = negLogLikelihood(myxs);
              monitor->Fill(newx, nll);
            }
          }
        }
        // renable likelihood monitoring for next fit.
        enableLikelihoodMonitoring_ = true;
      } // end monitor scan

      applyParameters(currentDiTauSolution_, minuitFittedParameterValues_);

      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        double likelihoodFunctionValue = (**likelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
        currentDiTauSolution_.negLogLikelihoods_.insert(std::make_pair((*likelihoodFunction)->name(), likelihoodFunctionValue));
      }

      currentDiTauSolution_.minuitStatus_ = minuitStatus_;

      if ( numSamplings_ > 0 ) compErrorEstimates();
    }

    void applyParameters(SVfitDiTauSolution& diTauSolution, const std::vector<double>& x) const
    {
//--- set primary event vertex position (tau lepton production vertex)
      diTauSolution.eventVertexPositionShift_(0) = x[SVfit_namespace::kPrimaryVertexShiftX];
      diTauSolution.eventVertexPositionShift_(1) = x[SVfit_namespace::kPrimaryVertexShiftY];
      diTauSolution.eventVertexPositionShift_(2) = x[SVfit_namespace::kPrimaryVertexShiftZ];

//--- build first tau decay "leg"
      applyParametersToLeg(SVfit_namespace::kLeg1thetaRest, diTauSolution.leg1_,
          x, diTauSolution.eventVertexPosSVrefitted(),
          parameterizeVertexAlongTrackLeg1_);

//--- build second tau decay "leg"
      applyParametersToLeg(SVfit_namespace::kLeg2thetaRest, diTauSolution.leg2_,
          x, diTauSolution.eventVertexPosSVrefitted(),
          parameterizeVertexAlongTrackLeg2_);
    }

    void applyParametersToLeg(int index0, SVfitLegSolution& legSolution,
        const std::vector<double>& x,
        const AlgebraicVector3& eventVertexPos,
        bool parameterizeVertexAlongTrack) const
    {
      int legOffset = index0 - SVfit_namespace::kLeg1thetaRest;

      double gjAngle  = x[legOffset + SVfit_namespace::kLeg1thetaRest];
      double phiLab   = x[legOffset + SVfit_namespace::kLeg1phiLab];
      double massNuNu = x[legOffset + SVfit_namespace::kLeg1nuInvMass];

      const reco::Candidate::LorentzVector& p4Vis = legSolution.p4Vis();

      // Compute the tau momentum in the rest frame
      double pVisRestFrame = SVfit_namespace::pVisRestFrame(p4Vis.mass(), massNuNu);
      // Get the opening angle in the lab frame
      double angleVisLabFrame = SVfit_namespace::gjAngleToLabFrame(pVisRestFrame, gjAngle, p4Vis.P());
      // Compute the tau momentum in the lab frame
      double momentumLabFrame = SVfit_namespace::tauMomentumLabFrame(p4Vis.mass(), pVisRestFrame, gjAngle, p4Vis.P());

      reco::Candidate::Vector direction;

//--- if parameterizeVertexAlongTrack is enabled, compute approximate position
//    of tau lepton decay vertex as intersection of tau momentum vector with track
//    and take phiLab and decayDistanceToTrack parameters as corrections relative to that position,
//    else compute tau lepton decay vertex from thetaRest, phiLab and decayDistanceLab parameters directly
      if ( !parameterizeVertexAlongTrack ) {
	double flightDistance = SVfit_namespace::square(x[legOffset + SVfit_namespace::kLeg1decayDistanceLab]);
	//std::cout << " flightDistance = " << flightDistance << std::endl;

        // Determine the direction of the tau
        direction = SVfit_namespace::tauDirection(p4Vis.Vect().Unit(), angleVisLabFrame, phiLab);
        //std::cout << " direction: x = " << direction.x() << ", y = " << direction.y() << ", z = " << direction.z() << std::endl;

        // Set the flight path
        legSolution.decayVertexPos_(0) = eventVertexPos(0) + flightDistance*direction.x();
        legSolution.decayVertexPos_(1) = eventVertexPos(1) + flightDistance*direction.y();
        legSolution.decayVertexPos_(2) = eventVertexPos(2) + flightDistance*direction.z();
      } else {
        // Find the approximate decay vertex given the track information.
	double flightDistanceCorr = x[legOffset + SVfit_namespace::kLeg1decayDistanceLab];

        SVfit::track::VertexOnTrackFinder svFinder(legSolution);
        // Get the reconstructed vertex and the direction
        GlobalPoint secondaryVertexOnTrack = svFinder.decayVertex(
            eventVertexPos, angleVisLabFrame,
            phiLab, flightDistanceCorr, &direction);

        legSolution.decayVertexPos_(0) = secondaryVertexOnTrack.x();
        legSolution.decayVertexPos_(1) = secondaryVertexOnTrack.y();
        legSolution.decayVertexPos_(2) = secondaryVertexOnTrack.z();

        if ( verbosity_ ) std::cout << "Found vertex via track @ " << secondaryVertexOnTrack << std::endl;
      }

      // Throw a warning if the tau direction is opposite to that of the visible
      // momentum.
      if (direction.Unit().Dot(p4Vis.Vect().Unit()) < 0) {
        edm::LogWarning("BadTauDirection")
          << "The computed tau direction points opposite to the vis momentum!"
          << " tau: " << direction.Unit() << " vis: " << p4Vis.Vect().Unit();
      }

      // Compute the tau P4
      reco::Candidate::LorentzVector tauP4 = SVfit_namespace::tauP4(direction.Unit(), momentumLabFrame);
      if ( verbosity_ ) std::cout << "tauP4: p = " << tauP4.P() << ","
				  << " eta = " << tauP4.eta() << ", phi = " << tauP4.phi() << std::endl;

      // Build the tau four vector. By construction, the neutrino is tauP4 - visP4
      legSolution.p4Invis_ = tauP4 - p4Vis;

      // Build boost vector and compute the rest frame quanitites
      legSolution.p4VisRestFrame_ = SVfit_namespace::boostToCOM(tauP4, legSolution.p4Vis_);
      legSolution.p4InvisRestFrame_ = SVfit_namespace::boostToCOM(tauP4, legSolution.p4Invis_);

      // Set meson decay angles for tau- --> rho- nu --> pi- pi0 nu
      // and tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decay modes
      // (needed in case likelihood functions for decays of polarized tau leptons are used)
      legSolution.thetaVMrho_ = x[legOffset + SVfit_namespace::kLeg1thetaVMrho];
      legSolution.thetaVMa1_  = x[legOffset + SVfit_namespace::kLeg1thetaVMa1];
      legSolution.thetaVMa1r_ = x[legOffset + SVfit_namespace::kLeg1thetaVMa1r];
      legSolution.phiVMa1r_   = x[legOffset + SVfit_namespace::kLeg1phiVMa1r];
    }

    void compErrorEstimates() const
    {
      //std::cout << "<SVfitAlgorithm::compErrorEstimates>:" << std::endl;

//--- compute estimate for error matrix using Minuit's HESSE algorithm
      minuit_->ExecuteCommand("HESSE", 0, 0);

//--- copy error matrix estimated by Minuit into TMatrixDSym object;
//    only copy elements corresponding to "free" fit parameters,
//    to avoid problems with Cholesky decomposition
      std::vector<unsigned> lutFreeToMinuitParameters(minuitNumFreeParameters_);
      std::vector<int> lutMinuitToFreeParameters(minuitNumParameters_);
      unsigned freeParameter_index = 0;
      for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
        if ( !minuit_->IsFixed(iParameter) ) {
          lutMinuitToFreeParameters[iParameter] = freeParameter_index;
          lutFreeToMinuitParameters[freeParameter_index] = iParameter;
          ++freeParameter_index;
        } else {
          lutMinuitToFreeParameters[iParameter] = -1;
        }
      }

      TMatrixDSym freeErrorMatrix(minuitNumFreeParameters_);
      for ( unsigned iRow = 0; iRow < minuitNumFreeParameters_; ++iRow ) {
        //unsigned minuitRow_index = lutFreeToMinuitParameters[iRow];
        //std::cout << "iRow = " << iRow << ": minuitRow_index = " << minuitRow_index << std::endl;
        for ( unsigned iColumn = 0; iColumn < minuitNumFreeParameters_; ++iColumn ) {
          //unsigned minuitColumn_index = lutFreeToMinuitParameters[iColumn];
          //std::cout << "iColumn = " << iColumn << ": minuitColumn_index = " << minuitColumn_index << std::endl;
          //--- CV: note that Minuit error matrix
          //        contains "free" fit parameters only
          //       (i.e. is of dimension nF x nF, where nF is the number of "free" fit parameters)
          //freeErrorMatrix(iRow, iColumn) = minuit_->GetCovarianceMatrixElement(minuitRow_index, minuitColumn_index);
          freeErrorMatrix(iRow, iColumn) = minuit_->GetCovarianceMatrixElement(iRow, iColumn);
        }
      }
      //freeErrorMatrix.Print();

//--- decompose "physical" error matrix A
//    into "square-root" matrix A, with U * U^T = A
      TDecompChol choleskyDecompAlgorithm;
      choleskyDecompAlgorithm.SetMatrix(freeErrorMatrix);
      choleskyDecompAlgorithm.Decompose();
      TMatrixD freeErrorMatrix_sqrt = choleskyDecompAlgorithm.GetU();
      //freeErrorMatrix_sqrt.Print();

//--- generate random variables distributed according to N-dimensional normal distribution
//    for mean vector m (vector of best fit paramaters determined by Minuit)
//    and covariance V (error matrix estimated by Minuit)
//
//    NB: correlations in the N-dimensional normal distribution
//        are generated by the affine transformation
//          rndCorrelated = mu + U * rndUncorrelated
//        described in section "Drawing values from the distribution"
//        of http://en.wikipedia.org/wiki/Multivariate_normal_distribution
//
      std::vector<double> massValues(numSamplings_);
      std::vector<double> x1Values(numSamplings_);
      std::vector<double> x2Values(numSamplings_);
      TVectorD rndFreeParameterValues(minuitNumFreeParameters_);
      std::vector<double> rndParameterValues(minuitNumParameters_);
      SVfitDiTauSolution rndDiTauSolution(currentDiTauSolution_);
      int iSampling = 0;
      while ( iSampling < numSamplings_ ) {
        for ( unsigned iFreeParameter = 0; iFreeParameter < minuitNumFreeParameters_; ++iFreeParameter ) {
          rndFreeParameterValues(iFreeParameter) = rnd_.Gaus(0., 1.);
        }

        rndFreeParameterValues *= freeErrorMatrix_sqrt;

        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          int freeParameter_index = lutMinuitToFreeParameters[iParameter];
          if ( freeParameter_index != -1 ) {
            rndParameterValues[iParameter] = minuitFittedParameterValues_[iParameter] + rndFreeParameterValues(freeParameter_index);
          } else {
            rndParameterValues[iParameter] = minuitFittedParameterValues_[iParameter];
          }
        }

        applyParameters(rndDiTauSolution, rndParameterValues);

        double mass = rndDiTauSolution.mass();
        double x1 = rndDiTauSolution.leg1().x();
        double x2 = rndDiTauSolution.leg2().x();

        if ( !(TMath::IsNaN(mass) || TMath::IsNaN(x1) || TMath::IsNaN(x2)) ) {
          massValues[iSampling] = mass;
          x1Values[iSampling] = x1;
          x2Values[iSampling] = x2;
          ++iSampling;
        }
      }

      std::sort(massValues.begin(), massValues.end());
      std::sort(x1Values.begin(), x1Values.end());
      std::sort(x2Values.begin(), x2Values.end());

      int median_index = TMath::Nint(0.50*numSamplings_);
      int oneSigmaUp_index = TMath::Nint(0.84*numSamplings_);
      int oneSigmaDown_index = TMath::Nint(0.16*numSamplings_);

      currentDiTauSolution_.hasErrorEstimates_ = true;
      currentDiTauSolution_.massErrUp_ = massValues[oneSigmaUp_index] - massValues[median_index];
      currentDiTauSolution_.massErrDown_ = massValues[median_index] - massValues[oneSigmaDown_index];

      currentDiTauSolution_.leg1_.hasErrorEstimates_ = true;
      currentDiTauSolution_.leg1_.xErrUp_ = x1Values[oneSigmaUp_index] - x1Values[median_index];
      currentDiTauSolution_.leg1_.xErrDown_ = x1Values[median_index] - x1Values[oneSigmaDown_index];

      currentDiTauSolution_.leg2_.hasErrorEstimates_ = true;
      currentDiTauSolution_.leg2_.xErrUp_ = x2Values[oneSigmaUp_index] - x2Values[median_index];
      currentDiTauSolution_.leg2_.xErrDown_ = x2Values[median_index] - x2Values[oneSigmaDown_index];
    }

    std::string name_;
    bool parameterizeVertexAlongTrackLeg1_;
    bool parameterizeVertexAlongTrackLeg2_;

    SVfitEventVertexRefitter* eventVertexRefitAlgorithm_;
    SVfitLegTrackExtractor<T1> leg1TrackExtractor_;
    SVfitLegTrackExtractor<T2> leg2TrackExtractor_;
    SVfitLegHasNeutralsExtractor<T1> leg1NeutralActivity_;
    SVfitLegHasNeutralsExtractor<T2> leg2NeutralActivity_;

    std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> likelihoodFunctions_;
    bool likelihoodsSupportPolarization_;

    mutable const CompositePtrCandidateT1T2MEt<T1,T2>* currentDiTau_;
    mutable SVfitDiTauSolution currentDiTauSolution_;

    mutable TFitterMinuit* minuit_;
    SVfitMinuitFCNadapter<T1,T2> minuitFCNadapter_;
    unsigned minuitMaxInterations_;
    const static unsigned minuitNumParameters_ = 19;
    mutable unsigned minuitNumFreeParameters_;
    mutable unsigned minuitNumFixedParameters_;
    mutable std::vector<double> minuitFittedParameterValues_;

    int numSamplings_;
    mutable TRandom3 rnd_;

    mutable long indexFitFunctionCall_;

    // Class level variables for extracting monitoring from minuit.
    int minuitStatus_;
    Double_t minuitEdm_, minuitNll_;

    int verbosity_;

    // Instrumentation
    // flag to enable monitoring in negLogLikelihood()
    bool enableLikelihoodMonitoring_;
    DQMStore* dqmStore_;
    // Keep track of the current event and how many tau pairs we have processed.
    mutable edm::EventID currentEvent_;
    mutable unsigned int tauPairIndex_;
    edm::ParameterSet monitorCfg_;

    typedef std::map<std::string, MonitorElement*> MEMap;
    // Monitors for the current fit
    MEMap monitors_;
};

#endif
