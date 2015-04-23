import FWCore.ParameterSet.Config as cms

import copy
import os
import re

process = cms.Process("makeToyTaus")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.load('Configuration.Geometry.GeometryExtended2017Reco_cff') 
#process.load('Configuration.Geometry.GeometryExtended2017_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
##from Configuration.AlCa.autoCond import autoCond
##print "autoCond = %s" % autoCond
##process.GlobalTag.globaltag = cms.string( autoCond[ 'upgrade2017' ] )
process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(                                
        'file:/Volumes/HDD/tuples/ZprimeToTauTau_M-1000_Tune4C_13TeV-pythia8-PU40bx25_PHYS14_25_V1-v2.root'                                
    ),
)

process.options   = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

process.dumpPATTauSequence = cms.Sequence()

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.dumpPATTauSequence += process.tauGenJets
process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
    'oneProng2Pi0', 
    'oneProngOther',
    'threeProng0Pi0', 
    'threeProng1Pi0', 
    'threeProngOther', 
    'rare'
)
process.tauGenJetsSelectorAllHadrons.filter = cms.bool(False)
process.dumpPATTauSequence += process.tauGenJetsSelectorAllHadrons

process.selectedTauGenJets = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    cut = cms.string("pt > 20. & abs(eta) < 2.3"),
    filter = cms.bool(False)
)
process.dumpPATTauSequence += process.selectedTauGenJets

#process.dumpGenTaus = cms.EDAnalyzer("DumpGenTauSummary",
#    src = cms.InputTag('genParticles')
#)
#process.dumpPATTauSequence += process.dumpGenTaus

#--------------------------------------------------------------------------------
# produce collection of "toy" PFCandidates 

process.toyPFCandidates = cms.EDProducer("ToyPFCandidateProducer",
    src = cms.InputTag('genParticles')
)
process.dumpPATTauSequence += process.toyPFCandidates
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run jet reconstruction on "toy" PFCandidates

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets

process.ak4ToyPFJets = ak4PFJets.clone(
   src = cms.InputTag('toyPFCandidates', 'pfCandidates')
)
process.dumpPATTauSequence += process.ak4ToyPFJets
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run tau reconstruction on "toy" PFCandidates

process.genTauMatchedPFJets = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak4ToyPFJets'),
    srcNotToBeFiltered = cms.VInputTag('selectedTauGenJets'),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
#process.dumpPATTauSequence += process.genTauMatchedPFJets

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag 
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('particleFlow'), cms.InputTag('toyPFCandidates', 'pfCandidates'))
#massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak4PFJets'), cms.InputTag('genTauMatchedPFJets'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak4PFJets'), cms.InputTag('ak4ToyPFJets'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('generalTracks'), cms.InputTag('toyPFCandidates', 'tracks'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('offlinePrimaryVertices'), cms.InputTag('toyPFCandidates', 'vertices'))

process.dumpPATTauSequence += process.PFTau
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.cleanPatTaus.preselection = cms.string('')
process.cleanPatTaus.checkOverlaps = cms.PSet()
process.cleanPatTaus.finalCut = cms.string("")

process.dumpPATTauSequence += process.makePatTaus
process.dumpPATTauSequence += process.selectedPatTaus
process.dumpPATTauSequence += process.cleanPatTaus
#--------------------------------------------------------------------------------

process.dumpPFCandidates = cms.EDAnalyzer("DumpPFCandidatesByROI",
    srcROIs = cms.InputTag('genTauMatchedPFJets'),
    srcPFCandidates = cms.InputTag('toyPFCandidates', 'pfCandidates'),
    srcVertex = cms.InputTag('toyPFCandidates', 'vertices'),
    dRcone = cms.double(0.3),
    minPt = cms.double(-1.)
)
##process.dumpPATTauSequence += process.dumpPFCandidates

process.dumpPATTaus = cms.EDAnalyzer("DumpPATTauSummary",
    src = cms.InputTag('patTaus')
)
#process.dumpPATTauSequence += process.dumpPATTaus

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("EventTree.root"),
  closeFileFast = cms.untracked.bool(True)
)
import UserCode.ICHiggsTauTau.default_producers_cfi as producers
import UserCode.ICHiggsTauTau.default_selectors_cfi as selectors

process.selectedVertices = cms.EDFilter("VertexRefSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("ndof >= 4 & abs(z) <= 24 & abs(position.Rho) <= 2")
)

process.selectedPFCandidates = cms.EDFilter("PFCandidateRefSelector",
  src = cms.InputTag('particleFlow'),
  cut = cms.string("pt > 0.0")
)

process.selectedToyPFCandidates = cms.EDFilter("PFCandidateRefSelector",
  src = cms.InputTag('toyPFCandidates', 'pfCandidates'),
  cut = cms.string("pt > 0.0")
)

process.icPFProducer = cms.EDProducer('ICPFProducer',
  branch  = cms.string("pfCandidates"),
  input   = cms.InputTag("selectedPFCandidates"),
  requestTracks       = cms.bool(False),
  requestGsfTracks    = cms.bool(False)
)

process.icToyPFProducer = cms.EDProducer('ICPFProducer',
  branch  = cms.string("genPFCandidates"),
  input   = cms.InputTag("selectedToyPFCandidates"),
  requestTracks       = cms.bool(False),
  requestGsfTracks    = cms.bool(False)
)


process.ak4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet')
)
process.ak4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)
process.ak4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)
process.ak4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)

process.selectedPFJets = cms.EDFilter("PFJetRefSelector",
    src = cms.InputTag("ak4PFJets"),
    cut = cms.string("pt > 15")
    )

process.icPFJetProducer = producers.icPFJetProducer.clone(
    branch                    = cms.string("pfJetsPFlow"),
    input                     = cms.InputTag("selectedPFJets"),
    srcConfig = cms.PSet(
      includeJetFlavour         = cms.bool(False),
      inputJetFlavour           = cms.InputTag(""),
      applyJECs                 = cms.bool(True),
      includeJECs               = cms.bool(False),
      JECs                      = cms.PSet(
        L1FastJet  = cms.string("ak4PFL1Fastjet"),
        L2Relative = cms.string("ak4PFL2Relative"),
        L3Absolute = cms.string("ak4PFL3Absolute")
      ),
      applyCutAfterJECs         = cms.bool(False),
      cutAfterJECs              = cms.string(""),
      inputSVInfo               = cms.InputTag(""),
      requestSVInfo             = cms.bool(False),
      BTagDiscriminators        = cms.PSet(
      )
    ),
    destConfig = cms.PSet(
      includePileupID       = cms.bool(False),
      inputPileupID         = cms.InputTag(""),
      includeTrackBasedVars = cms.bool(False),
      inputTracks           = cms.InputTag(""),
      inputVertices         = cms.InputTag(""),
      requestTracks         = cms.bool(False)
    )
)

process.selectedPFTaus = cms.EDFilter("PFTauRefSelector",
  src = cms.InputTag("hpsPFTauProducer", "", "RECO"),
  cut = cms.string("pt > 10.0 & abs(eta) < 3.0")
)

process.selectedToyPFTaus = cms.EDFilter("PFTauRefSelector",
  src = cms.InputTag("hpsPFTauProducer"),
  cut = cms.string("pt > 10.0 & abs(eta) < 3.0")
)

import UserCode.ICHiggsTauTau.tau_discriminators_cfi as tauIDs

process.icTauProducer = producers.icTauProducer.clone(
  branch                  = cms.string("taus"),
  input                   = cms.InputTag("selectedPFTaus"),
  inputVertices           = cms.InputTag("selectedVertices"),
  includeVertexIP         = cms.bool(True),
  requestTracks           = cms.bool(False),
  tauIDs = cms.PSet(
        decayModeFindingNewDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "RECO"),
        decayModeFindingOldDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs", "", "RECO"),
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding", "", "RECO"),
        byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation", "", "RECO"),
        byCombinedIsolationDeltaBetaCorrRaw = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr", "", "RECO"),
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", "", "RECO"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum", "", "RECO"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum", "", "RECO"),
        puCorrPtSum = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum", "", "RECO")
    )
)

process.icToyTauProducer = producers.icTauProducer.clone(
  branch                  = cms.string("toyTaus"),
  input                   = cms.InputTag("selectedToyPFTaus"),
  inputVertices           = cms.InputTag('toyPFCandidates', 'vertices'),
  includeVertexIP         = cms.bool(True),
  requestTracks           = cms.bool(False),
  tauIDs = cms.PSet(
        decayModeFindingNewDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "makeToyTaus"),
        decayModeFindingOldDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingOldDMs", "", "makeToyTaus"),
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding", "", "makeToyTaus"),
        byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation", "", "makeToyTaus"),
        byCombinedIsolationDeltaBetaCorrRaw = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr", "", "makeToyTaus"),
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", "", "makeToyTaus"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationChargedIsoPtSum", "", "makeToyTaus"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauMVA3IsolationNeutralIsoPtSum", "", "makeToyTaus"),
        puCorrPtSum = cms.InputTag("hpsPFTauMVA3IsolationPUcorrPtSum", "", "makeToyTaus")
    )

)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer("ICGenParticlePruner",
  src = cms.InputTag("genParticles", "", "HLT"),
  select = cms.vstring(
    "drop  *",
    "keep status == 3 || status == 22 || status == 23",  # all status 3
    "keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15",  # all charged leptons
    "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",  # all neutrinos
    "keep++ abs(pdgId) == 15",  # keep full tau decay chain
    "keep (4 <= abs(pdgId) <= 5)", # keep heavy flavour quarks
    "keep (400 <= abs(pdgId) < 600) || (4000 <= abs(pdgId) < 6000)", # keep b and c hadrons
    "keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433", # additional c hadrons for jet fragmentation studies
    "keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543" # additional b hadrons for jet fragmentation studies
  )
)

process.icGenParticleProducer = producers.icGenParticleProducer.clone(
  input   = cms.InputTag("prunedGenParticles"),
  includeMothers = cms.bool(True),
  includeDaughters = cms.bool(True)
)


process.icEventInfoProducer = producers.icEventInfoProducer.clone(
  includeJetRho       = cms.bool(False),
  inputJetRho         = cms.InputTag("kt6PFJets", "rho"),
  includeLeptonRho    = cms.bool(False),
  inputLeptonRho      = cms.InputTag("kt6PFJets", "rho"),
  includeVertexCount  = cms.bool(True),
  inputVertices       = cms.InputTag("selectedVertices"),
  includeCSCFilter    = cms.bool(False),
  inputCSCFilter      = cms.InputTag("BeamHaloSummary"),
)

process.icPileupInfoProducer = producers.icPileupInfoProducer.clone()

process.icEventProducer = producers.icEventProducer.clone()

process.icSequence = cms.Sequence(
  process.dumpPATTauSequence+
  process.selectedVertices+
  process.selectedPFCandidates+
  process.selectedToyPFCandidates+
  process.icPFProducer+
  process.icToyPFProducer+
  process.selectedPFJets+
  process.icPFJetProducer+
  process.selectedPFTaus+
  process.selectedToyPFTaus+
  process.icTauProducer+
  process.icToyTauProducer+
  process.prunedGenParticles+
  process.icGenParticleProducer+
  process.icPileupInfoProducer+
  process.icEventInfoProducer+
  process.icEventProducer
)

process.p = cms.Path(process.icSequence)

# processDumpFile = open('makeToyTaus.dump', 'w')
# print >> processDumpFile, process.dumpPython()



