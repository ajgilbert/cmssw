import FWCore.ParameterSet.Config as cms

# iso deposits and isolation values, defined by the Muon POG

# this file is temporarily moved to this package.
# in 4_4_X, we load this file from the RecoMuon/MuonIsolation package
from  CommonTools.ParticleFlow.Isolation.muonPFIsolation_cff import *

# computing the isolation for the muons produced by PF2PAT, and not for reco muons
sourceMuons = 'pfSelectedMuons'



muPFIsoDepositCharged.src = sourceMuons
muPFIsoDepositChargedAll.src = sourceMuons
muPFIsoDepositNeutral.src = sourceMuons
muPFIsoDepositGamma.src = sourceMuons
muPFIsoDepositPU.src = sourceMuons

# muPFIsoDepositCharged.ExtractorPSet.inputCandView = 'pfAllChargedHadronsIso'
# muPFIsoDepositChargedAll.ExtractorPSet.inputCandView = 'pfAllChargedParticlesIso'
# muPFIsoDepositNeutral.ExtractorPSet.inputCandView = 'pfAllNeutralHadronsIso'
# muPFIsoDepositGamma.ExtractorPSet.inputCandView = 'pfAllPhotonsIso'
# muPFIsoDepositPU.ExtractorPSet.inputCandView = 'pfPileUpAllChargedParticlesIso'

pfMuonIsolationSequence = cms.Sequence(
    muonPFIsolationSequence 
    )
