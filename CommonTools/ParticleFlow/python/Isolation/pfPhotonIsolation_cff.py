import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.photonPFIsolationDeposits_cff import *
from CommonTools.ParticleFlow.Isolation.photonPFIsolationValues_cff import *

pfPhotonIsolationSequence = cms.Sequence(
    photonPFIsolationDepositsSequence +
    photonPFIsolationValuesSequence
    )

