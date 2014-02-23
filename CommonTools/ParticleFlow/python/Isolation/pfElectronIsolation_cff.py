import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.electronPFIsolationDeposits_cff import *
from CommonTools.ParticleFlow.Isolation.electronPFIsolationValues_cff import *

pfElectronIsolationSequence = cms.Sequence(
    electronPFIsolationDepositsSequence +
    electronPFIsolationValuesSequence
    )


