import FWCore.ParameterSet.Config as cms


from CommonTools.ParticleFlow.Isolation.muonPFIsolationDeposits_cff import *
from CommonTools.ParticleFlow.Isolation.muonPFIsolationValues_cff import *

muonPFIsolationSequence =  cms.Sequence(
    muonPFIsolationDepositsSequence + 
    muonPFIsolationValuesSequence
)                                         





                 

