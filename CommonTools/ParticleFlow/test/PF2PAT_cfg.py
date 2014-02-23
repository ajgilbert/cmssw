# official example for PF2PAT

import FWCore.ParameterSet.Config as cms

process = cms.Process("PF2PAT")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)


# the following line does not work anymore:
# process.load("CommonTools.ParticleFlow.Sources/Data/source_124120_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      # you need to be a CAF user to be able to use these files. They're good to test pile-up
      '/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-PromptSkim-v6/0000/FE65D2DC-67C5-E011-AF6C-002618943894.root'
      # but you can always use a relval like this one:
      # '/store/relval/CMSSW_4_2_3/RelValZTT/GEN-SIM-RECO/START42_V12-v2/0062/4CEA9C47-287B-E011-BAB7-00261894396B.root')
      ))
    
     
                   
print process.source

# path ---------------------------------------------------------------


process.load("CommonTools.ParticleFlow.PF2PAT_cff")

from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT

# the following is advocated by JetMET, but leads to include very far tracks in the no pile up collection
enablePileUpCorrectionInPF2PAT( process, postfix='')

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.PF2PAT 
    )

# output ------------------------------------------------------------

#process.load("FastSimulation.Configuration.EventContent_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.pf2pat = cms.OutputModule("PoolOutputModule",
                                  #    process.AODSIMEventContent,
                                  outputCommands = cms.untracked.vstring('drop *'),
                                  fileName = cms.untracked.string('PF2PAT.root')
)
process.aod = cms.OutputModule("PoolOutputModule",
                               process.AODSIMEventContent,
                               #outputCommands = cms.untracked.vstring('drop *'),
                               fileName = cms.untracked.string('aod.root')
)
process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
process.pf2pat.outputCommands.extend( process.PF2PATEventContent.outputCommands )
process.pf2pat.outputCommands.extend( process.PF2PATStudiesEventContent.outputCommands )


process.outpath = cms.EndPath(
    process.pf2pat 
#    process.aod
    )


# other stuff

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# the following are necessary for taus:

# process.load("Configuration.StandardSequences.GeometryPilot2_cff")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# do not forget to set the global tag according to the
# release you are using and the sample you are reading (data or MC)
# global tags can be found here:
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Valid_Global_Tags_by_Release
# process.GlobalTag.globaltag = cms.string('GR09_R_34X_V2::All')
