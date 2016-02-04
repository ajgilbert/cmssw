import FWCore.ParameterSet.Config as cms

MAPAnalyzer =cms.EDAnalyzer('MAPAnalyzer',
                           srcVariables = cms.InputTag("MVAMET"),
                           srcVariableNames = cms.InputTag("MVAMET"),
                           variableNamesToSave = cms.vstring(
                                                            "Jet0_Eta",
                                                            "Jet0_M",
                                                            "Jet0_Phi",
                                                            "Jet0_Pt",
                                                            "Jet1_Eta",
                                                            "Jet1_M",
                                                            "Jet1_Phi",
                                                            "Jet1_Pt",
                                                            "Jet2_Eta",
                                                            "Jet2_Pt",
                                                            "NCleanedJets",
                                                            "NVertex",
                                                            "PhiCorrectedRecoil_Phi",
                                                            "PhiCorrectedRecoil_dPhi",
                                                            "PhiCorrectedRecoil_Pt",
                                                            "PhiCorrectedRecoil_sumEt",
                                                            "PhiCorrectedRecoil_sumEtFraction",
                                                            "recoilpatpfMET_Cov00",
                                                            "recoilpatpfMET_Cov11",
                                                            "recoilpatpfMET_Phi",
                                                            "recoilpatpfMET_dPhi",
                                                            "recoilpatpfMET_Pt",
                                                            "recoilpatpfMET_sumEt",
                                                            "recoilpatpfMET_sumEtFraction",
                                                            "recoilpatpfNoPUMET_Cov00",
                                                            "recoilpatpfNoPUMET_Cov11",
                                                            "recoilpatpfNoPUMET_Phi",
                                                            "recoilpatpfNoPUMET_dPhi",
                                                            "recoilpatpfNoPUMET_Pt",
                                                            "recoilpatpfNoPUMET_sumEt",
                                                            "recoilpatpfNoPUMET_sumEtFraction",
                                                            "recoilpatpfPUCorrectedMET_Cov00",
                                                            "recoilpatpfPUCorrectedMET_Cov11",
                                                            "recoilpatpfPUCorrectedMET_Phi",
                                                            "recoilpatpfPUCorrectedMET_dPhi",
                                                            "recoilpatpfPUCorrectedMET_Pt",
                                                            "recoilpatpfPUCorrectedMET_sumEt",
                                                            "recoilpatpfPUCorrectedMET_sumEtFraction",
                                                            "recoilpatpfPUMET_Cov00",
                                                            "recoilpatpfPUMET_Cov11",
                                                            "recoilpatpfPUMET_Phi",
                                                            "recoilpatpfPUMET_dPhi",
                                                            "recoilpatpfPUMET_Pt",
                                                            "recoilpatpfPUMET_sumEt",
                                                            "recoilpatpfPUMET_sumEtFraction",
                                                            "recoilpatpfTrackMET_Cov00",
                                                            "recoilpatpfTrackMET_Cov11",
                                                            "recoilpatpfTrackMET_Phi",
                                                            "recoilpatpfTrackMET_dPhi",
                                                            "recoilpatpfTrackMET_Pt",
                                                            "recoilpatpfTrackMET_sumEt",
                                                            "recoilpatpfTrackMET_sumEtFraction",
                                                            "recoilslimmedMETsPuppi_Cov00",
                                                            "recoilslimmedMETsPuppi_Cov11",
                                                            "recoilslimmedMETsPuppi_Phi",
                                                            "recoilslimmedMETsPuppi_dPhi",
                                                            "recoilslimmedMETsPuppi_Pt",
                                                            "recoilslimmedMETsPuppi_sumEt",
                                                            "recoilslimmedMETsPuppi_sumEtFraction",
                                                            "recoilslimmedMETs_Cov00",
                                                            "recoilslimmedMETs_Cov11",
                                                            "recoilslimmedMETs_Phi",
                                                            "recoilslimmedMETs_dPhi",
                                                            "recoilslimmedMETs_Pt",
                                                            "recoilslimmedMETs_sumEt",
                                                            "recoilslimmedMETs_sumEtFraction",
                                                            "recoilpatpfChargedPUMET_Pt",
                                                            "recoilpatpfChargedPUMET_Phi",
                                                            "recoilpatpfChargedPUMET_dPhi",
                                                            "recoilpatpfChargedPUMET_sumEt",
                                                            "recoilpatpfChargedPUMET_sumEtFraction",
                                                            "recoilpatpfNeutralPUMET_Pt",
                                                            "recoilpatpfNeutralPUMET_Phi",
                                                            "recoilpatpfNeutralPUMET_dPhi",
                                                            "recoilpatpfNeutralPUMET_sumEt",
                                                            "recoilpatpfNeutralPUMET_sumEtFraction",
                                                            "recoilpatpfNeutralPVMET_Pt",
                                                            "recoilpatpfNeutralPVMET_Phi",
                                                            "recoilpatpfNeutralPVMET_dPhi",
                                                            "recoilpatpfNeutralPVMET_sumEt",
                                                            "recoilpatpfNeutralPVMET_sumEtFraction",
                                                            "recoilpatpfNeutralUnclusteredMET_Pt",
                                                            "recoilpatpfNeutralUnclusteredMET_Phi",
                                                            "recoilpatpfNeutralUnclusteredMET_dPhi",
                                                            "recoilpatpfNeutralUnclusteredMET_sumEt",
                                                            "recoilpatpfNeutralUnclusteredMET_sumEtFraction",
                                                            "Boson_Pt", "Boson_Phi", "Boson_M", "Boson_Eta", "Boson_sumET", "Boson_daughter1", "Boson_daughter2",
                                                            "nCombinations"
                                                               ) )
