
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('CondCore.DBCommon.CondDBCommon_cfi');
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.AwaySideAnalysis = cms.EDAnalyzer('AwayAnalyzer',
                              verbose = cms.untracked.bool(True),
                              vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
                              trackSrc = cms.InputTag("generalTracks"),
                              etaMin = cms.double(-3.0),
                              etaMax = cms.double(3.0),
                              ptMin = cms.double(0.4),
                              vertexZMax = cms.double(15.0),
                              qualityString = cms.string("highPurity"),
                              etaBins = cms.vdouble(-2.4, 0.0, 2.4),
                              vzBins = cms.vdouble(-15.0, -13.0, -11.0, -9.0, -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0),
                              NptBins = cms.vdouble(1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0),
                              cutMultMin = cms.double(0.0),
                              cutMultMax = cms.double(1000.0),
                              cutDzErrMax = cms.untracked.double(3.0),
                              cutDxyErrMax = cms.untracked.double(3.0),
                              cutPtErrMax = cms.untracked.double(0.1),
                              etaMinTrg = cms.double(0.0),
                              etaMaxTrg = cms.double(1.2),
                              etaMinAsso = cms.double(0.0),
                              etaMaxAsso = cms.double(1.2),
                              ptMinTrg = cms.double(3.0),
                              ptMaxTrg = cms.double(10.0),
                              ptMinAsso = cms.double(0.4),
                              ptMaxAsso = cms.double(3.0)
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("pPb_Pbp_CombinedEfficiency_DeltaZ2cm.root")
                                   )

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
        centralityVariable = cms.string("HFtowersPlusTrunc"),
        nonDefaultGlauberModel = cms.string(""),
        centralitySrc = cms.InputTag("pACentrality"),
        pPbRunFlip = cms.untracked.uint32(99999999)
        )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('AwaySideAnalysis.AwayAnalyzer.PAPileUpVertexFilter_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.hltMult100 = process.hltHighLevel.clone()
process.hltMult100.HLTPaths = ["HLT_PAPixelTracks_Multiplicity100_v1",
                               "HLT_PAPixelTracks_Multiplicity100_v2"]

process.hltMult130 = process.hltHighLevel.clone()
process.hltMult130.HLTPaths = ["HLT_PAPixelTracks_Multiplicity130_v1",
                               "HLT_PAPixelTracks_Multiplicity130_v2"]

process.hltMult160 = process.hltHighLevel.clone()
process.hltMult160.HLTPaths = ["HLT_PAPixelTracks_Multiplicity160_v1",
                               "HLT_PAPixelTracks_Multiplicity160_v2"]

process.hltMult190 = process.hltHighLevel.clone()
process.hltMult190.HLTPaths = ["HLT_PAPixelTracks_Multiplicity190_v1",
                               "HLT_PAPixelTracks_Multiplicity190_v2"]

process.hltMult100.andOr = cms.bool(True)
process.hltMult100.throw = cms.bool(False)

process.hltMult130.andOr = cms.bool(True)
process.hltMult130.throw = cms.bool(False)

process.hltMult160.andOr = cms.bool(True)
process.hltMult160.throw = cms.bool(False)

process.hltMult190.andOr = cms.bool(True)
process.hltMult190.throw = cms.bool(False)

process.AwaySideAnalysisMult100 = process.AwaySideAnalysis.clone(
                                                       cutMultMin = cms.double(120),
                                                       cutMultMax = cms.double(150)
                                                       )

process.AwaySideAnalysisMult130 = process.AwaySideAnalysis.clone(
                                                       cutMultMin = cms.double(150),
                                                       cutMultMax = cms.double(185)
                                                       )

process.AwaySideAnalysisMult160 = process.AwaySideAnalysis.clone(
                                                       cutMultMin = cms.double(185),
                                                       cutMultMax = cms.double(220)
                                                       )

process.AwaySideAnalysisMult190 = process.AwaySideAnalysis.clone(
                                                       cutMultMin = cms.double(220),
                                                       cutMultMax = cms.double(260)
                                                       )

process.Mult100 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.AwaySideAnalysisMult100
                           )

process.Mult130 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.AwaySideAnalysisMult130
                           )

process.Mult160 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #     process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.AwaySideAnalysisMult160
                           )

process.Mult190 = cms.Path(process.hltSingleTrigger *
                           process.PAcollisionEventSelection *
                           process.siPixelRecHits *
                           #      process.pileupVertexFilterCutGplus *
                           process.pACentrality *
                           process.AwaySideAnalysisMult190
                           )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                        '/store/hidata/HIRun2013/PAHighPt/RECO/FlowCorrPA-PromptSkim-v2/00005/F0A4F448-A477-E211-A550-003048F316C8.root'
                                                              )
                            )


