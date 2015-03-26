
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
                              qualityString = cms.string("highPurity"),
                              cutDzErrMax = cms.untracked.double(3.0),
                              cutDxyErrMax = cms.untracked.double(3.0),
                              cutPtErrMax = cms.untracked.double(0.1),
                              vertexZMax = cms.double(15.0),
                              cutMultMin = cms.double(0.0),
                              cutMultMax = cms.double(1000.0),
                              cutMinTrack = cms.double(0.4),
                              ptBins = cms.vdouble(1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0,
                                                   5.5, 6.0, 7.0, 8.0, 10.0, 12.0),
                              NptBins = cms.vdouble(1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0),
                              trackSrc = cms.InputTag("generalTracks"),
                              vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
                              jetSrc = cms.InputTag('akPu3PFpatJets'),
                              etaMin = cms.double(-3.0),
                              etaMax = cms.double(3.0),
                              ptMin = cms.double(0.0),
                              etaBins = cms.vdouble(-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
                                                    0.4, 0.8, 1.2, 1.6, 2.0, 2.4),
                              occByLeadingJetEt = cms.bool(False),
                              applyCuts = cms.bool(True),
                              dxyErrMax = cms.double(3.0),
                              dzErrMax = cms.double(3.0),
                              ptErrMax = cms.double(0.1),
                              applyJetCuts = cms.bool(False),
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("pPb_ntuple.root")
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
process.load('AwaySideRun2015Ana.AwayAnalyzer.PAPileUpVertexFilter_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.path = cms.Path(process.hltSingleTrigger *
                        process.PAcollisionEventSelection *
                        process.siPixelRecHits *
                        # process.pileupVertexFilterCutGplus *
                        process.pACentrality *
                        process.AwaySideAnalysis
                        )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                        '/store/hidata/HIRun2013/PAHighPt/RECO/FlowCorrPA-PromptSkim-v2/00005/F0A4F448-A477-E211-A550-003048F316C8.root'
                                                              )
                            )


