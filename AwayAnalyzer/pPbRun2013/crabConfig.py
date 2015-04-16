from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'AwayAnalysisEtaAsso_0_1_2_Pt3_10'
config.General.workArea = 'pPbRun2013'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'awaySideAnalyzer_PAHighPt_cfg.py'
config.JobType.outputFiles = ['pPb_Pbp_CombinedDiHadronCorrFunc_NoCutetaTrg_EtaAsso_0_1_2_Pt3_10.root']

config.section_('Data')
config.Data.inputDataset = '/PAHighPt/HIRun2013-PromptReco-v1/RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'Cert_pPb_Pbpcombined_HI_PromptReco_Collisions13_JSON.txt'
config.Data.publication = False
#config.Data.runRange = '193093-193999' # '193093-194075'
#config.Data.outLFN = '/store/user/<subdir>' # or '/store/group/<subdir>'
#config.Data.publishDataName = 'CRAB3_tutorial_Data_analysis_test5'

config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'