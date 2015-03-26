from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = False
config.section_('JobType')
config.JobType.psetName = 'ecalflowntpNew_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['pPb_ntuple.root']
config.section_('Data')
config.Data.inputDataset = '/PAMinBiasUPC/HIRun2013-PromptReco-v1/RECO'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_210676_HI_PromptReco_Collisions13_JSON.txt'
config.section_('User')
config.section_('Site')
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_US_Vanderbilt'
