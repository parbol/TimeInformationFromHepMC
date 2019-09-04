from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = '1500-494-160_DisplacedSUSY_NTuple'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runLongLivedAnalysis_cfg.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['output.root']
config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP2_13TeV_80X_18052019-0257/fernance-DisplacedSUSY_RunIISummer16MiniAODv3_24052019-0256-bd3e7bcff6c9bcad356ea4ed7e4f08b4/USER'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fernance/' 

config.section_('Site')
config.Site.storageSite = 'T2_ES_IFCA'
