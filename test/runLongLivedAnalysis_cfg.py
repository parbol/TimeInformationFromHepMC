import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalyzer.IFCALongLivedAnalysis.IFCALongLivedAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v12'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root"
       #'file:/afs/cern.ch/work/p/pablom/public/1686A035-14E9-E811-BCC8-0242AC130002.root'
       # 'file:test/merged.root'
       [
	'file:/afs/cern.ch/work/p/pablom/private/TimingFixForCandidate/CMSSW_10_4_0_mtd5/src/MyAnalyzer/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1000_MChi_148_ctau_60mm_TuneCP2_13TeV_pythia8_80X_994.root' 
       ]
    )
)

process.p = cms.Path(process.longlivedanalyzer)


