import FWCore.ParameterSet.Config as cms

longlivedanalyzer = cms.EDAnalyzer('LongLivedAnalysis',
		HepMCProduct = cms.InputTag("generatorSmeared")    
)




