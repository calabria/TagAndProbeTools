import FWCore.ParameterSet.Config as cms

oneTp = cms.EDFilter("CandViewCountFilter",
                     	src = cms.InputTag("tnpAnyMuAnyTau"),
                     	minNumber = cms.uint32(1)
                     	)
