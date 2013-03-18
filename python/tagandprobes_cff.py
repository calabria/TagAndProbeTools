import FWCore.ParameterSet.Config as cms

tnpAnyMuAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                          decay = cms.string("tagAnyMu probeAnyTau"),
			  roles = cms.vstring('muon', 'tau'),
                          cut   = cms.string("30 < mass < 120"),
                          checkCharge = cms.bool(False),   
                          )
