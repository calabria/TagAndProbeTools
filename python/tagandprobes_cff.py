import FWCore.ParameterSet.Config as cms

tnpAnyMuAnyTau = cms.EDProducer("CandViewShallowCloneCombiner",
                          decay = cms.string("tagAnyMu probeAnyTau"),
			  roles = cms.vstring('muon', 'tau'),
                          cut   = cms.string("(abs(daughter('muon').vz - daughter('tau').vz) < 0.14) && (deltaR(daughter('muon').eta,daughter('muon').phi,daughter('tau').eta,daughter('tau').phi) > 0.5)"),
                          checkCharge = cms.bool(False),   
                          )
