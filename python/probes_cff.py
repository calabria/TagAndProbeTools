import FWCore.ParameterSet.Config as cms

probeAnyTau = cms.EDFilter("PATTauRefSelector",
                         src = cms.InputTag("selectedPatTaus"),
                         cut = cms.string('pt > 20.0 && abs(eta) < 2.3'),
                         filter = cms.bool(False)
                         )

passingIsoLooseMuonVetoLoose = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLoose") > 0.5')
 	)

passingIsoLooseMuonVetoMedium = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonMedium") > 0.5')
	)

passingIsoLooseMuonVetoTight = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTight") > 0.5')
	)

passingIsoLooseMuonVetoLoose2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLoose2") > 0.5')
 	)

passingIsoLooseMuonVetoMedium2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonMedium2") > 0.5')
	)

passingIsoLooseMuonVetoTight2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTight2") > 0.5')
	)

passingProbes = cms.Sequence(probeAnyTau* 
			    (passingIsoLooseMuonVetoLoose +
			     passingIsoLooseMuonVetoMedium +
			     passingIsoLooseMuonVetoTight 
			    )*
			    (passingIsoLooseMuonVetoLoose2 +
			     passingIsoLooseMuonVetoMedium2 +
			     passingIsoLooseMuonVetoTight2 
			    )
			   )
