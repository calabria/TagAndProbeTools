import FWCore.ParameterSet.Config as cms

probeAnyTau = cms.EDFilter("PATTauRefSelector",
                         src = cms.InputTag("selectedPatTaus"),
                         cut = cms.string('pt > 20.0 && abs(eta) < 2.3 && tauID("decayModeFinding")'),
                         filter = cms.bool(False)
                         )

passingIsoLooseMuonVetoLoose = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLoose") > 0.5')
 	)

passingIsoLooseMuonVetoMedium = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonMedium") > 0.5')
	)

passingIsoLooseMuonVetoTight = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTight") > 0.5')
	)

passingIsoLooseMuonVetoLoose2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLoose2") > 0.5')
 	)

passingIsoLooseMuonVetoMedium2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonMedium2") > 0.5')
	)

passingIsoLooseMuonVetoTight2 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTight2") > 0.5')
	)

passingIsoLooseMuonVetoLoose3 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLoose3") > 0.5')
 	)

passingIsoLooseMuonVetoTight3 = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTight3") > 0.5')
	)

passingIsoLooseMuonVetoLooseMVA = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonLooseMVA") > 0.5')
 	)

passingIsoLooseMuonVetoMediumMVA = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonMediumMVA") > 0.5')
	)

passingIsoLooseMuonVetoTightMVA = probeAnyTau.clone(
	cut = cms.string(probeAnyTau.cut.value() + '&& tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstMuonTightMVA") > 0.5')
	)

passingProbes = cms.Sequence(probeAnyTau* 
			    (passingIsoLooseMuonVetoLoose +
			     passingIsoLooseMuonVetoMedium +
			     passingIsoLooseMuonVetoTight 
			    )*
			    (passingIsoLooseMuonVetoLoose2 +
			     passingIsoLooseMuonVetoMedium2 +
			     passingIsoLooseMuonVetoTight2 
			    )*
			    (passingIsoLooseMuonVetoLoose3 +
			     passingIsoLooseMuonVetoTight3 
			    )*
			    (passingIsoLooseMuonVetoLooseMVA +
			     passingIsoLooseMuonVetoMediumMVA +
			     passingIsoLooseMuonVetoTightMVA 
			    )
			   )
