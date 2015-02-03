import FWCore.ParameterSet.Config as cms

vecMC = (2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06)

vecData = (12261.2, 32854.9, 89717.9, 262466, 546468, 3.03611e+06, 1.74973e+07, 5.15968e+07, 1.21961e+08, 2.4679e+08, 4.34103e+08, 6.74775e+08, 8.72426e+08, 9.86769e+08, 1.05924e+09, 1.11387e+09, 1.14671e+09, 1.15622e+09, 1.14973e+09, 1.13116e+09, 1.10319e+09, 1.06971e+09, 1.03229e+09, 9.86238e+08, 9.23288e+08, 8.38926e+08, 7.36123e+08, 6.22188e+08, 5.04938e+08, 3.92125e+08, 2.90994e+08, 2.06419e+08, 1.3992e+08, 9.04343e+07, 5.55578e+07, 3.23672e+07, 1.78798e+07, 9.39093e+06, 4.71572e+06, 2.28233e+06, 1.07582e+06, 500386, 233316, 111004, 54799.4, 28395.5, 15488.2, 8844.91, 5236.19, 3180.09, 1964.04, 1225.15, 767.779, 481.279, 300.644, 186.558, 114.687, 69.6938, 41.7929, 24.6979)

addUserVariables = cms.EDProducer("UserDefinedVariables",
    	objects = cms.InputTag("muonVariables"),
    	triggerResults = cms.InputTag("patTriggerEvent"),
    	met = cms.InputTag("patPFMETsTypeIcorrected"),
    	isMC = cms.bool(False), #####
    	TrueDist2011A = cms.vdouble(vecData),
    	TrueDist2011B = cms.vdouble(vecData),
    	MCDist = cms.vdouble(vecMC),
    	)

tagAnyMu = cms.EDFilter("PATMuonRefSelector",
	src = cms.InputTag("muonVariables"),
	cut = cms.string("pt > 24.0 && abs(eta) < 2.1 && abs(userFloat('dxyWrtPV')) < 0.045 && abs(userFloat('dzWrtPV')) < 0.2" +
                         " && ("+
                         "(isGlobalMuon"+
                         " && globalTrack.isNonnull "+
                         " && globalTrack.normalizedChi2 < 10"+
                         " && globalTrack.hitPattern.numberOfValidMuonHits > 0"+
                         " && numberOfMatchedStations > 1"+
                         " && innerTrack.hitPattern.numberOfValidPixelHits > 0"+
                         " && track.hitPattern.trackerLayersWithMeasurement > 5)"+
                         " || userInt('isPFMuon') > 0.5)"),
	filter = cms.bool(False)
 	)
