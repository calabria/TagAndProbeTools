#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooLognormal.h"
#include "RooProdPdf.h"

#include <vector>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"


using namespace std;
using namespace RooFit;

RooDataHist dataSetProducer(TTree * fullTreeSgnCut, RooRealVar mass, float ScaleFactorSgn, float weight){

	TH1F* hMass = new TH1F("hMass","",50,70,120); 
  	if(ScaleFactorSgn != 1) {
  		fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("mass>>hMass");
        cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;
  	RooDataHist sgnDataHist("sgnDataHist", "", RooArgSet(mass), hMass, weight);

	return sgnDataHist;

}

float getCorrectedEntries(TTree * fullTreeSgnCut, RooRealVar mass, float ScaleFactorSgn, int a = 1, int b = 50){

	TH1F* hMass = new TH1F("hMass","",50,70,120);
  	if(ScaleFactorSgn != 1) {
  		fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("mass>>hMass");
        //cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        //cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	return hMass->Integral(a,b);

}

float getError(TTree * fullTreeSgnCut, RooRealVar mass, float ScaleFactorSgn, int a = 1, int b = 50){

	TH1F* hMass = new TH1F("hMass","",50,70,120);
  	if(ScaleFactorSgn != 1) {
  		fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
		//hMass->Scale(ScaleFactorSgn*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("mass>>hMass");
        //cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        //cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	float numRaw = hMass->Integral(a,b);
	float errRaw = sqrt(numRaw);
	float errScaled = errRaw*ScaleFactorSgn*0.987883333*0.985533333*0.9548436313;

	return errScaled;

}

void fitStudyTemplatesFromMC(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	double nBins_                    = 50,
	const double binCenter_          = 0.75,
	const double binWidth_           = 0.75,
	double xLow_                     = 70,
	double xHigh_                    = 120,
	bool doBinned_                   = true,
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	//const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 > 0.5 && tag_muPFIsolation < 0.1",
	//const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 > 0.5 && tag_muPFIsolation < 0.1",
	//const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_IsoMu24_eta2p1 > 0.5 && tag_muPFIsolation < 0.1"
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5"
	){
  
  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  string path = "./InputFileNew/";
  string pathData = "./InputFile/";

  // signal
  TFile fsgn((path + "testTagAndProbe_DYToMuMu_TNP_TNP.root").c_str());
  fsgn.cd("counter");
  TH1F* totalEventsSgn = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsSgn = totalEventsSgn->GetBinContent(1);
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());

  // WJets
  TFile fWJets((path + "testTagAndProbe_WJets_TNP_TNP.root").c_str());
  fWJets.cd("counter");
  TH1F* totalEventsWJets = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWJets = totalEventsWJets->GetBinContent(1);
  TTree *fullTreeWJets = (TTree*)fWJets.Get((tnp_+"/fitter_tree").c_str());

  // Ztautau
  TFile fZtt((path + "testTagAndProbe_DYToTauTau_TNP_TNP.root").c_str());
  fZtt.cd("counter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());

  // TTJets
  TFile fTTJets((path + "testTagAndProbe_TTJets_TNP_TNP.root").c_str());
  fTTJets.cd("counter");
  TH1F* totalEventsTTJets = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsTTJets = totalEventsTTJets->GetBinContent(1);
  TTree *fullTreeTTJets = (TTree*)fTTJets.Get((tnp_+"/fitter_tree").c_str());

  // WW
  TFile fWW((path + "testTagAndProbe_WW_TNP_TNP.root").c_str());
  fWW.cd("counter");
  TH1F* totalEventsWW = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWW = totalEventsWW->GetBinContent(1);
  TTree *fullTreeWW = (TTree*)fWW.Get((tnp_+"/fitter_tree").c_str());

  // WZ
  TFile fWZ((path + "testTagAndProbe_WZ_TNP_TNP.root").c_str());
  fWZ.cd("counter");
  TH1F* totalEventsWZ = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWZ = totalEventsWZ->GetBinContent(1);
  TTree *fullTreeWZ = (TTree*)fWZ.Get((tnp_+"/fitter_tree").c_str());

  // ZZ
  TFile fZZ((path + "testTagAndProbe_ZZ_TNP_TNP.root").c_str());
  fZZ.cd("counter");
  TH1F* totalEventsZZ = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsZZ = totalEventsZZ->GetBinContent(1);
  TTree *fullTreeZZ = (TTree*)fZZ.Get((tnp_+"/fitter_tree").c_str());

  // Data
  TFile fData((pathData + "testTagAndProbe_SingleMu_ABCD.root").c_str());
  TTree *fullTreeData = (TTree*)fData.Get((tnp_+"/fitter_tree").c_str());

  TFile *McP = new TFile("dummy1.root","RECREATE");

  // MC truth Efficiency

  TH1F* hS           = new TH1F("hS","",1,0,150);
  TH1F* hSP          = new TH1F("hSP","",1,0,150);

  fullTreeSgn->Draw("mass>>hS",Form("tag_puMCWeightRun2012*(%s && mass>%f && mass<%f && mcTrue && pair_charge==0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_triggerBitSingleMu > 0.5 )",bin_.c_str(),xLow_,xHigh_));
  double SGNtrue = hS->Integral(0,2);
  fullTreeSgn->Draw("mass>>hSP",Form("tag_puMCWeightRun2012*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && pair_charge==0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_triggerBitSingleMu > 0.5 )",bin_.c_str(),category_.c_str(),cutValue_,xLow_,xHigh_));
  double SGNtruePass = hSP->Integral(0,2);

  double McTruthEff    = SGNtruePass/SGNtrue;
  double BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  delete hS; delete hSP;

  // Create trees with cuts: passing

  std::cout<<"Processing Passing, OS, Low MT"<<std::endl;

  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutTempJet = fullTreeSgn->CopyTree( Form("(mcTrue == 0 && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  std::cout<<"Processing Passing, SS, Low MT"<<std::endl;

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

  std::cout<<"Processing Passing, OS, High MT"<<std::endl;

  //High MT
  TTree* fullTreeDataHiMtCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeSgnHiMtCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsHiMtCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttHiMtCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsHiMtCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWHiMtCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZHiMtCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZHiMtCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );

  //Normalization passing

  float Lumi_ = 19484.55;

  //Cross sections
  float SigmaWJets = 37509.0;
  float SigmaZtt = 1915.083;
  float SigmaTTJets = 234;
  float SigmaSgn = 1915.083;
  float SigmaWW = 56.7532;
  float SigmaWZ = 33.85;
  float SigmaZZ = 8.297;

  float ScaleFactorWJets = Lumi_/(readEventsWJets/SigmaWJets);
  float ScaleFactorZtt = Lumi_/(readEventsZtt/SigmaZtt);
  float ScaleFactorTTJets = Lumi_/(readEventsTTJets/SigmaTTJets);
  float ScaleFactorSgn = Lumi_/(readEventsSgn/SigmaSgn);
  float ScaleFactorWW = Lumi_/(readEventsWW/SigmaWW);
  float ScaleFactorWZ = Lumi_/(readEventsWZ/SigmaWZ);
  float ScaleFactorZZ = Lumi_/(readEventsZZ/SigmaZZ);

  std::cout<<"WJets Ztt TTJets Sgn WW WZ ZZ"<<std::endl;
  std::cout<<readEventsWJets<<" "<<readEventsZtt<<" "<<readEventsTTJets<<" "<<readEventsSgn<<" "<<readEventsWW<<" "<<readEventsWZ<<" "<<readEventsZZ<<std::endl;
  std::cout<<ScaleFactorWJets<<" "<<ScaleFactorZtt<<" "<<ScaleFactorTTJets<<" "<<ScaleFactorSgn<<" "<<ScaleFactorWW<<" "<<ScaleFactorWZ<<" "<<ScaleFactorZZ<<std::endl;

  McP->cd();

  // mass variable

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  //Roofit datasets of data passing and failing cut

  RooDataSet DataDataSetP("DataDataSetP", "dataset for Data pass", RooArgSet(mass), Import( *fullTreeDataCut ) );
  //std::cout << "data dataset Pass " << DataDataSetP.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistP("DataDataHistP", "", RooArgSet(mass), DataDataSetP, 1.0);

  RooDataSet dataDataSetSS("dataDataSetSS","dataset for signal-pass template SS", RooArgSet(mass), Import( *fullTreeDataSSCut ) );
  RooDataHist dataDataHistSS("dataDataHistSS", "", RooArgSet(mass), dataDataSetSS, 1.0);

  RooDataSet dataDataSetHiMt("dataDataSetHiMt","dataset for signal-pass template HiMt", RooArgSet(mass), Import( *fullTreeDataHiMtCut ) );
  RooDataHist dataDataHistHiMt("dataDataHistHiMt", "", RooArgSet(mass), dataDataSetHiMt, 1.0);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist sgnDataHist = dataSetProducer(fullTreeSgnCut, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdf("sgnTemplatePdf", "", RooArgSet(mass), sgnDataHist, 4);
  float numSgnScaled = getCorrectedEntries(fullTreeSgnCut, mass, ScaleFactorSgn);
  cout<<"SgnSumEntries "<<numSgnScaled<<endl;
  RooDataHist sgnDataHistTemp = dataSetProducer(fullTreeSgnCutTemp, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfTemp("sgnTemplatePdfTemp", "", RooArgSet(mass), sgnDataHistTemp, 4);
  float numSgnScaledTemp = getCorrectedEntries(fullTreeSgnCutTemp, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistTempJet = dataSetProducer(fullTreeSgnCutTempJet, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfTempJet("sgnTemplatePdfTempJet", "", RooArgSet(mass), sgnDataHistTempJet, 4);
  float numSgnScaledTempJet = getCorrectedEntries(fullTreeSgnCutTempJet, mass, ScaleFactorSgn);

  RooDataHist sgnDataHistSS = dataSetProducer(fullTreeSgnSSCut, mass, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMt = dataSetProducer(fullTreeSgnHiMtCut, mass, ScaleFactorSgn, -1.0);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHist = dataSetProducer(fullTreeWJetsCut, mass, ScaleFactorWJets, 1.0);
  RooHistPdf wjetsHistPdf("wjetsHistPdf", "", RooArgSet(mass), wjetsDataHist, 4);
  float numWJetsScaled = getCorrectedEntries(fullTreeWJetsCut, mass, ScaleFactorWJets,0,51);

  RooDataHist wjetsDataHistSS = dataSetProducer(fullTreeWJetsSSCut, mass, ScaleFactorWJets, -1.0);

  RooDataHist wjetsDataHistHiMt = dataSetProducer(fullTreeWJetsHiMtCut, mass, ScaleFactorWJets, -1.0);
  float numWJetsScaledHiMt = getCorrectedEntries(fullTreeWJetsHiMtCut, mass, ScaleFactorWJets,0,51);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zttDataHist = dataSetProducer(fullTreeZttCut, mass, ScaleFactorZtt, 1.0);
  RooHistPdf zttPdf("zttPdf", "", RooArgSet(mass), zttDataHist, 4);  
  float numZttScaled = getCorrectedEntries(fullTreeZttCut, mass, ScaleFactorZtt);

  RooDataHist zttDataHistSS = dataSetProducer(fullTreeZttSSCut, mass, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMt = dataSetProducer(fullTreeZttHiMtCut, mass, ScaleFactorZtt, -1.0);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHist = dataSetProducer(fullTreeTTJetsCut, mass, ScaleFactorTTJets, 1.0);
  RooHistPdf ttjetsPdf("ttjetsPdf", "", RooArgSet(mass), ttjetsDataHist,4);
  float numTTJetsScaled = getCorrectedEntries(fullTreeTTJetsCut, mass, ScaleFactorTTJets);

  RooDataHist ttjetsDataHistSS = dataSetProducer(fullTreeTTJetsSSCut, mass, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMt = dataSetProducer(fullTreeTTJetsHiMtCut, mass, ScaleFactorTTJets, -1.0);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHist = dataSetProducer(fullTreeWWCut, mass, ScaleFactorWW, 1.0);
  RooHistPdf wwPdf("wwPdf", "", RooArgSet(mass), wwDataHist,4);
  float numWWScaled = getCorrectedEntries(fullTreeWWCut, mass, ScaleFactorWW);

  mass.setBins( 50 );
  RooDataHist wwDataHistSS = dataSetProducer(fullTreeWWSSCut, mass, ScaleFactorWW, -1.0);

  mass.setBins( 50 );
  RooDataHist wwDataHistHiMt = dataSetProducer(fullTreeWWHiMtCut, mass, ScaleFactorWW, -1.0);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wzDataHist = dataSetProducer(fullTreeWZCut, mass, ScaleFactorWZ, 1.0);
  RooHistPdf wzPdf("wzPdf", "", RooArgSet(mass), wzDataHist,4);
  float numWZScaled = getCorrectedEntries(fullTreeWZCut, mass, ScaleFactorWZ);

  mass.setBins( 50 );
  RooDataHist wzDataHistSS = dataSetProducer(fullTreeWZSSCut, mass, ScaleFactorWZ, -1.0);

  mass.setBins( 50 );
  RooDataHist wzDataHistHiMt = dataSetProducer(fullTreeWZHiMtCut, mass, ScaleFactorWZ, -1.0);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zzDataHist = dataSetProducer(fullTreeZZCut, mass, ScaleFactorZZ, 1.0);
  RooHistPdf zzPdf("zzPdf", "", RooArgSet(mass), zzDataHist,4);
  float numZZScaled = getCorrectedEntries(fullTreeZZCut, mass, ScaleFactorZZ);

  mass.setBins( 50 );
  RooDataHist zzDataHistSS = dataSetProducer(fullTreeZZSSCut, mass, ScaleFactorZZ, -1.0);

  mass.setBins( 50 );
  RooDataHist zzDataHistHiMt = dataSetProducer(fullTreeZZHiMtCut, mass, ScaleFactorZZ, -1.0);

  //WJets from SideBand
  //RooPlot* massFrame1 = mass.frame();
  //dataDataHistHiMt.plotOn(massFrame1,MarkerColor(kRed));
  dataDataHistHiMt.add(sgnDataHistHiMt);
  dataDataHistHiMt.add(zttDataHistHiMt);
  dataDataHistHiMt.add(ttjetsDataHistHiMt);
  dataDataHistHiMt.add(wwDataHistHiMt);
  dataDataHistHiMt.add(wzDataHistHiMt);
  dataDataHistHiMt.add(zzDataHistHiMt);
  //dataDataHistHiMt.plotOn(massFrame1,MarkerColor(kBlue));
  float exWJetsP = numWJetsScaled/numWJetsScaledHiMt;
  cout<<"exWJetsP "<<exWJetsP<<endl;
  RooDataHist dataDataHistHiMtScaled("dataDataHistHiMtScaled", "", RooArgSet(mass), dataDataHistHiMt, exWJetsP);
  float numWJetsDataDriven = dataDataHistHiMtScaled.sumEntries();
  //dataDataHistHiMtScaled.plotOn(massFrame1,MarkerColor(kViolet));
  //massFrame1->Draw();
  RooHistPdf WJetsDataDrivenHistPdfP("WJetsDataDrivenHistPdfP", "", RooArgSet(mass), dataDataHistHiMtScaled, 4);

  // QCD
  dataDataHistSS.add(sgnDataHistSS);
  dataDataHistSS.add(wjetsDataHistSS);
  dataDataHistSS.add(ttjetsDataHistSS);
  dataDataHistSS.add(zttDataHistSS);
  dataDataHistSS.add(wwDataHistSS);
  dataDataHistSS.add(wzDataHistSS);
  dataDataHistSS.add(zzDataHistSS);
  RooHistPdf QCDHistPdfP("QCDHistPdfP", "", RooArgSet(mass), dataDataHistSS, 4);
  float numQCDDataDriven = dataDataHistSS.sumEntries();

  delete fullTreeSgnCut;
  delete fullTreeSgnCutTemp;
  delete fullTreeWJetsCut;
  delete fullTreeZttCut;
  delete fullTreeTTJetsCut;
  delete fullTreeWWCut;
  delete fullTreeWZCut;
  delete fullTreeZZCut;
  delete fullTreeDataCut;

  delete fullTreeSgnSSCut;
  delete fullTreeWJetsSSCut;
  delete fullTreeZttSSCut;
  delete fullTreeTTJetsSSCut;
  delete fullTreeWWSSCut;
  delete fullTreeWZSSCut;
  delete fullTreeZZSSCut;
  delete fullTreeDataSSCut;

  delete fullTreeSgnHiMtCut;
  delete fullTreeWJetsHiMtCut;
  delete fullTreeZttHiMtCut;
  delete fullTreeTTJetsHiMtCut;
  delete fullTreeWWHiMtCut;
  delete fullTreeWZHiMtCut;
  delete fullTreeZZHiMtCut;
  delete fullTreeDataHiMtCut;
  
  ///////////////////////////////////// Failing /////////////////////////////////////

  TFile *McP2 = new TFile("dummy2.root","RECREATE");

  std::cout<<"Processing Failing, OS, Low MT"<<std::endl;

  // Create trees with cuts: failing

  TTree* fullTreeSgnCutF = fullTreeSgn->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutFTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutFTempJet = fullTreeSgn->CopyTree( Form("(mcTrue == 0 && %s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  std::cout<<"Processing Failing, SS, Low MT"<<std::endl;

  //SS
  TTree* fullTreeDataSSCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCutF = fullTreeSgn->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

  std::cout<<"Processing Failing, OS, High MT"<<std::endl;

  //High MT
  TTree* fullTreeDataHiMtCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeSgnHiMtCutF = fullTreeSgn->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsHiMtCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttHiMtCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsHiMtCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWHiMtCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZHiMtCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZHiMtCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );

  //Normalization failing

  McP2->cd();

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist DataDataHistF = dataSetProducer(fullTreeDataCutF, mass, 1.0, 1.0);

  RooDataHist dataDataHistSSF = dataSetProducer(fullTreeDataSSCutF, mass, 1.0, 1.0);

  RooDataHist dataDataHistHiMtF = dataSetProducer(fullTreeDataHiMtCutF, mass, 1.0, 1.0);

  delete fullTreeDataCutF;
  delete fullTreeDataSSCutF;
  delete fullTreeDataHiMtCutF;

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist sgnDataHistF = dataSetProducer(fullTreeSgnCutF, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfF("sgnTemplatePdfF", "", RooArgSet(mass), sgnDataHistF, 4);
  float numSgnScaledF = getCorrectedEntries(fullTreeSgnCutF, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistFTemp = dataSetProducer(fullTreeSgnCutFTemp, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfFTemp("sgnTemplatePdfFTemp", "", RooArgSet(mass), sgnDataHistFTemp, 4);
  float numSgnScaledTempF = getCorrectedEntries(fullTreeSgnCutFTemp, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistFTempJet = dataSetProducer(fullTreeSgnCutFTempJet, mass, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfFTempJet("sgnTemplatePdfFTempJet", "", RooArgSet(mass), sgnDataHistFTempJet, 4);
  float numSgnScaledTempFJet = getCorrectedEntries(fullTreeSgnCutFTempJet, mass, ScaleFactorSgn);

  RooDataHist sgnDataHistSSF = dataSetProducer(fullTreeSgnSSCutF, mass, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMtF = dataSetProducer(fullTreeSgnHiMtCutF, mass, ScaleFactorSgn, -1.0);

  delete fullTreeSgnCutF;
  delete fullTreeSgnSSCutF;
  delete fullTreeSgnHiMtCutF;

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHistF = dataSetProducer(fullTreeWJetsCutF, mass, ScaleFactorWJets, 1.0);
  RooHistPdf wjetsHistPdfF("wjetsHistPdfF", "", RooArgSet(mass), wjetsDataHistF, 4);
  float numWJetsScaledF = getCorrectedEntries(fullTreeWJetsCutF, mass, ScaleFactorWJets,0,51);

  RooDataHist wjetsDataHistSSF = dataSetProducer(fullTreeWJetsSSCutF, mass, ScaleFactorWJets, -1.0);

  RooDataHist wjetsDataHistHiMtF = dataSetProducer(fullTreeWJetsHiMtCutF, mass, ScaleFactorWJets, -1.0);
  float numWJetsScaledFHiMt = getCorrectedEntries(fullTreeWJetsHiMtCutF, mass, ScaleFactorWJets,0,51);

  delete fullTreeWJetsCutF;
  delete fullTreeWJetsSSCutF;
  delete fullTreeWJetsHiMtCutF;

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zttDataHistF = dataSetProducer(fullTreeZttCutF, mass, ScaleFactorZtt, 1.0);
  RooHistPdf zttPdfF("zttPdfF", "", RooArgSet(mass), zttDataHistF, 4); 
  float numZttScaledF = getCorrectedEntries(fullTreeZttCutF, mass, ScaleFactorZtt);

  RooDataHist zttDataHistSSF = dataSetProducer(fullTreeZttSSCutF, mass, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMtF = dataSetProducer(fullTreeZttHiMtCutF, mass, ScaleFactorZtt, -1.0);

  delete fullTreeZttCutF;
  delete fullTreeZttSSCutF;
  delete fullTreeZttHiMtCutF;

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHistF = dataSetProducer(fullTreeTTJetsCutF, mass, ScaleFactorTTJets, 1.0);
  RooHistPdf ttjetsPdfF("ttjetsPdfF", "", RooArgSet(mass), ttjetsDataHistF, 4);
  float numTTJetsScaledF = getCorrectedEntries(fullTreeTTJetsCutF, mass, ScaleFactorTTJets);

  RooDataHist ttjetsDataHistSSF = dataSetProducer(fullTreeTTJetsSSCutF, mass, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMtF = dataSetProducer(fullTreeTTJetsHiMtCutF, mass, ScaleFactorTTJets, -1.0);

  delete fullTreeTTJetsCutF;
  delete fullTreeTTJetsSSCutF;
  delete fullTreeTTJetsHiMtCutF;

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHistF = dataSetProducer(fullTreeWWCutF, mass, ScaleFactorWW, 1.0);
  RooHistPdf wwPdfF("wwPdfF", "", RooArgSet(mass), wwDataHistF, 4);
  float numWWScaledF = getCorrectedEntries(fullTreeWWCutF, mass, ScaleFactorWW); 

  mass.setBins( 50 );
  RooDataHist wwDataHistSSF = dataSetProducer(fullTreeWWSSCutF, mass, ScaleFactorWW, -1.0);

  mass.setBins( 50 );
  RooDataHist wwDataHistHiMtF = dataSetProducer(fullTreeWWHiMtCutF, mass, ScaleFactorWW, -1.0);

  delete fullTreeWWCutF;
  delete fullTreeWWSSCutF;
  delete fullTreeWWHiMtCutF;

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wzDataHistF = dataSetProducer(fullTreeWZCutF, mass, ScaleFactorWZ, 1.0);
  RooHistPdf wzPdfF("wzPdfF", "", RooArgSet(mass), wzDataHistF, 4);
  float numWZScaledF = getCorrectedEntries(fullTreeWZCutF, mass, ScaleFactorWZ);

  mass.setBins( 50 );
  RooDataHist wzDataHistSSF = dataSetProducer(fullTreeWZSSCutF, mass, ScaleFactorWZ, -1.0);

  mass.setBins( 50 );
  RooDataHist wzDataHistHiMtF = dataSetProducer(fullTreeWZHiMtCutF, mass, ScaleFactorWZ, -1.0);

  delete fullTreeWZCutF;
  delete fullTreeWZSSCutF;
  delete fullTreeWZHiMtCutF;

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zzDataHistF = dataSetProducer(fullTreeZZCutF, mass, ScaleFactorZZ, 1.0);
  RooHistPdf zzPdfF("zzPdfF", "", RooArgSet(mass), zzDataHistF, 4);
  float numZZScaledF = getCorrectedEntries(fullTreeZZCutF, mass, ScaleFactorZZ);

  mass.setBins( 50 );
  RooDataHist zzDataHistSSF = dataSetProducer(fullTreeZZSSCutF, mass, ScaleFactorZZ, -1.0);

  mass.setBins( 50 );
  RooDataHist zzDataHistHiMtF = dataSetProducer(fullTreeZZHiMtCutF, mass, ScaleFactorZZ, -1.0);

  delete fullTreeZZCutF;
  delete fullTreeZZSSCutF;
  delete fullTreeZZHiMtCutF;

  //WJets from SideBand
  //RooPlot* massFrame2 = mass.frame();
  //dataDataHistHiMtF.plotOn(massFrame2,MarkerColor(kRed));
  dataDataHistHiMtF.add(sgnDataHistHiMtF);
  dataDataHistHiMtF.add(zttDataHistHiMtF);
  dataDataHistHiMtF.add(ttjetsDataHistHiMtF);
  dataDataHistHiMtF.add(wwDataHistHiMtF);
  dataDataHistHiMtF.add(wzDataHistHiMtF);
  dataDataHistHiMtF.add(zzDataHistHiMtF);
  //dataDataHistHiMtF.plotOn(massFrame2,MarkerColor(kBlue));
  float exWJetsF = numWJetsScaledF/numWJetsScaledFHiMt;
  cout<<"exWJetsF "<<exWJetsF<<endl;
  RooDataHist dataDataHistHiMtScaledF("dataDataHistHiMtScaledF", "", RooArgSet(mass), dataDataHistHiMtF, exWJetsF);
  float numWJetsDataDrivenF = dataDataHistHiMtScaledF.sumEntries(); 
  //dataDataHistHiMtScaledF.plotOn(massFrame2,MarkerColor(kViolet));
  RooHistPdf WJetsDataDrivenHistPdfF("WJetsDataDrivenHistPdfF", "", RooArgSet(mass), dataDataHistHiMtScaledF, 4); //Shape WJets Fail da Sideband

  // QCD
  dataDataHistSSF.add(sgnDataHistSSF);
  dataDataHistSSF.add(wjetsDataHistSSF);
  dataDataHistSSF.add(ttjetsDataHistSSF);
  dataDataHistSSF.add(zttDataHistSSF);
  dataDataHistSSF.add(wwDataHistSSF);
  dataDataHistSSF.add(wzDataHistSSF);
  dataDataHistSSF.add(zzDataHistSSF);

  RooHistPdf QCDHistPdfF("QCDHistPdfF", "", RooArgSet(mass), dataDataHistSSF, 4);
  float numQCDDataDrivenF = dataDataHistSSF.sumEntries(); 

  //////////////////////////////////////////////////////////////// Datasets & Pdfs //////////////////////////////////////////////////////////////// 

  // Constraints 
  // -----------------------------------------
  RooRealVar f("f","f",0.5,0.,1.) ;
  RooRealVar fW("fW","fW",0.5,0.,1.) ;
  // Construct Gaussian constraint p.d.f on parameter f at 1.05 with resolution of 0.1
  RooGaussian fconstraint("fconstraint","fconstraint",f,RooConst(1.05),RooConst(0.1));
  RooGaussian fconstraintW("fconstraintW","fconstraintW",fW,RooConst(1.05),RooConst(0.1));

  // Multiply constraint term with regular p.d.f using RooProdPdf
  // Specify in fitTo() that internal constraints on parameter f should be used

  // Multiply constraint with p.d.f
  RooProdPdf QCDHistPdfF_C("QCDHistPdfF_C","QCD model with constraint",RooArgSet(QCDHistPdfF,fconstraint));
  RooProdPdf QCDHistPdfP_C("QCDHistPdfP_C","QCD model with constraint",RooArgSet(QCDHistPdfP,fconstraint));
  RooProdPdf WJetsHistPdfF_C("WJetsHistPdfF_C","WJets model with constraint",RooArgSet(WJetsDataDrivenHistPdfF,fconstraintW));
  RooProdPdf WJetsHistPdfP_C("WJetsHistPdfP_C","WJets model with constraint",RooArgSet(WJetsDataDrivenHistPdfP,fconstraintW));

  //////////////////////////////////////////////////////////////// SimFit ////////////////////////////////////////////////////////////////

  //Create pdfs for passing
  RooRealVar CoeffSgnP("CoeffSgnP","",0,10000000);
  RooRealVar CoeffSgnJetP("CoeffSgnJetP","",numSgnScaledTempJet,0.8*numSgnScaledTempJet,1.2*numSgnScaledTempJet);
  RooRealVar CoeffZttP("CoeffZttP","",numZttScaled,0,10000000);
  RooRealVar CoeffWJetsP("CoeffWJetsP","",numWJetsDataDriven,0,10000000);
  RooRealVar CoeffTTJetsP("CoeffTTJetsP","",numTTJetsScaled,0,10000000);
  RooRealVar CoeffQCDP("CoeffQCDP","",numQCDDataDriven,0,10000000);
  RooRealVar CoeffWWP("CoeffWWP","",numWWScaled,0,10000000);
  RooRealVar CoeffWZP("CoeffWZP","",numWZScaled,0,10000000);
  RooRealVar CoeffZZP("CoeffZZP","",numZZScaled,0,10000000);

  RooAddPdf DataModelP("DataModelP", "", RooArgList(sgnTemplatePdfTemp,sgnTemplatePdfTempJet,zttPdf,WJetsHistPdfP_C,ttjetsPdf,QCDHistPdfP_C,wwPdf,wzPdf,zzPdf), RooArgList(CoeffSgnP/*DataNumSgnP*/,CoeffSgnJetP,CoeffZttP,CoeffWJetsP,CoeffTTJetsP,CoeffQCDP,CoeffWWP,CoeffWZP,CoeffZZP));

  float errWZF = getError(fullTreeWZCutF, mass, ScaleFactorWZ, 1.0);
  float errZZF = getError(fullTreeZZCutF, mass, ScaleFactorZZ, 1.0);

  //Create pdfs for failing
  RooRealVar CoeffSgnF("CoeffSgnF","",0,10000000);
  RooRealVar CoeffSgnJetF("CoeffSgnJetF","",numSgnScaledTempFJet,0.8*numSgnScaledTempFJet,1.2*numSgnScaledTempFJet);
  RooRealVar CoeffZttF("CoeffZttF","",numZttScaledF,0,10000000);
  RooRealVar CoeffWJetsF("CoeffWJetsF","",numWJetsDataDrivenF,0,10000000);
  RooRealVar CoeffTTJetsF("CoeffTTJetsF","",numTTJetsScaledF,0,10000000);
  RooRealVar CoeffQCDF("CoeffQCDF","",numQCDDataDrivenF,0,10000000);
  RooRealVar CoeffWWF("CoeffWWF","",numWWScaledF,0,10000000);
  RooRealVar CoeffWZF("CoeffWZF","",numWZScaledF,numWZScaledF-2*errWZF,numWZScaledF+2*errWZF);
  RooRealVar CoeffZZF("CoeffZZF","",numZZScaledF,numZZScaledF-2*errZZF,numZZScaledF+2*errZZF);

  RooAddPdf DataModelF("DataModelF", "", RooArgList(sgnTemplatePdfFTemp,sgnTemplatePdfFTempJet,zttPdfF,wjetsHistPdfF,ttjetsPdfF,QCDHistPdfF_C,wwPdfF,wzPdfF,zzPdfF), RooArgList(CoeffSgnF/*DataNumSgnP*/,CoeffSgnJetF,CoeffZttF,CoeffWJetsF,CoeffTTJetsF,CoeffQCDF,CoeffWWF,CoeffWZF,CoeffZZF));

  mass.setBins(nBins_);

  RooCategory category("category","category") ;
  category.defineType("pass") ;
  category.defineType("fail") ;

  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataHistP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataHistF.createHistogram("histoDataF",mass))), Weight(0.5) ) ;
  // unbinned combined dataset
  //RooDataSet DataCombDataUnBinned("DataCombDataUnBinned","combined data",mass,Index(category),Import("pass", DataDataSetP) ,Import("fail",DataDataSetF), Weight(0.5) ) ;

  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;

  //mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
  RooFitResult* ResDataCombinedFit =  0;
  if(doBinned_)  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4), /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/  SumW2Error(1));
  //else ResDataCombinedFit = DataSimPdf.fitTo(DataCombDataUnBinned , Extended(1), Minos(1), Save(1), NumCPU(4),  /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/ SumW2Error(1));

  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());

  string theSample = "Data";

  RooPlot* DataFrameP = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012 #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}: passing probe",theSample.c_str(),Lumi_)));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"),Name("dataP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlack),Name("modelP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnTemplatePdfTemp"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnTemplatePdfTempJet"), LineColor(kBlue-10), LineStyle(kSolid),Name("signal onlyP (jet)"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("zttPdf"), LineColor(kRed), LineStyle(kSolid),Name("ZttP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("WJetsDataDrivenHistPdfP"), LineColor(kGreen), LineStyle(kSolid),Name("WJetsP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("ttjetsPdf"), LineColor(kMagenta), LineStyle(kSolid),Name("TTJetsP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("QCDHistPdfP"), LineColor(kOrange), LineStyle(kSolid),Name("QCDP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("wwPdf"), LineColor(kViolet), LineStyle(kSolid),Name("WWP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("wzPdf"), LineColor(kCyan), LineStyle(kSolid),Name("WZP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("zzPdf"), LineColor(kYellow), LineStyle(kSolid),Name("ZZP"));
  DataFrameP->SetTitleOffset(1,"Y");
  DataFrameP->SetTitleSize(1,"Y");
  float binWidth = (xHigh_ - xLow_) / nBins_;
  DataFrameP->GetYaxis()->SetTitle(Form("Entries/%.0f GeV/c^{2}", binWidth));

  RooPlot* DataFrameF = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012 #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}: failing probe",theSample.c_str(),Lumi_)));
  DataCombData.plotOn(DataFrameF,Cut("category==category::fail"),Name("dataF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), LineColor(kBlack),Name("modelF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnTemplatePdfFTemp"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnTemplatePdfFTempJet"), LineColor(kBlue-10), LineStyle(kSolid),Name("signal onlyF (jet)"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("zttPdfF"), LineColor(kRed), LineStyle(kSolid),Name("ZttF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("wjetsHistPdfF"), LineColor(kGreen), LineStyle(kSolid),Name("WJetsF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("ttjetsPdfF"), LineColor(kMagenta), LineStyle(kSolid),Name("TTJetsF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("QCDHistPdfF"), LineColor(kOrange), LineStyle(kSolid),Name("QCDF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("wwPdfF"), LineColor(kViolet), LineStyle(kSolid),Name("WWF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("wzPdfF"), LineColor(kCyan), LineStyle(kSolid),Name("WZF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("zzPdfF"), LineColor(kYellow), LineStyle(kSolid),Name("ZZF"));
  DataFrameF->SetTitleOffset(1,"Y");
  DataFrameF->SetTitleSize(1,"Y");
  DataFrameF->GetYaxis()->SetTitle(Form("Entries/%.0f GeV/c^{2}", binWidth));

  TCanvas *cPass = new TCanvas("fitCanvasP","canvas",10,30,650,600);
  cPass->SetGrid(0,0);
  cPass->SetFillStyle(4000);
  cPass->SetFillColor(10);
  cPass->SetTicky();
  cPass->SetObjectStat(0);

  cPass->cd();
  DataFrameP->Draw();
  TLegend *leg1 = new TLegend(0.55,0.55,0.8,0.8);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry("dataP","Data", "P");
  leg1->AddEntry("modelP","Signal + bkg","L");
  //leg1->AddEntry("backgroundP","Background only", "L");
  leg1->AddEntry("signal onlyP","Signal (#mu#rightarrow#tau)", "L");
  leg1->AddEntry("signal onlyP (jet)","Signal (jet#rightarrow#tau)", "L");
  leg1->AddEntry("ZttP","Z#rightarrow#tau#tau", "L");
  leg1->AddEntry("WJetsP","WJets", "L");
  leg1->AddEntry("TTJetsP","TTJets", "L");
  leg1->AddEntry("QCDP","QCD", "L");
  leg1->AddEntry("WWP","WW", "L");
  leg1->AddEntry("WZP","WZ", "L");
  leg1->AddEntry("ZZP","ZZ", "L");
  leg1->Draw();
  //return;

  string fileNameP = "fitCanvasPassMuToTau_"+tnp_+"_"+category_;
  cPass->SaveAs(Form("%s_%.2f.png",fileNameP.c_str(), binCenter_));

  TCanvas *cFail = new TCanvas("fitCanvasF","canvas",10,30,650,600);
  cFail->SetGrid(0,0);
  cFail->SetFillStyle(4000);
  cFail->SetFillColor(10);
  cFail->SetTicky();
  cFail->SetObjectStat(0);

  cFail->cd();
  DataFrameF->Draw();
  TLegend *leg2 = new TLegend(0.55,0.55,0.8,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);
  leg2->AddEntry("dataF","Data", "P");
  leg2->AddEntry("modelF","Signal + bkg","L");
  //leg2->AddEntry("backgroundF","Background only", "L");
  leg2->AddEntry("signal onlyF","Signal (#mu#rightarrow#tau)", "L");
  leg2->AddEntry("signal onlyF (jet)","Signal (jet#rightarrow#tau)", "L");
  leg2->AddEntry("ZttF","Z#rightarrow#tau#tau", "L");
  leg2->AddEntry("WJetsF","WJets", "L");
  leg2->AddEntry("TTJetsF","TTJets", "L");
  leg2->AddEntry("QCDF","QCD", "L");
  leg2->AddEntry("WWF","WW", "L");
  leg2->AddEntry("WZF","WZ", "L");
  leg2->AddEntry("ZZF","ZZ", "L");
  leg2->Draw();

  string fileNameF = "fitCanvasFailMuToTau_"+tnp_+"_"+category_;
  cFail->SaveAs(Form("%s_%.2f.png",fileNameF.c_str(), binCenter_));

  ResDataCombinedFit->printArgs(std::cout);
  cout << endl;
  ResDataCombinedFit->printValue(std::cout);
  cout << endl;

  RooFormulaVar DataEffFit("DataEffFit","DataEffFit","CoeffSgnP/(CoeffSgnP+CoeffSgnF)", RooArgList(CoeffSgnP,CoeffSgnP,CoeffSgnF));
  float funcError = DataEffFit.getPropagatedError(*ResDataCombinedFit) ; 

  std::vector<float> results;
  results.push_back(CoeffSgnP.getVal());
  results.push_back(CoeffSgnF.getVal());
  results.push_back(CoeffSgnP.getError());
  results.push_back(CoeffSgnF.getError());
  results.push_back(DataEffFit.getVal());
  results.push_back(funcError);
  results.push_back(McTruthEff);
  results.push_back(BinomialError);

  cout<<"Sgn passing : "<< CoeffSgnP.getVal()<<" +/- "<<results[2]<<" Sgn failing : "<<CoeffSgnF.getVal()<< " +/- "<<results[3]<<endl;
  cout<<"Efficiency : "<< CoeffSgnP.getVal()/(CoeffSgnP.getVal()+CoeffSgnF.getVal()) << " +/- " <<results[5]<<endl;

  string TreeName = "treeForSaving_"+tnp_+"_"+category_;
  TFile *Save = new TFile(Form("%s_%.2f.root",TreeName.c_str(), binCenter_),"RECREATE");
  ResDataCombinedFit->Write();
  DataFrameP->Write();
  DataFrameF->Write();
  Save->Close();
  string SaveFileName = "fileForSaving_"+tnp_+"_"+category_;
  DataFitParam.writeToFile(Form("%s_%.2f.txt",SaveFileName.c_str(), binCenter_));

  std::ofstream outputSys(Form("%s_%.2f.txt",SaveFileName.c_str(), binCenter_), ios::app);
  outputSys<<"\n"<< bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError <<endl;
  outputSys<<"Sgn passing : "<< CoeffSgnP.getVal()<<" +/- "<<results[2]<<" Sgn failing : "<<CoeffSgnF.getVal()<< " +/- "<<results[3]<<endl;
  outputSys<<"Efficiency : "<< CoeffSgnP.getVal()/(CoeffSgnP.getVal()+CoeffSgnF.getVal()) << " +/- " <<results[5]<<endl;
  float eff = CoeffSgnP.getVal()/(CoeffSgnP.getVal()+CoeffSgnF.getVal());
  float DataMC = eff/McTruthEff;
  float errDataMC = sqrt(results[5]*results[5] + DataMC*DataMC*BinomialError*BinomialError)/McTruthEff;
  cout<<"Data/MC : "<< DataMC << " +/- " <<errDataMC<<endl;
  outputSys<<"Data/MC : "<< DataMC << " +/- " <<errDataMC<<endl;
  outputSys.close();

}

void calcutateFit(){

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7",50,2.0,0.3);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7",50,2.0,0.3);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7",50,2.0,0.3);

	cout<<"passingIsoLooseMuonVetoLoose2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta>1.7",50,2.0,0.3);

	cout<<"passingIsoLooseMuonVetoMedium2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta>1.7",50,2.0,0.3);

	cout<<"passingIsoLooseMuonVetoTight2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta<1.2",50,0.6,0.6);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta>1.2 && abseta<1.7",50,1.45,0.25);
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta>1.7",50,2.0,0.3);

}

