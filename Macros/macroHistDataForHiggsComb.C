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

RooDataHist dataSetProducer(TTree * fullTreeSgnCut, RooRealVar mass, float ScaleFactorSgn, float weight, std::string name = ""){

	TH1F* hMass = new TH1F("hMass","",50,70,120); 
  	if(ScaleFactorSgn != 1) {
  		fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("mass>>hMass");
        cout<<name.c_str()<<" "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	TH1F * hMass2 = (TH1F*)hMass->Clone();

	for(int i=1; i<hMass->GetSize()-1; i++){

		int val = (int)hMass->GetBinContent(i);
		hMass2->SetBinContent(i,val);

	}

	if(name != "") {

		hMass->Write(name.c_str());
		hMass2->Write((name + "Int").c_str());

	}

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
	std::string binCenter_           = "Barrel",
	double xLow_                     = 70,
	double xHigh_                    = 120,
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	//const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_muPFIsolation < 0.1",
	//const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_muPFIsolation < 0.1",
	//const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_muPFIsolation < 0.1"
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"
	){
  
  /*TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);*/

  string path = "./Input2014/";
  string pathData = "./Input2014/";

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

  fullTreeSgn->Draw("mass>>hS",Form("tag_puMCWeightRun2012*(%s && mass>%f && mass<%f && mcTrue && pair_charge==0 && tag_Mt < 40 && tag_muTriggerMatching && tag_triggerBitSingleMu > 0.5 )",bin_.c_str(),xLow_,xHigh_));
  double SGNtrue = hS->Integral(0,2);
  fullTreeSgn->Draw("mass>>hSP",Form("tag_puMCWeightRun2012*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && pair_charge==0 && tag_Mt < 40 && tag_muTriggerMatching && tag_triggerBitSingleMu > 0.5 )",bin_.c_str(),category_.c_str(),cutValue_,xLow_,xHigh_));
  double SGNtruePass = hSP->Integral(0,2);

  cout<<SGNtruePass<<" "<<SGNtrue<<endl;
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

  McP->cd();

  //Normalization passing

  //float Lumi_ = 19484.55;
  float Lumi_ = 15829.435;

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

  // mass variable

  RooRealVar mass("mass","m_{tp} (GeV/c^{2})",xLow_,xHigh_);
  mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );

  string TreeNameP = "./HistosForHiggsCombine_2014/histForSaving_"+tnp_+"_"+category_;
  TFile *SaveP = new TFile(Form("%s%s.root",TreeNameP.c_str(), binCenter_.c_str()),"RECREATE");

  SaveP->mkdir("pass");
  SaveP->cd("pass");

  TH1F* hDummy1 = new TH1F("hDummy","",50,70,120); 
  hDummy1->SetBinContent(1,1);
  hDummy1->Write("signal");

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  //Roofit datasets of data passing and failing cut

  RooDataSet DataDataSetP("DataDataSetP", "dataset for Data pass", RooArgSet(mass), Import( *fullTreeDataCut ) );
  //std::cout << "data dataset Pass " << DataDataSetP.numEntries() << "  " << std::endl;
  RooDataHist DataDataHistP("DataDataHistP", "", RooArgSet(mass), DataDataSetP, 1.0);
  TH1F * dataP = (TH1F*)DataDataHistP.createHistogram("mass",50);
  dataP->Write("data_obs");

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
  //float numSgnScaled = getCorrectedEntries(fullTreeSgnCut, mass, ScaleFactorSgn);
  //cout<<"SgnSumEntries "<<numSgnScaled<<endl;
  RooDataHist sgnDataHistTemp = dataSetProducer(fullTreeSgnCutTemp, mass, ScaleFactorSgn, 1.0, "zmumu");
  RooHistPdf sgnTemplatePdfTemp("sgnTemplatePdfTemp", "", RooArgSet(mass), sgnDataHistTemp, 4);
  //float numSgnScaledTemp = getCorrectedEntries(fullTreeSgnCutTemp, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistTempJet = dataSetProducer(fullTreeSgnCutTempJet, mass, ScaleFactorSgn, 1.0, "zmumuJet");
  RooHistPdf sgnTemplatePdfTempJet("sgnTemplatePdfTempJet", "", RooArgSet(mass), sgnDataHistTempJet, 4);
  //float numSgnScaledTempJet = getCorrectedEntries(fullTreeSgnCutTempJet, mass, ScaleFactorSgn);

  RooDataHist sgnDataHistSS = dataSetProducer(fullTreeSgnSSCut, mass, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMt = dataSetProducer(fullTreeSgnHiMtCut, mass, ScaleFactorSgn, -1.0);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHist = dataSetProducer(fullTreeWJetsCut, mass, ScaleFactorWJets, 1.0, "wjets");
  RooHistPdf wjetsHistPdf("wjetsHistPdf", "", RooArgSet(mass), wjetsDataHist, 4);
  float numWJetsScaled = getCorrectedEntries(fullTreeWJetsCut, mass, ScaleFactorWJets,0,51);
  TH1F * wjets = (TH1F*)wjetsDataHist.createHistogram("mass",50);

  RooDataHist wjetsDataHistSS = dataSetProducer(fullTreeWJetsSSCut, mass, ScaleFactorWJets, -1.0);

  RooDataHist wjetsDataHistHiMt = dataSetProducer(fullTreeWJetsHiMtCut, mass, ScaleFactorWJets, -1.0);
  float numWJetsScaledHiMt = getCorrectedEntries(fullTreeWJetsHiMtCut, mass, ScaleFactorWJets,0,51);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zttDataHist = dataSetProducer(fullTreeZttCut, mass, ScaleFactorZtt, 1.0, "ztautau");
  RooHistPdf zttPdf("zttPdf", "", RooArgSet(mass), zttDataHist, 4);  
  //float numZttScaled = getCorrectedEntries(fullTreeZttCut, mass, ScaleFactorZtt);

  RooDataHist zttDataHistSS = dataSetProducer(fullTreeZttSSCut, mass, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMt = dataSetProducer(fullTreeZttHiMtCut, mass, ScaleFactorZtt, -1.0);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHist = dataSetProducer(fullTreeTTJetsCut, mass, ScaleFactorTTJets, 1.0, "ttjets");
  RooHistPdf ttjetsPdf("ttjetsPdf", "", RooArgSet(mass), ttjetsDataHist,4);
  //float numTTJetsScaled = getCorrectedEntries(fullTreeTTJetsCut, mass, ScaleFactorTTJets);

  RooDataHist ttjetsDataHistSS = dataSetProducer(fullTreeTTJetsSSCut, mass, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMt = dataSetProducer(fullTreeTTJetsHiMtCut, mass, ScaleFactorTTJets, -1.0);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHist = dataSetProducer(fullTreeWWCut, mass, ScaleFactorWW, 1.0, "ww");
  RooHistPdf wwPdf("wwPdf", "", RooArgSet(mass), wwDataHist,4);
  //float numWWScaled = getCorrectedEntries(fullTreeWWCut, mass, ScaleFactorWW);

  mass.setBins( 50 );
  RooDataHist wwDataHistSS = dataSetProducer(fullTreeWWSSCut, mass, ScaleFactorWW, -1.0);

  mass.setBins( 50 );
  RooDataHist wwDataHistHiMt = dataSetProducer(fullTreeWWHiMtCut, mass, ScaleFactorWW, -1.0);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wzDataHist = dataSetProducer(fullTreeWZCut, mass, ScaleFactorWZ, 1.0, "wz");
  RooHistPdf wzPdf("wzPdf", "", RooArgSet(mass), wzDataHist,4);
  //float numWZScaled = getCorrectedEntries(fullTreeWZCut, mass, ScaleFactorWZ);

  mass.setBins( 50 );
  RooDataHist wzDataHistSS = dataSetProducer(fullTreeWZSSCut, mass, ScaleFactorWZ, -1.0);

  mass.setBins( 50 );
  RooDataHist wzDataHistHiMt = dataSetProducer(fullTreeWZHiMtCut, mass, ScaleFactorWZ, -1.0);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zzDataHist = dataSetProducer(fullTreeZZCut, mass, ScaleFactorZZ, 1.0, "zz");
  RooHistPdf zzPdf("zzPdf", "", RooArgSet(mass), zzDataHist,4);
  //float numZZScaled = getCorrectedEntries(fullTreeZZCut, mass, ScaleFactorZZ);

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
  TH1F * wjetsDataDrivenP = (TH1F*)dataDataHistHiMtScaled.createHistogram("mass",50);
  wjetsDataDrivenP->Write("wjetsDataDriven");

  float numWJetsDataDriven = dataDataHistHiMtScaled.sumEntries();
  float normWJetsMC = numWJetsDataDriven / wjets->Integral(1,50);
  TH1F * wjetsNorm = (TH1F*)wjets->Clone();
  wjetsNorm->Scale(normWJetsMC);
  wjetsNorm->Write("wjetsNorm");

  TH1F * wjetsDataDrivenP2 = (TH1F*)wjetsDataDrivenP->Clone();

  for(int i=1; i<wjetsDataDrivenP->GetSize()-1; i++){

		int val = (int)wjetsDataDrivenP->GetBinContent(i);
		wjetsDataDrivenP2->SetBinContent(i,val);

  }

  wjetsDataDrivenP2->Write("wjetsDataDrivenInt");

  //float numWJetsDataDriven = dataDataHistHiMtScaled.sumEntries();
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
  TH1F * qcdDataDrivenP = (TH1F*)dataDataHistSS.createHistogram("mass",50);
  qcdDataDrivenP->Write("qcdDataDriven");
  RooHistPdf QCDHistPdfP("QCDHistPdfP", "", RooArgSet(mass), dataDataHistSS, 4);
  //float numQCDDataDriven = dataDataHistSS.sumEntries();

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

  McP->Close();
  SaveP->Close();
  
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

  //string TreeNameF = "./HistosForHiggsCombine/histForSaving_"+tnp_+"_"+category_ + "_fail";
  TFile *SaveF = new TFile(Form("%s%s.root",TreeNameP.c_str(), binCenter_.c_str()),"UPDATE");

  SaveF->mkdir("fail");
  SaveF->cd("fail");

  TH1F* hDummy2 = new TH1F("hDummy","",50,70,120); 
  hDummy2->SetBinContent(1,1);
  hDummy2->Write("signal");

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist DataDataHistF = dataSetProducer(fullTreeDataCutF, mass, 1.0, 1.0);
  TH1F * dataF = (TH1F*)DataDataHistF.createHistogram("mass",50);
  dataF->Write("data_obs");

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
  //float numSgnScaledF = getCorrectedEntries(fullTreeSgnCutF, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistFTemp = dataSetProducer(fullTreeSgnCutFTemp, mass, ScaleFactorSgn, 1.0,"zmumu");
  RooHistPdf sgnTemplatePdfFTemp("sgnTemplatePdfFTemp", "", RooArgSet(mass), sgnDataHistFTemp, 4);
  //float numSgnScaledTempF = getCorrectedEntries(fullTreeSgnCutFTemp, mass, ScaleFactorSgn);
  RooDataHist sgnDataHistFTempJet = dataSetProducer(fullTreeSgnCutFTempJet, mass, ScaleFactorSgn, 1.0,"zmumuJet");
  RooHistPdf sgnTemplatePdfFTempJet("sgnTemplatePdfFTempJet", "", RooArgSet(mass), sgnDataHistFTempJet, 4);
  //float numSgnScaledTempFJet = getCorrectedEntries(fullTreeSgnCutFTempJet, mass, ScaleFactorSgn);

  RooDataHist sgnDataHistSSF = dataSetProducer(fullTreeSgnSSCutF, mass, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMtF = dataSetProducer(fullTreeSgnHiMtCutF, mass, ScaleFactorSgn, -1.0);

  delete fullTreeSgnCutF;
  delete fullTreeSgnSSCutF;
  delete fullTreeSgnHiMtCutF;

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHistF = dataSetProducer(fullTreeWJetsCutF, mass, ScaleFactorWJets, 1.0,"wjets");
  RooHistPdf wjetsHistPdfF("wjetsHistPdfF", "", RooArgSet(mass), wjetsDataHistF, 4);
  float numWJetsScaledF = getCorrectedEntries(fullTreeWJetsCutF, mass, ScaleFactorWJets,0,51);
  TH1F * wjetsF = (TH1F*)wjetsDataHistF.createHistogram("mass",50);

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
  RooDataHist zttDataHistF = dataSetProducer(fullTreeZttCutF, mass, ScaleFactorZtt, 1.0,"ztautau");
  RooHistPdf zttPdfF("zttPdfF", "", RooArgSet(mass), zttDataHistF, 4); 
  //float numZttScaledF = getCorrectedEntries(fullTreeZttCutF, mass, ScaleFactorZtt);

  RooDataHist zttDataHistSSF = dataSetProducer(fullTreeZttSSCutF, mass, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMtF = dataSetProducer(fullTreeZttHiMtCutF, mass, ScaleFactorZtt, -1.0);

  delete fullTreeZttCutF;
  delete fullTreeZttSSCutF;
  delete fullTreeZttHiMtCutF;

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHistF = dataSetProducer(fullTreeTTJetsCutF, mass, ScaleFactorTTJets, 1.0,"ttjets");
  RooHistPdf ttjetsPdfF("ttjetsPdfF", "", RooArgSet(mass), ttjetsDataHistF, 4);
  //float numTTJetsScaledF = getCorrectedEntries(fullTreeTTJetsCutF, mass, ScaleFactorTTJets);

  RooDataHist ttjetsDataHistSSF = dataSetProducer(fullTreeTTJetsSSCutF, mass, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMtF = dataSetProducer(fullTreeTTJetsHiMtCutF, mass, ScaleFactorTTJets, -1.0);

  delete fullTreeTTJetsCutF;
  delete fullTreeTTJetsSSCutF;
  delete fullTreeTTJetsHiMtCutF;

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHistF = dataSetProducer(fullTreeWWCutF, mass, ScaleFactorWW, 1.0,"ww");
  RooHistPdf wwPdfF("wwPdfF", "", RooArgSet(mass), wwDataHistF, 4);
  //float numWWScaledF = getCorrectedEntries(fullTreeWWCutF, mass, ScaleFactorWW); 

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
  RooDataHist wzDataHistF = dataSetProducer(fullTreeWZCutF, mass, ScaleFactorWZ, 1.0,"wz");
  RooHistPdf wzPdfF("wzPdfF", "", RooArgSet(mass), wzDataHistF, 4);
  //float numWZScaledF = getCorrectedEntries(fullTreeWZCutF, mass, ScaleFactorWZ);

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
  RooDataHist zzDataHistF = dataSetProducer(fullTreeZZCutF, mass, ScaleFactorZZ, 1.0,"zz");
  RooHistPdf zzPdfF("zzPdfF", "", RooArgSet(mass), zzDataHistF, 4);
  //float numZZScaledF = getCorrectedEntries(fullTreeZZCutF, mass, ScaleFactorZZ);

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
  TH1F * wjetsDataDrivenF = (TH1F*)dataDataHistHiMtScaledF.createHistogram("mass",50);
  wjetsDataDrivenF->Write("wjetsDataDriven");

  float numWJetsDataDrivenF = dataDataHistHiMtScaledF.sumEntries();
  float normWJetsMCF = numWJetsDataDrivenF / wjetsF->Integral(1,50);
  TH1F * wjetsNormF = (TH1F*)wjetsF->Clone();
  wjetsNormF->Scale(normWJetsMCF);
  wjetsNormF->Write("wjetsNorm");

  TH1F * wjetsDataDrivenF2 = (TH1F*)wjetsDataDrivenF->Clone();

  for(int i=1; i<wjetsDataDrivenF->GetSize()-1; i++){

		int val = (int)wjetsDataDrivenF->GetBinContent(i);
		wjetsDataDrivenF2->SetBinContent(i,val);

  }

  wjetsDataDrivenF2->Write("wjetsDataDrivenInt");

  //float numWJetsDataDrivenF = dataDataHistHiMtScaledF.sumEntries(); 
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
  TH1F * qcdDataDrivenF = (TH1F*)dataDataHistSSF.createHistogram("mass",50);
  qcdDataDrivenF->Write("qcdDataDriven");
  float numQCDDataDrivenF = dataDataHistSSF.sumEntries(); 
  std::cout<<"QCDFail "<<numQCDDataDrivenF<<std::endl;

  SaveF->Close();

}

void calculateFit(){

	cout<<"passingIsoLooseMuonVeto_LooseV1"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2",50,"Barrel");
	/*fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_MediumV1"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_TightV1"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_LooseV2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose2","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_MediumV2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMedium2","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_TightV2"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight2","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_LooseV3"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose3","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose3","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLoose3","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_TightV3"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight3","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight3","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTight3","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_LooseMVA"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLooseMVA","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLooseMVA","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoLooseMVA","abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_MediumMVA"<<endl;
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMediumMVA","abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMediumMVA","abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoMediumMVA","abseta>1.7",50,"Endcap");*/

	//cout<<"passingIsoLooseMuonVeto_TightMVA"<<endl;
	//fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTightMVA","abseta<1.2",50,"Barrel");
	//fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTightMVA","abseta>1.2 && abseta<1.7",50,"Overlap");
	//fitStudyTemplatesFromMC("muToTau","passingIsoLooseMuonVetoTightMVA","abseta>1.7",50,"Endcap");

}

