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

#include "TFile.h"
#include "TKey.h"
#include "TStyle.h"
#include "TClass.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TGaxis.h"
#include <TStyle.h>
#include "TText.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include <string>
#include <fstream>
#include <iostream>

#include <vector>
#include <sstream>
#include "TGraphAsymmErrors.h"

using namespace std;
using namespace RooFit;

RooDataHist dataSetProducer(TTree * fullTreeSgnCut, RooRealVar mass, float NumSgnP, float ScaleFactorSgn, float weight){

  	//RooDataSet sgnDataSet("sgnDataSet", "dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  	//if(NumSgnP*ScaleFactorSgn != 1) sgnDataSet.reduce(EventRange(1,(int)( NumSgnP*ScaleFactorSgn)));
	TH1F* hMass = new TH1F("hMass","",50,70,120); 
  	fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
	hMass->Scale(ScaleFactorSgn);
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;
  	RooDataHist sgnDataHist("sgnDataHist", "", RooArgSet(mass), hMass, weight);

	return sgnDataHist;

}

TH1F * histoMtProducer(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight){

	TH1F* hMass = new TH1F("hMass","",50,70,120); 
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw("mass>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("mass>>hMass");
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	return hMass;

}

TH1F * histoMtProducer2(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight){

	TH1F* hMt = new TH1F("hMt","",30,0,120); 
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw("tag_Mt>>hMt","tag_puMCWeightRun2012");
		hMt->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("tag_Mt>>hMt");
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMt->Integral(0,121)<<endl;

	return hMt;

}

TH1F * histoVtxProducer(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight){

	TH1F* hMass = new TH1F("hMass","",40,0,40); 
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw("event_nPV>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("event_nPV>>hMass");
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,61)<<endl;

	return hMass;

}

void checkVisMass(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	double nBins_                    = 50,
	double xLow_                     = 70,
	double xHigh_                    = 120,
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
  //float readEventsTreeSgn = fullTreeSgn->GetEntries();
  //float effSgn = readEventsTreeSgn / readEventsSgn;

  // WJets
  TFile fWJets((path + "testTagAndProbe_WJets_TNP_TNP.root").c_str());
  fWJets.cd("counter");
  TH1F* totalEventsWJets = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWJets = totalEventsWJets->GetBinContent(1);
  TTree *fullTreeWJets = (TTree*)fWJets.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeWJets = fullTreeWJets->GetEntries();
  //float effWJets = readEventsTreeWJets / readEventsWJets;

  // Ztautau
  TFile fZtt((path + "testTagAndProbe_DYToTauTau_TNP_TNP.root").c_str());
  fZtt.cd("counter");
  TH1F* totalEventsZtt = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsZtt = totalEventsZtt->GetBinContent(1);
  TTree *fullTreeZtt = (TTree*)fZtt.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeZtt = fullTreeZtt->GetEntries();
  //float effZtt = readEventsTreeZtt / readEventsZtt;

  // TTJets
  TFile fTTJets((path + "testTagAndProbe_TTJets_TNP_TNP.root").c_str());
  fTTJets.cd("counter");
  TH1F* totalEventsTTJets = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsTTJets = totalEventsTTJets->GetBinContent(1);
  TTree *fullTreeTTJets = (TTree*)fTTJets.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeTTJets = fullTreeTTJets->GetEntries();
  //float effTTJets = readEventsTreeTTJets / readEventsTTJets;

  // WW
  TFile fWW((path + "testTagAndProbe_WW_TNP_TNP.root").c_str());
  fWW.cd("counter");
  TH1F* totalEventsWW = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWW = totalEventsWW->GetBinContent(1);
  TTree *fullTreeWW = (TTree*)fWW.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeWW = fullTreeWW->GetEntries();
  //float effWW = readEventsTreeWW / readEventsWW;

  // WZ
  TFile fWZ((path + "testTagAndProbe_WZ_TNP_TNP.root").c_str());
  fWZ.cd("counter");
  TH1F* totalEventsWZ = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsWZ = totalEventsWZ->GetBinContent(1);
  TTree *fullTreeWZ = (TTree*)fWZ.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeWZ = fullTreeWZ->GetEntries();
  //float effWZ = readEventsTreeWZ / readEventsWZ;

  // ZZ
  TFile fZZ((path + "testTagAndProbe_ZZ_TNP_TNP.root").c_str());
  fZZ.cd("counter");
  TH1F* totalEventsZZ = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsZZ = totalEventsZZ->GetBinContent(1);
  TTree *fullTreeZZ = (TTree*)fZZ.Get((tnp_+"/fitter_tree").c_str());
  //float readEventsTreeZZ = fullTreeZZ->GetEntries();
  //float effZZ = readEventsTreeZZ / readEventsZZ;

  // Data
  TFile fData((pathData + "testTagAndProbe_SingleMu_ABCD.root").c_str());
  TTree *fullTreeData = (TTree*)fData.Get((tnp_+"/fitter_tree").c_str());

  TFile *McP = new TFile("dummy3.root","RECREATE");

  // Create trees with cuts: passing

  std::cout<<"Processing Passing, OS, Low MT"<<std::endl;

  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );


  TH1F* hMet = new TH1F("hMet","",1,0,1500); 
  fullTreeDataCut->Draw("event_met_pfmet>>hMet");
  //float NumDataP = hMet->Integral();
  hMet->Reset();
  fullTreeSgnCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumSgnP = hMet->Integral();
  hMet->Reset();
  fullTreeWJetsCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWJetsP = hMet->Integral();
  hMet->Reset();
  fullTreeZttCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZttP = hMet->Integral();
  hMet->Reset();
  fullTreeTTJetsCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumTTJetsP = hMet->Integral();
  hMet->Reset();
  fullTreeWWCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWWP = hMet->Integral();
  hMet->Reset();
  fullTreeWZCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWZP = hMet->Integral();
  hMet->Reset();
  fullTreeZZCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZZP = hMet->Integral();
  hMet->Reset();

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

  fullTreeDataSSCut->Draw("event_met_pfmet>>hMet");
  //float NumDataSSP = hMet->Integral();
  hMet->Reset();
  fullTreeSgnSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumSgnSSP = hMet->Integral();
  hMet->Reset();
  fullTreeWJetsSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWJetsSSP = hMet->Integral();
  hMet->Reset();
  fullTreeZttSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZttSSP = hMet->Integral();
  hMet->Reset();
  fullTreeTTJetsSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumTTJetsSSP = hMet->Integral();
  hMet->Reset();
  fullTreeWWSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWWSSP = hMet->Integral();
  hMet->Reset();
  fullTreeWZSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWZSSP = hMet->Integral();
  hMet->Reset();
  fullTreeZZSSCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZZSSP = hMet->Integral();
  hMet->Reset();

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

  fullTreeDataHiMtCut->Draw("event_met_pfmet>>hMet");
  float NumDataHiMtP = hMet->Integral();
  cout<<"NumDataHiMtP "<<NumDataHiMtP<<endl;
  hMet->Reset();
  fullTreeSgnHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumSgnHiMtP = hMet->Integral();
  cout<<"NumSgnHiMtP "<<NumSgnHiMtP<<endl;
  hMet->Reset();
  fullTreeWJetsHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWJetsHiMtP = hMet->Integral();
  hMet->Reset();
  fullTreeZttHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZttHiMtP = hMet->Integral();
  hMet->Reset();
  fullTreeTTJetsHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumTTJetsHiMtP = hMet->Integral();
  hMet->Reset();
  fullTreeWWHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWWHiMtP = hMet->Integral();
  hMet->Reset();
  fullTreeWZHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumWZHiMtP = hMet->Integral();
  hMet->Reset();
  fullTreeZZHiMtCut->Draw("event_met_pfmet>>hMet","tag_puMCWeightRun2012");
  float NumZZHiMtP = hMet->Integral();
  hMet->Reset();

  delete hMet;

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
  //std::cout<<readEventsTreeWJets<<" "<<readEventsTreeZtt<<" "<<readEventsTreeTTJets<<" "<<readEventsTreeSgn<<" "<<readEventsTreeWW<<" "<<readEventsTreeWZ<<" "<<readEventsTreeZZ<<std::endl;
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
  TH1 * hDataDataHistP = DataDataHistP.createHistogram("mass",50);
  //float nPass = DataDataHistP.sum(false);

  RooDataSet dataDataSetSS("dataDataSetSS","dataset for signal-pass template SS", RooArgSet(mass), Import( *fullTreeDataSSCut ) );
  RooDataHist dataDataHistSS("dataDataHistSS", "", RooArgSet(mass), dataDataSetSS, 1.0);

  RooDataSet dataDataSetHiMt("dataDataSetHiMt","dataset for signal-pass template HiMt", RooArgSet(mass), Import( *fullTreeDataHiMtCut ) );
  RooDataHist dataDataHistHiMt("dataDataHistHiMt", "", RooArgSet(mass), dataDataSetHiMt, 1.0);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist sgnDataHist = dataSetProducer(fullTreeSgnCut, mass, NumSgnP, ScaleFactorSgn, 1.0);
  TH1 * hsgnDataHist = sgnDataHist.createHistogram("mass",50);
  RooHistPdf sgnTemplatePdf("sgnTemplatePdf", "", RooArgSet(mass), sgnDataHist,4);

  RooDataHist sgnDataHistSS = dataSetProducer(fullTreeSgnSSCut, mass, NumSgnSSP, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMt = dataSetProducer(fullTreeSgnHiMtCut, mass, NumSgnHiMtP, ScaleFactorSgn, -1.0);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHist = dataSetProducer(fullTreeWJetsCut, mass, NumWJetsP, ScaleFactorWJets, 1.0);
  RooHistPdf wjetsHistPdf("wjetsHistPdf", "", RooArgSet(mass), wjetsDataHist, 4);

  RooDataHist wjetsDataHistSS = dataSetProducer(fullTreeWJetsSSCut, mass, NumWJetsSSP, ScaleFactorWJets, -1.0);

  RooDataHist wjetsDataHistHiMt = dataSetProducer(fullTreeWJetsHiMtCut, mass, NumWJetsHiMtP, ScaleFactorWJets, -1.0);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zttDataHist = dataSetProducer(fullTreeZttCut, mass, NumZttP, ScaleFactorZtt, 1.0);
  TH1 * hzttDataHist = zttDataHist.createHistogram("mass",50);
  RooHistPdf zttPdf("zttPdf", "", RooArgSet(mass), zttDataHist, 4);  

  RooDataHist zttDataHistSS = dataSetProducer(fullTreeZttSSCut, mass, NumZttSSP, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMt = dataSetProducer(fullTreeZttHiMtCut, mass, NumZttHiMtP, ScaleFactorZtt, -1.0);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHist = dataSetProducer(fullTreeTTJetsCut, mass, NumTTJetsP, ScaleFactorTTJets, 1.0);
  TH1 * httjetsDataHist = ttjetsDataHist.createHistogram("mass",50);
  RooHistPdf ttjetsPdf("ttjetsPdf", "", RooArgSet(mass), ttjetsDataHist,4);

  RooDataHist ttjetsDataHistSS = dataSetProducer(fullTreeTTJetsSSCut, mass, NumTTJetsSSP, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMt = dataSetProducer(fullTreeTTJetsHiMtCut, mass, NumTTJetsHiMtP, ScaleFactorTTJets, -1.0);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHist = dataSetProducer(fullTreeWWCut, mass, NumWWP, ScaleFactorWW, 1.0);
  TH1 * hwwDataHist = wwDataHist.createHistogram("mass",50);
  RooHistPdf wwPdf("wwPdf", "", RooArgSet(mass), wwDataHist,4);

  mass.setBins( 50 );
  RooDataHist wwDataHistSS = dataSetProducer(fullTreeWWSSCut, mass, NumWWSSP, ScaleFactorWW, -1.0);

  mass.setBins( 50 );
  RooDataHist wwDataHistHiMt = dataSetProducer(fullTreeWWHiMtCut, mass, NumWWHiMtP, ScaleFactorWW, -1.0);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wzDataHist = dataSetProducer(fullTreeWZCut, mass, NumWZP, ScaleFactorWZ, 1.0);
  TH1 * hwzDataHist = wzDataHist.createHistogram("mass",50);
  RooHistPdf wzPdf("wzPdf", "", RooArgSet(mass), wzDataHist,4);

  mass.setBins( 50 );
  RooDataHist wzDataHistSS = dataSetProducer(fullTreeWZSSCut, mass, NumWZSSP, ScaleFactorWZ, -1.0);

  mass.setBins( 50 );
  RooDataHist wzDataHistHiMt = dataSetProducer(fullTreeWZHiMtCut, mass, NumWZHiMtP, ScaleFactorWZ, -1.0);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zzDataHist = dataSetProducer(fullTreeZZCut, mass, NumZZP, ScaleFactorZZ, 1.0);
  TH1 * hzzDataHist = zzDataHist.createHistogram("mass",50);
  RooHistPdf zzPdf("zzPdf", "", RooArgSet(mass), zzDataHist,4);

  mass.setBins( 50 );
  RooDataHist zzDataHistSS = dataSetProducer(fullTreeZZSSCut, mass, NumZZSSP, ScaleFactorZZ, -1.0);

  mass.setBins( 50 );
  RooDataHist zzDataHistHiMt = dataSetProducer(fullTreeZZHiMtCut, mass, NumZZHiMtP, ScaleFactorZZ, -1.0);

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
  float exWJetsP = NumWJetsP/NumWJetsHiMtP;
  cout<<"exWJetsP "<<exWJetsP<<endl;
  RooDataHist dataDataHistHiMtScaled("dataDataHistHiMtScaled", "", RooArgSet(mass), dataDataHistHiMt, exWJetsP);
  TH1 * hdataDataHistHiMtScaled = dataDataHistHiMtScaled.createHistogram("mass",50);
  //dataDataHistHiMtScaled.plotOn(massFrame1,MarkerColor(kViolet));
  //massFrame1->Draw();
  RooHistPdf WJetsDataDrivenHistPdfP("WJetsDataDrivenHistPdfP", "", RooArgSet(mass), dataDataHistHiMtScaled, 4);
  //WJetsDataDrivenHistPdfP.plotOn(massFrame1,MarkerColor(kRed));

  // QCD
  dataDataHistSS.add(sgnDataHistSS);
  dataDataHistSS.add(wjetsDataHistSS);
  dataDataHistSS.add(ttjetsDataHistSS);
  dataDataHistSS.add(zttDataHistSS);
  dataDataHistSS.add(wwDataHistSS);
  dataDataHistSS.add(wzDataHistSS);
  dataDataHistSS.add(zzDataHistSS);
  TH1 * hdataDataHistSS = dataDataHistSS.createHistogram("mass",50);
  RooHistPdf QCDHistPdfP("QCDHistPdfP", "", RooArgSet(mass), dataDataHistSS, 4);
  //QCDHistPdfP.plotOn(massFrame1,MarkerColor(kOrange));
  //massFrame1->Draw();

  delete fullTreeSgnCut;
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

  THStack * hs = new THStack("VisMass","VisMass");
  TLegend *leg2 = new TLegend(0.15,0.5,0.35,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);
  hdataDataHistHiMtScaled->SetLineColor(1);
  hdataDataHistHiMtScaled->SetFillColor(8);
  hs->Add(hdataDataHistHiMtScaled);
  hdataDataHistSS->SetLineColor(1);
  hdataDataHistSS->SetFillColor(6);
  hs->Add(hdataDataHistSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.5 fb^{-1}");
  hs->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "P");
  leg2->AddEntry(hsgnDataHist, "Signal", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  leg2->AddEntry(hdataDataHistHiMtScaled, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  leg2->AddEntry(hdataDataHistSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

}

void checkMass(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
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

  TFile *McP = new TFile("dummy4.root","RECREATE");

  //Low MT
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

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

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoMtProducer(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoMtProducer(fullTreeDataSSCut, 1, 1);
  TH1F * hDataDataHistPHiMt = histoMtProducer(fullTreeDataHiMtCut, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoMtProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoMtProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);
  TH1F * hsgnDataHistHiMt = histoMtProducer(fullTreeSgnHiMtCut, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoMtProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoMtProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);
  TH1F * hwjetsDataHistHiMt = histoMtProducer(fullTreeWJetsHiMtCut, ScaleFactorWJets, 1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoMtProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoMtProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);
  TH1F * hzttDataHistHiMt = histoMtProducer(fullTreeZttHiMtCut, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoMtProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoMtProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);
  TH1F * httjetsDataHistHiMt = histoMtProducer(fullTreeTTJetsHiMtCut, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoMtProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoMtProducer(fullTreeWWSSCut, ScaleFactorWW, -1);
  TH1F * hwwDataHistHiMt = histoMtProducer(fullTreeWWHiMtCut, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoMtProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoMtProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);
  TH1F * hwzDataHistHiMt = histoMtProducer(fullTreeWZHiMtCut, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoMtProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoMtProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);
  TH1F * hzzDataHistHiMt = histoMtProducer(fullTreeZZHiMtCut, ScaleFactorZZ, -1);

  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  hDataDataHistPHiMt->Add(hsgnDataHistHiMt);
  hDataDataHistPHiMt->Add(hzttDataHistHiMt);
  hDataDataHistPHiMt->Add(httjetsDataHistHiMt);
  hDataDataHistPHiMt->Add(hwwDataHistHiMt);
  hDataDataHistPHiMt->Add(hwzDataHistHiMt);
  hDataDataHistPHiMt->Add(hzzDataHistHiMt);
  float expW = hwjetsDataHist->Integral(0,51) / hwjetsDataHistHiMt->Integral(0,51);
  cout<<"expW "<<expW<<endl;
  hDataDataHistPHiMt->Scale(expW);

  THStack * hs = new THStack("VisMass","VisMass");
  TLegend *leg2 = new TLegend(0.6,0.5,0.85,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);

  //hwjetsDataHist->SetLineColor(1);
  //hwjetsDataHist->SetFillColor(8);
  //hs->Add(hwjetsDataHist);
  hDataDataHistPHiMt->SetLineColor(1);
  hDataDataHistPHiMt->SetFillColor(8);
  hs->Add(hDataDataHistPHiMt);

  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.5 fb^{-1}");
  hs->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  //hDataDataHistP->SetMarkerSize(1.5);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "p");
  leg2->AddEntry(hsgnDataHist, "Signal", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  //leg2->AddEntry(hwjetsDataHist, "WJets", "f");
  leg2->AddEntry(hDataDataHistPHiMt, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

  string fileName = "mass_"+tnp_+"_"+category_;
  string antiMu;
  string region;
  if(condition_ == ">=" && cutValue_ == 0.5) antiMu = "passing";
  else if(condition_ == "<" && cutValue_ == 0.5) antiMu = "failing";
  else if(condition_ == ">=" && cutValue_ == 0.0) antiMu = "all";
  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";
  c2->SaveAs(Form("%s_%s_%s.png",fileName.c_str(), antiMu.c_str(), region.c_str()));

}

void checkVtx(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.0,
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

  TFile *McP = new TFile("dummy4.root","RECREATE");

  // Create trees with cuts: passing

  //Low MT
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

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

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoVtxProducer(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoVtxProducer(fullTreeDataSSCut, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoVtxProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoVtxProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoVtxProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoVtxProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoVtxProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoVtxProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoVtxProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoVtxProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoVtxProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoVtxProducer(fullTreeWWSSCut, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoVtxProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoVtxProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoVtxProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoVtxProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);

  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  THStack * hs = new THStack("VisMass","VisMass");
  TLegend *leg2 = new TLegend(0.6,0.5,0.85,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);
  hwjetsDataHist->SetLineColor(1);
  hwjetsDataHist->SetFillColor(8);
  hs->Add(hwjetsDataHist);
  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.5 fb^{-1}");
  hs->GetXaxis()->SetTitle("# PV");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "P");
  leg2->AddEntry(hsgnDataHist, "Signal", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  leg2->AddEntry(hwjetsDataHist, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

  string fileName = "npv_"+tnp_+"_"+category_;
  string antiMu;
  string region;
  if(condition_ == ">=" && cutValue_ == 0.5) antiMu = "passing";
  else if(condition_ == "<" && cutValue_ == 0.5) antiMu = "failing";
  else if(condition_ == ">=" && cutValue_ == 0.0) antiMu = "all";
  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";
  c2->SaveAs(Form("%s_%s_%s.png",fileName.c_str(), antiMu.c_str(), region.c_str()));

}

void checkMassHiMt(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt > 60 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
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

  TFile *McP = new TFile("dummy4.root","RECREATE");

  // Create trees with cuts: passing

  //Low MT
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

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

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoMtProducer(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoMtProducer(fullTreeDataSSCut, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoMtProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoMtProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoMtProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoMtProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoMtProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoMtProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoMtProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoMtProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoMtProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoMtProducer(fullTreeWWSSCut, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoMtProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoMtProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoMtProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoMtProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);

  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  THStack * hs = new THStack("Mt","Mt");
  TLegend *leg2 = new TLegend(0.6,0.5,0.85,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);
  hwjetsDataHist->SetLineColor(1);
  hwjetsDataHist->SetFillColor(8);
  hs->Add(hwjetsDataHist);
  //hDataDataHistPSS->SetLineColor(1);
  //hDataDataHistPSS->SetFillColor(6);
  //hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.5 fb^{-1}");
  hs->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "P");
  leg2->AddEntry(hsgnDataHist, "Signal", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  leg2->AddEntry(hwjetsDataHist, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  //leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

  string fileName = "massHiMt_"+tnp_+"_"+category_;
  string antiMu;
  string region;
  if(condition_ == ">=" && cutValue_ == 0.5) antiMu = "passing";
  else if(condition_ == "<" && cutValue_ == 0.5) antiMu = "failing";
  else if(condition_ == ">=" && cutValue_ == 0.0) antiMu = "all";
  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";
  c2->SaveAs(Form("%s_%s_%s.png",fileName.c_str(), antiMu.c_str(), region.c_str()));

}

void checkMt(
	const string tnp_                = "muToTau",
	const string category_           = "passingIsoLooseMuonVetoLoose",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 400 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 400 && tag_IsoMu24_eta2p1 > 0.5 && tag_triggerBitSingleMu > 0.5",
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

  TFile *McP = new TFile("dummy4.root","RECREATE");

  //Low MT
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

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

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoMtProducer2(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoMtProducer2(fullTreeDataSSCut, 1, 1);
  TH1F * hDataDataHistPHiMt = histoMtProducer2(fullTreeDataHiMtCut, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoMtProducer2(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoMtProducer2(fullTreeSgnSSCut, ScaleFactorSgn, -1);
  TH1F * hsgnDataHistHiMt = histoMtProducer2(fullTreeSgnHiMtCut, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoMtProducer2(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoMtProducer2(fullTreeWJetsSSCut, ScaleFactorWJets, -1);
  TH1F * hwjetsDataHistHiMt = histoMtProducer2(fullTreeWJetsHiMtCut, ScaleFactorWJets, 1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoMtProducer2(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoMtProducer2(fullTreeZttSSCut, ScaleFactorZtt, -1);
  TH1F * hzttDataHistHiMt = histoMtProducer2(fullTreeZttHiMtCut, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoMtProducer2(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoMtProducer2(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);
  TH1F * httjetsDataHistHiMt = histoMtProducer2(fullTreeTTJetsHiMtCut, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoMtProducer2(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoMtProducer2(fullTreeWWSSCut, ScaleFactorWW, -1);
  TH1F * hwwDataHistHiMt = histoMtProducer2(fullTreeWWHiMtCut, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoMtProducer2(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoMtProducer2(fullTreeWZSSCut, ScaleFactorWZ, -1);
  TH1F * hwzDataHistHiMt = histoMtProducer2(fullTreeWZHiMtCut, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoMtProducer2(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoMtProducer2(fullTreeZZSSCut, ScaleFactorZZ, -1);
  TH1F * hzzDataHistHiMt = histoMtProducer2(fullTreeZZHiMtCut, ScaleFactorZZ, -1);

  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  hDataDataHistPHiMt->Add(hsgnDataHistHiMt);
  hDataDataHistPHiMt->Add(hzttDataHistHiMt);
  hDataDataHistPHiMt->Add(httjetsDataHistHiMt);
  hDataDataHistPHiMt->Add(hwwDataHistHiMt);
  hDataDataHistPHiMt->Add(hwzDataHistHiMt);
  hDataDataHistPHiMt->Add(hzzDataHistHiMt);
  float expW = hwjetsDataHist->Integral(0,121) / hwjetsDataHistHiMt->Integral(0,121);
  cout<<"expW "<<expW<<endl;
  hDataDataHistPHiMt->Scale(expW);

  THStack * hs = new THStack("Mt","Mt");
  TLegend *leg2 = new TLegend(0.6,0.5,0.85,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);

  hwjetsDataHist->SetLineColor(1);
  hwjetsDataHist->SetFillColor(8);
  hs->Add(hwjetsDataHist);
  //hDataDataHistPHiMt->SetLineColor(1);
  //hDataDataHistPHiMt->SetFillColor(8);
  //hs->Add(hDataDataHistPHiMt);

  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.5 fb^{-1}");
  hs->GetXaxis()->SetTitle("M_{T}(#mu,MET) [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  //hDataDataHistP->SetMarkerSize(1.5);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "p");
  leg2->AddEntry(hsgnDataHist, "Signal", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  leg2->AddEntry(hwjetsDataHist, "WJets", "f");
  //leg2->AddEntry(hDataDataHistPHiMt, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

  string fileName = "mt_"+tnp_+"_"+category_;
  string antiMu;
  string region;
  if(condition_ == ">=" && cutValue_ == 0.5) antiMu = "passing";
  else if(condition_ == "<" && cutValue_ == 0.5) antiMu = "failing";
  else if(condition_ == ">=" && cutValue_ == 0.0) antiMu = "all";
  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";
  c2->SaveAs(Form("%s_%s_%s.png",fileName.c_str(), antiMu.c_str(), region.c_str()));

}

void calcutateMass(){

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7","<",0.5);

	checkMass("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2",">=",0.0);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7",">=",0.5);
	checkMass("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7","<",0.5);

}

void calcutateMassHiMt(){

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta<1.2","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.2 && abseta<1.7","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoLoose","abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta<1.2","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.2 && abseta<1.7","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoMedium","abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta<1.2","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta>1.2 && abseta<1.7","<",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau","passingIsoLooseMuonVetoTight","abseta>1.7","<",0.5);

}

