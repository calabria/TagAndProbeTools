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
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
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
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
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
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
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
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,61)<<endl;

	return hMass;

}

TH1F * histoMuPtProducer(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight){

	TH1F* hMass = new TH1F("hMass","",50,0,100);
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw("tag_pt>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("tag_pt>>hMass");
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	return hMass;

}

TH1F * histoTauPtProducer(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight){

	TH1F* hMass = new TH1F("hMass","",50,0,100);
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw("pt>>hMass","tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw("pt>>hMass");
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	return hMass;

}

/*void openCreateTrees(std::string path, std::string dir, std::string nameIn, std::string nameOut, std::string cut, std::string bin, std::string additionalCut){

  // signal
  TFile fsgn((path + nameIn).c_str());
  fsgn.cd("counter");
  TH1F* totalEventsSgn = (TH1F*)gDirectory->Get("N_eventi_Tot");
  float readEventsSgn = totalEventsSgn->GetBinContent(1);
  TTree *fullTreeSgn = (TTree*)fsgn.Get((dir+"/fitter_tree").c_str());
  //float readEventsTreeSgn = fullTreeSgn->GetEntries();
  //float effSgn = readEventsTreeSgn / readEventsSgn;

  TFile *McP = new TFile(nameOut.c_str(),"RECREATE");
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",cut.c_str(),bin.c_str(),additionalCut.c_str()) );

}

void createAllNewTrees(std::string path, std::string dir, std::string nameIn, std::string nameOut, std::string cut, std::string bin){

	std::string additionalCut = "tag_pt > 25 && tag_muPFIsolation < 0.1 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	std::string additionalCutSS = "tag_pt > 25 && tag_muPFIsolation < 0.1 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	std::string additionalCutHiMt = "tag_pt > 25 && tag_muPFIsolation < 0.1 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"

	openCreateTrees(path, dir, nameIn, nameOut, cut, bin, additionalCut);
	openCreateTrees(path, dir, nameIn, nameOut, cut, bin, additionalCutSS);
	openCreateTrees(path, dir, nameIn, nameOut, cut, bin, additionalCutHiMt);

}*/

void checkVisMass(
	const string tnp_                = "muToTau",
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	double nBins_                    = 50,
	double xLow_                     = 70,
	double xHigh_                    = 120,
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"
	){
  
  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  string path = "./Input2014/";
  string pathData = "./Input2014/";

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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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
  TH1F * hdataDataHistHiMtScaled = (TH1F*)dataDataHistHiMtScaled.createHistogram("mass",50);
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
  TH1F * hdataDataHistSS = (TH1F*)dataDataHistSS.createHistogram("mass",50);
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
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
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
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"
	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,700,700);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  TFile *McP = new TFile("dummy4.root","RECREATE");

  //Low MT
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutTempNo = fullTreeSgn->CopyTree( Form("(!mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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
  TH1F * hsgnDataHistTemp = histoMtProducer(fullTreeSgnCutTemp, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistTempNo = histoMtProducer(fullTreeSgnCutTempNo, ScaleFactorSgn, 1);
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
  THStack * hs2 = new THStack("VisMass","VisMass");  
  TLegend *leg2 = new TLegend(0.6,0.5,0.85,0.8);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  hzttDataHist->SetLineColor(1);
  hzttDataHist->SetFillColor(95);
  hs->Add(hzttDataHist);
  hs2->Add(hzttDataHist);
  httjetsDataHist->SetLineColor(1);
  httjetsDataHist->SetFillColor(4);
  hs->Add(httjetsDataHist);
  hs2->Add(httjetsDataHist);
  hwwDataHist->SetLineColor(1);
  hwwDataHist->SetFillColor(5);
  hs->Add(hwwDataHist);
  hs2->Add(hwwDataHist);
  hwzDataHist->SetLineColor(1);
  hwzDataHist->SetFillColor(50);
  hs->Add(hwzDataHist);
  hs2->Add(hwzDataHist);
  hzzDataHist->SetLineColor(1);
  hzzDataHist->SetFillColor(7);
  hs->Add(hzzDataHist);
  hs2->Add(hzzDataHist);

  //hwjetsDataHist->SetLineColor(1);
  //hwjetsDataHist->SetFillColor(8);
  //hs->Add(hwjetsDataHist);
  hDataDataHistPHiMt->SetLineColor(1);
  hDataDataHistPHiMt->SetFillColor(8);
  hs->Add(hDataDataHistPHiMt);
  hs2->Add(hDataDataHistPHiMt);

  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hs2->Add(hDataDataHistPSS);
  hsgnDataHistTemp->SetLineColor(1);
  hsgnDataHistTemp->SetFillColor(2);
  hs->Add(hsgnDataHistTemp);
  hs2->Add(hsgnDataHistTemp);
  hsgnDataHistTempNo->SetLineColor(1);
  hsgnDataHistTempNo->SetFillColor(40);
  hs->Add(hsgnDataHistTempNo);
  hs2->Add(hsgnDataHistTempNo);

  //gPad->SetLogy();
  //hs->SetMinimum(1);

  TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,10);
  TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,10);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  pad1->SetTickx();
  pad1->SetTicky();

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
  hs->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  //hDataDataHistP->SetMarkerSize(1.5);
  hDataDataHistP->Draw("SAME,ep");

  leg2->AddEntry(hsgnDataHist, "Data", "p");
  leg2->AddEntry(hsgnDataHistTemp, "Signal (#mu#rightarrow#tau)", "f");
  leg2->AddEntry(hsgnDataHistTempNo, "Signal (jet#rightarrow#tau)", "f");
  leg2->AddEntry(hzttDataHist, "Z#rightarrow#tau#tau", "f");
  //leg2->AddEntry(hwjetsDataHist, "WJets", "f");
  leg2->AddEntry(hDataDataHistPHiMt, "WJets", "f");
  leg2->AddEntry(httjetsDataHist, "TTJets", "f");
  leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
  leg2->AddEntry(hwwDataHist, "WW", "f");
  leg2->AddEntry(hwzDataHist, "WZ", "f");
  leg2->AddEntry(hzzDataHist, "ZZ", "f");
  leg2->Draw();

  pad2->cd();

  pad2->SetTickx();
  pad2->SetTicky();

  TH1F * stack = (TH1F*)hs2->GetStack()->Last()->Clone();
  TH1F * hDataDataHistPClone = (TH1F*)hDataDataHistP->Clone();
  hDataDataHistPClone->Divide(stack);
  hDataDataHistPClone->SetStats(kFALSE);
  hDataDataHistPClone->SetMarkerStyle(20);
  hDataDataHistPClone->SetMarkerColor(1);
  hDataDataHistPClone->SetMinimum(0);
  hDataDataHistPClone->SetMaximum(+2);
  hDataDataHistPClone->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  hDataDataHistPClone->GetYaxis()->SetTitle("obs/exp");
  hDataDataHistPClone->Draw("ep");

  string fileName = "mass_"+tnp_+"_"+category_;
  string antiMu;
  string region;
  if(condition_ == ">=" && cutValue_ == 0.5) antiMu = "passing";
  else if(condition_ == "<" && cutValue_ == 0.5) antiMu = "failing";
  else if(condition_ == ">=" && cutValue_ == 0.0) antiMu = "all";
  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";
  else if(bin_ == "abseta<0.4") region = "barrel04";
  else if(bin_ == "abseta>0.4 && abseta<0.8") region = "barrel08";
  else if(bin_ == "abseta>0.8 && abseta<1.2") region = "barrel12";
  else if(bin_ == "abseta<0.2") region = "barrel02b";
  else if(bin_ == "abseta>0.2 && abseta<0.4") region = "barrel04b";
  else if(bin_ == "abseta>0.4 && abseta<0.6") region = "barrel06b";
  else if(bin_ == "abseta>0.6 && abseta<0.8") region = "barrel08b";
  else if(bin_ == "abseta>0.8 && abseta<1.0") region = "barrel10b";
  else if(bin_ == "abseta>1.0 && abseta<1.2") region = "barrel12b";
  c2->SaveAs(Form("%s_%s_%s.png",fileName.c_str(), antiMu.c_str(), region.c_str()));

}

void checkVtx(
	const string tnp_                = "muToTau",
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.0,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"
	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
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
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,

	const std::string additionalCut_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)",
	const std::string additionalCutSS_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)",
	const std::string additionalCutHiMt_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)"

	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  //Low MT No Discriminator
  TTree* fullTreeSgnCut2 = fullTreeSgn->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsCut2 = fullTreeWJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttCut2 = fullTreeZtt->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsCut2 = fullTreeTTJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWCut2 = fullTreeWW->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZCut2 = fullTreeWZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZCut2 = fullTreeZZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeDataCut2 = fullTreeData->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );

  //SS No Cut
  TTree* fullTreeDataSSCut2 = fullTreeData->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut2 = fullTreeSgn->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut2 = fullTreeWJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut2 = fullTreeZtt->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut2 = fullTreeTTJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut2 = fullTreeWW->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut2 = fullTreeWZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut2 = fullTreeZZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );

  //Normalization passing

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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

  TH1F * hDataDataHistP2 = histoMtProducer(fullTreeDataCut2, 1, 1);
  TH1F * hDataDataHistPSS2 = histoMtProducer(fullTreeDataSSCut2, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoMtProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoMtProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);

  TH1F * hsgnDataHist2 = histoMtProducer(fullTreeSgnCut2, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS2 = histoMtProducer(fullTreeSgnSSCut2, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoMtProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoMtProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);

  TH1F * hwjetsDataHist2 = histoMtProducer(fullTreeWJetsCut2, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS2 = histoMtProducer(fullTreeWJetsSSCut2, ScaleFactorWJets, -1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoMtProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoMtProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);

  TH1F * hzttDataHist2 = histoMtProducer(fullTreeZttCut2, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS2 = histoMtProducer(fullTreeZttSSCut2, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoMtProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoMtProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);

  TH1F * httjetsDataHist2 = histoMtProducer(fullTreeTTJetsCut2, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS2 = histoMtProducer(fullTreeTTJetsSSCut2, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoMtProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoMtProducer(fullTreeWWSSCut, ScaleFactorWW, -1);

  TH1F * hwwDataHist2 = histoMtProducer(fullTreeWWCut2, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS2 = histoMtProducer(fullTreeWWSSCut2, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoMtProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoMtProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);

  TH1F * hwzDataHist2 = histoMtProducer(fullTreeWZCut2, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS2 = histoMtProducer(fullTreeWZSSCut2, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoMtProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoMtProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);

  TH1F * hzzDataHist2 = histoMtProducer(fullTreeZZCut2, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS2 = histoMtProducer(fullTreeZZSSCut2, ScaleFactorZZ, -1);


  hDataDataHistP->Sumw2();
  hDataDataHistPSS->Sumw2();

  hDataDataHistP2->Sumw2(); 
  hDataDataHistPSS2->Sumw2(); 

  hsgnDataHist->Sumw2(); 
  hsgnDataHistSS->Sumw2(); 

  hsgnDataHist2->Sumw2(); 
  hsgnDataHistSS2->Sumw2(); 

  hwjetsDataHist->Sumw2(); 
  hwjetsDataHistSS->Sumw2(); 

  hwjetsDataHist2->Sumw2(); 
  hwjetsDataHistSS2->Sumw2(); 

  hzttDataHist->Sumw2(); 
  hzttDataHistSS->Sumw2(); 

  hzttDataHist2->Sumw2(); 
  hzttDataHistSS2->Sumw2();

  httjetsDataHist->Sumw2();
  httjetsDataHistSS->Sumw2();

  httjetsDataHist2->Sumw2(); 
  httjetsDataHistSS2->Sumw2();

  hwwDataHist->Sumw2();
  hwwDataHistSS->Sumw2(); 

  hwwDataHist2->Sumw2(); 
  hwwDataHistSS2->Sumw2(); 

  hwzDataHist->Sumw2();
  hwzDataHistSS->Sumw2(); 

  hwzDataHist2->Sumw2(); 
  hwzDataHistSS2->Sumw2();

  hzzDataHist->Sumw2(); 
  hzzDataHistSS->Sumw2(); 

  hzzDataHist2->Sumw2(); 
  hzzDataHistSS2->Sumw2();


  hDataDataHistPSS2->Add(hsgnDataHistSS2);
  hDataDataHistPSS2->Add(hwjetsDataHistSS2);
  hDataDataHistPSS2->Add(hzttDataHistSS2);
  hDataDataHistPSS2->Add(httjetsDataHistSS2);
  hDataDataHistPSS2->Add(hwwDataHistSS2);
  hDataDataHistPSS2->Add(hwzDataHistSS2);
  hDataDataHistPSS2->Add(hzzDataHistSS2);

  hsgnDataHist2->Add(hwjetsDataHist2);
  hsgnDataHist2->Add(hzttDataHist2);
  hsgnDataHist2->Add(httjetsDataHist2);
  hsgnDataHist2->Add(hwwDataHist2);
  hsgnDataHist2->Add(hwzDataHist2);
  hsgnDataHist2->Add(hzzDataHist2);

  hDataDataHistP2->Add(hDataDataHistPSS2, -1);
  hDataDataHistP2->Divide(hsgnDataHist2);

  hsgnDataHist->Multiply(hDataDataHistP2); //rescale of the signal

  for(int i = 0; i < hDataDataHistP2->GetSize(); i++){

 	std::cout<<"Bin "<<i<<" "<<hDataDataHistP2->GetBinContent(i)<<std::endl;
 	std::cout<<"Error "<<i<<" "<<hDataDataHistP2->GetBinError(i)<<std::endl;

  }

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
  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
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
  leg2->AddEntry(hDataDataHistPSS, "QCD", "f");
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
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"
	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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
  TLegend *leg2 = new TLegend(0.7,0.6,0.95,0.9);
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
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
  hs->GetXaxis()->SetTitle("M_{T}(#mu,MET) [GeV/c^{2}]");
  hs->GetYaxis()->SetTitle("Entries");
  hDataDataHistP->SetMarkerStyle(20);
  hDataDataHistP->SetMarkerColor(1);
  hDataDataHistP->SetMarkerSize(1.5);
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

void checkMuPtHiMt(
	const string tnp_                = "muToTau",
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"

	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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

  TH1F * hDataDataHistP = histoMuPtProducer(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoMuPtProducer(fullTreeDataSSCut, 1, 1);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoMuPtProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoMuPtProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoMuPtProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoMuPtProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoMuPtProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoMuPtProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoMuPtProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoMuPtProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoMuPtProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoMuPtProducer(fullTreeWWSSCut, ScaleFactorWW, -1);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoMuPtProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoMuPtProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoMuPtProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoMuPtProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);

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
  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
  hs->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
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

  string fileName = "muPt_"+tnp_+"_"+category_;
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

void checkTauPtHiMt(
	const string tnp_                = "muToTau",
	const string category_           = "",
	const string bin_                = "abseta<1.2",
	const string condition_          = ">=",
	double cutValue_                 = 0.5,
	const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5",
	const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5"

	){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,650,600);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

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

  float Lumi_ = 19988;

  //Cross sections
  float SigmaWJets = 36267.2;
  float SigmaZtt = 1966.7;
  float SigmaTTJets = 245.8;
  float SigmaSgn = 1966.7;
  float SigmaWW = 56.0;
  float SigmaWZ = 33.6;
  float SigmaZZ = 17.0;

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

  TH1F * hDataDataHistP = histoTauPtProducer(fullTreeDataCut, 1, 1);
  TH1F * hDataDataHistPSS = histoTauPtProducer(fullTreeDataSSCut, 1, 1);

  hDataDataHistP->Sumw2();
  hDataDataHistPSS->Sumw2();

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoTauPtProducer(fullTreeSgnCut, ScaleFactorSgn, 1);
  TH1F * hsgnDataHistSS = histoTauPtProducer(fullTreeSgnSSCut, ScaleFactorSgn, -1);

  hsgnDataHist->Sumw2();
  hsgnDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoTauPtProducer(fullTreeWJetsCut, ScaleFactorWJets, 1);
  TH1F * hwjetsDataHistSS = histoTauPtProducer(fullTreeWJetsSSCut, ScaleFactorWJets, -1);

  hwjetsDataHist->Sumw2();
  hwjetsDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoTauPtProducer(fullTreeZttCut, ScaleFactorZtt, 1);
  TH1F * hzttDataHistSS = histoTauPtProducer(fullTreeZttSSCut, ScaleFactorZtt, -1);

  hzttDataHist->Sumw2();
  hzttDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoTauPtProducer(fullTreeTTJetsCut, ScaleFactorTTJets, 1);
  TH1F * httjetsDataHistSS = histoTauPtProducer(fullTreeTTJetsSSCut, ScaleFactorTTJets, -1);

  httjetsDataHist->Sumw2();
  httjetsDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoTauPtProducer(fullTreeWWCut, ScaleFactorWW, 1);
  TH1F * hwwDataHistSS = histoTauPtProducer(fullTreeWWSSCut, ScaleFactorWW, -1);

  hwwDataHist->Sumw2();
  hwwDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoTauPtProducer(fullTreeWZCut, ScaleFactorWZ, 1);
  TH1F * hwzDataHistSS = histoTauPtProducer(fullTreeWZSSCut, ScaleFactorWZ, -1);

  hwzDataHist->Sumw2();
  hwzDataHistSS->Sumw2();

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoTauPtProducer(fullTreeZZCut, ScaleFactorZZ, 1);
  TH1F * hzzDataHistSS = histoTauPtProducer(fullTreeZZSSCut, ScaleFactorZZ, -1);

  hzzDataHist->Sumw2();
  hzzDataHistSS->Sumw2();

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
  hDataDataHistPSS->SetLineColor(1);
  hDataDataHistPSS->SetFillColor(6);
  hs->Add(hDataDataHistPSS);
  hsgnDataHist->SetLineColor(1);
  hsgnDataHist->SetFillColor(2);
  hs->Add(hsgnDataHist);

  hs->Draw("HIST");
  hs->SetTitle("CMS Preliminary 2012 #int L = 19.9 fb^{-1}");
  hs->GetXaxis()->SetTitle("p_{T}^{#tau} [GeV/c]");
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

  string fileName = "tauPt_"+tnp_+"_"+category_;
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

	std::string looseV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose > 0.5";
	std::string mediumV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium > 0.5";
	std::string tightV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight > 0.5";

	std::string looseV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose2 > 0.5";
	std::string mediumV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium2 > 0.5";
	std::string tightV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight2 > 0.5";

	std::string looseV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose3 > 0.5";
	std::string tightV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight3 > 0.5";

	std::string looseMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLooseMVA > 0.5";
	std::string mediumMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMediumMVA > 0.5";
	std::string tightMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTightMVA > 0.5";

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkMass("muToTau",looseV1,"abseta<1.2",">=",0.5);
	checkMass("muToTau",looseV1,"abseta<1.2","<",0.5);
	checkMass("muToTau",looseV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau",looseV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau",looseV1,"abseta>1.7",">=",0.5);
	checkMass("muToTau",looseV1,"abseta>1.7","<",0.5);

	//checkMass("muToTau",looseV1,"abseta<1.2",">=",0.0);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMass("muToTau",mediumV1,"abseta<1.2",">=",0.5);
	checkMass("muToTau",mediumV1,"abseta<1.2","<",0.5);
	checkMass("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau",mediumV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau",mediumV1,"abseta>1.7",">=",0.5);
	checkMass("muToTau",mediumV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMass("muToTau",tightV1,"abseta<1.2",">=",0.5);
	checkMass("muToTau",tightV1,"abseta<1.2","<",0.5);
	checkMass("muToTau",tightV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMass("muToTau",tightV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMass("muToTau",tightV1,"abseta>1.7",">=",0.5);
	checkMass("muToTau",tightV1,"abseta>1.7","<",0.5);

}

void calcutateMassHiMt(){

	std::string looseV1 = "HpsAntiMuLoose";
	std::string mediumV1 = "HpsAntiMuMedium";
	std::string tightV1 = "HpsAntiMuTight";

	std::string looseV2 = "HpsAntiMuLoose2";
	std::string mediumV2 = "HpsAntiMuMedium2";
	std::string tightV2 = "HpsAntiMuTight2";

	std::string looseV3 = "HpsAntiMuLoose3";
	std::string tightV3 = "HpsAntiMuTight3";

	std::string looseMVA = "HpsAntiMuLooseMVA";
	std::string mediumMVA = "HpsAntiMuMediumMVA";
	std::string tightMVA = "HpsAntiMuTightMVA";

	cout<<"passingIsoLooseMuonVetoLoose Barrel"<<endl;
	//checkMassHiMt("muToTau",looseV1,"abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau",looseV1,"abseta<1.2","<",0.5);
	cout<<"passingIsoLooseMuonVetoLoose Overlap"<<endl;
	//checkMassHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7","<",0.5);
	cout<<"passingIsoLooseMuonVetoLoose Endcap"<<endl;
	//checkMassHiMt("muToTau",looseV1,"abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau",looseV1,"abseta>1.7","<",0.5);

	/*cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMassHiMt("muToTau",mediumV1,"abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau",mediumV1,"abseta<1.2","<",0.5);
	checkMassHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMassHiMt("muToTau",mediumV1,"abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau",mediumV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMassHiMt("muToTau",tightV1,"abseta<1.2",">=",0.5);
	checkMassHiMt("muToTau",tightV1,"abseta<1.2","<",0.5);
	checkMassHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMassHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMassHiMt("muToTau",tightV1,"abseta>1.7",">=",0.5);
	checkMassHiMt("muToTau",tightV1,"abseta>1.7","<",0.5);*/

}

void checkMT(){

	std::string looseV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose > 0.5";
	std::string mediumV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium > 0.5";
	std::string tightV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight > 0.5";

	std::string looseV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose2 > 0.5";
	std::string mediumV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium2 > 0.5";
	std::string tightV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight2 > 0.5";

	std::string looseV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose3 > 0.5";
	std::string tightV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight3 > 0.5";

	std::string looseMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLooseMVA > 0.5";
	std::string mediumMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMediumMVA > 0.5";
	std::string tightMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTightMVA > 0.5";

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkMt("muToTau",looseV1,"abseta<1.2",">=",0.5);
	checkMt("muToTau",looseV1,"abseta<1.2","<",0.5);
	checkMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMt("muToTau",looseV1,"abseta>1.7",">=",0.5);
	checkMt("muToTau",looseV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMt("muToTau",mediumV1,"abseta<1.2",">=",0.5);
	checkMt("muToTau",mediumV1,"abseta<1.2","<",0.5);
	checkMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMt("muToTau",mediumV1,"abseta>1.7",">=",0.5);
	checkMt("muToTau",mediumV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMt("muToTau",tightV1,"abseta<1.2",">=",0.5);
	checkMt("muToTau",tightV1,"abseta<1.2","<",0.5);
	checkMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMt("muToTau",tightV1,"abseta>1.7",">=",0.5);
	checkMt("muToTau",tightV1,"abseta>1.7","<",0.5);

}

void checkMU(){

	std::string looseV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose > 0.5";
	std::string mediumV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium > 0.5";
	std::string tightV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight > 0.5";

	std::string looseV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose2 > 0.5";
	std::string mediumV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium2 > 0.5";
	std::string tightV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight2 > 0.5";

	std::string looseV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose3 > 0.5";
	std::string tightV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight3 > 0.5";

	std::string looseMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLooseMVA > 0.5";
	std::string mediumMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMediumMVA > 0.5";
	std::string tightMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTightMVA > 0.5";

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkMuPtHiMt("muToTau",looseV1,"abseta<1.2",">=",0.5);
	checkMuPtHiMt("muToTau",looseV1,"abseta<1.2","<",0.5);
	checkMuPtHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMuPtHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMuPtHiMt("muToTau",looseV1,"abseta>1.7",">=",0.5);
	checkMuPtHiMt("muToTau",looseV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkMuPtHiMt("muToTau",mediumV1,"abseta<1.2",">=",0.5);
	checkMuPtHiMt("muToTau",mediumV1,"abseta<1.2","<",0.5);
	checkMuPtHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMuPtHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMuPtHiMt("muToTau",mediumV1,"abseta>1.7",">=",0.5);
	checkMuPtHiMt("muToTau",mediumV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkMuPtHiMt("muToTau",tightV1,"abseta<1.2",">=",0.5);
	checkMuPtHiMt("muToTau",tightV1,"abseta<1.2","<",0.5);
	checkMuPtHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkMuPtHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkMuPtHiMt("muToTau",tightV1,"abseta>1.7",">=",0.5);
	checkMuPtHiMt("muToTau",tightV1,"abseta>1.7","<",0.5);

}

void checkTAU(){

	std::string looseV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose > 0.5";
	std::string mediumV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium > 0.5";
	std::string tightV1 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight > 0.5";

	std::string looseV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose2 > 0.5";
	std::string mediumV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMedium2 > 0.5";
	std::string tightV2 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight2 > 0.5";

	std::string looseV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLoose3 > 0.5";
	std::string tightV3 = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTight3 > 0.5";

	std::string looseMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuLooseMVA > 0.5";
	std::string mediumMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuMediumMVA > 0.5";
	std::string tightMVA = "pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && HpsAntiMuTightMVA > 0.5";

	cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	checkTauPtHiMt("muToTau",looseV1,"abseta<1.2",">=",0.5);
	checkTauPtHiMt("muToTau",looseV1,"abseta<1.2","<",0.5);
	checkTauPtHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkTauPtHiMt("muToTau",looseV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkTauPtHiMt("muToTau",looseV1,"abseta>1.7",">=",0.5);
	checkTauPtHiMt("muToTau",looseV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	checkTauPtHiMt("muToTau",mediumV1,"abseta<1.2",">=",0.5);
	checkTauPtHiMt("muToTau",mediumV1,"abseta<1.2","<",0.5);
	checkTauPtHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkTauPtHiMt("muToTau",mediumV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkTauPtHiMt("muToTau",mediumV1,"abseta>1.7",">=",0.5);
	checkTauPtHiMt("muToTau",mediumV1,"abseta>1.7","<",0.5);

	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	checkTauPtHiMt("muToTau",tightV1,"abseta<1.2",">=",0.5);
	checkTauPtHiMt("muToTau",tightV1,"abseta<1.2","<",0.5);
	checkTauPtHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7",">=",0.5);
	checkTauPtHiMt("muToTau",tightV1,"abseta>1.2 && abseta<1.7","<",0.5);
	checkTauPtHiMt("muToTau",tightV1,"abseta>1.7",">=",0.5);
	checkTauPtHiMt("muToTau",tightV1,"abseta>1.7","<",0.5);

}

