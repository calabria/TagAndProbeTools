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

RooDataHist dataSetProducer(TTree * fullTreeSgnCut, RooRealVar mass, float NumSgnP, float ScaleFactorSgn, float weight){

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
  double SGNtrue = hS->Integral(0,151);
  fullTreeSgn->Draw("mass>>hSP",Form("tag_puMCWeightRun2012*(%s && %s>=%f && mass>%f && mass<%f && mcTrue && pair_charge==0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_triggerBitSingleMu > 0.5 )",bin_.c_str(),category_.c_str(),cutValue_,xLow_,xHigh_));
  double SGNtruePass = hSP->Integral(0,151);

  double McTruthEff    = SGNtruePass/SGNtrue;
  double BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  delete hS; delete hSP;

  // Create trees with cuts: passing

  std::cout<<"Processing Passing, OS, Low MT"<<std::endl;

  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  //TTree* fullTreeSgnCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
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
  //TTree* fullTreeSgnSSCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
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
  //TTree* fullTreeSgnHiMtCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s%s%f && %s && %s)",category_.c_str(),condition_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
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

  //float readEventsWJets = 76102995 -4*115419;
  float ScaleFactorWJets = Lumi_/(readEventsWJets/SigmaWJets);

  //float readEventsZtt = 19937479;
  float ScaleFactorZtt = Lumi_/(readEventsZtt/SigmaZtt);

  //float readEventsTTJets = 6923750 -2*13847;
  float ScaleFactorTTJets = Lumi_/(readEventsTTJets/SigmaTTJets);

  //float readEventsSgn = 29743564 - 6*97638;
  float ScaleFactorSgn = Lumi_/(readEventsSgn/SigmaSgn);

  //float readEventsWW = 10000431 - 20000;
  float ScaleFactorWW = Lumi_/(readEventsWW/SigmaWW);

  //float readEventsWZ = 10000283;
  float ScaleFactorWZ = Lumi_/(readEventsWZ/SigmaWZ);

  //float readEventsZZ = 9799908 -2*19599;
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
  RooHistPdf zttPdf("zttPdf", "", RooArgSet(mass), zttDataHist, 4);  

  RooDataHist zttDataHistSS = dataSetProducer(fullTreeZttSSCut, mass, NumZttSSP, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMt = dataSetProducer(fullTreeZttHiMtCut, mass, NumZttHiMtP, ScaleFactorZtt, -1.0);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHist = dataSetProducer(fullTreeTTJetsCut, mass, NumTTJetsP, ScaleFactorTTJets, 1.0);
  RooHistPdf ttjetsPdf("ttjetsPdf", "", RooArgSet(mass), ttjetsDataHist,4);

  RooDataHist ttjetsDataHistSS = dataSetProducer(fullTreeTTJetsSSCut, mass, NumTTJetsSSP, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMt = dataSetProducer(fullTreeTTJetsHiMtCut, mass, NumTTJetsHiMtP, ScaleFactorTTJets, -1.0);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHist = dataSetProducer(fullTreeWWCut, mass, NumWWP, ScaleFactorWW, 1.0);
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

  delete fullTreeSgnCut;
  //delete fullTreeSgnCutTemp;
  delete fullTreeWJetsCut;
  delete fullTreeZttCut;
  delete fullTreeTTJetsCut;
  delete fullTreeWWCut;
  delete fullTreeWZCut;
  delete fullTreeZZCut;
  delete fullTreeDataCut;

  delete fullTreeSgnSSCut;
  //delete fullTreeSgnSSCutTemp;
  delete fullTreeWJetsSSCut;
  delete fullTreeZttSSCut;
  delete fullTreeTTJetsSSCut;
  delete fullTreeWWSSCut;
  delete fullTreeWZSSCut;
  delete fullTreeZZSSCut;
  delete fullTreeDataSSCut;

  delete fullTreeSgnHiMtCut;
  //delete fullTreeSgnHiMtCutTemp;
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
  //TTree* fullTreeSgnCutFTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeDataCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCut_.c_str()) );


  TH1F* hMetF = new TH1F("hMetF","",1,0,1500); 
  fullTreeDataCutF->Draw("event_met_pfmet>>hMetF");
  float NumDataF = hMetF->Integral();
  cout<<"NumDataF "<<NumDataF<<endl;
  hMetF->Reset();
  fullTreeSgnCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumSgnF = hMetF->Integral();
  cout<<"NumSgnF "<<NumSgnF<<endl;
  hMetF->Reset();
  fullTreeWJetsCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWJetsF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZttCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZttF = hMetF->Integral();
  hMetF->Reset();
  fullTreeTTJetsCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumTTJetsF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWWCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWWF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWZCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWZF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZZCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZZF = hMetF->Integral();
  hMetF->Reset();

  std::cout<<"Processing Failing, SS, Low MT"<<std::endl;

  //SS
  TTree* fullTreeDataSSCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCutF = fullTreeSgn->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  //TTree* fullTreeSgnSSCutFTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutSS_.c_str()) );


  fullTreeDataSSCutF->Draw("event_met_pfmet>>hMetF");
  float NumDataSSF = hMetF->Integral();
  cout<<"NumDataSSF "<<NumDataSSF<<endl;
  hMetF->Reset();
  fullTreeSgnSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumSgnSSF = hMetF->Integral();
  cout<<"NumSgnSSF "<<NumSgnSSF<<endl;
  hMetF->Reset();
  fullTreeWJetsSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWJetsSSF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZttSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZttSSF = hMetF->Integral();
  hMetF->Reset();
  fullTreeTTJetsSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumTTJetsSSF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWWSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWWSSF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWZSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWZSSF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZZSSCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZZSSF = hMetF->Integral();
  hMetF->Reset();

  std::cout<<"Processing Failing, OS, High MT"<<std::endl;

  //High MT
  TTree* fullTreeDataHiMtCutF = fullTreeData->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeSgnHiMtCutF = fullTreeSgn->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  //TTree* fullTreeSgnHiMtCutFTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsHiMtCutF = fullTreeWJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttHiMtCutF = fullTreeZtt->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsHiMtCutF = fullTreeTTJets->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWHiMtCutF = fullTreeWW->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZHiMtCutF = fullTreeWZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZHiMtCutF = fullTreeZZ->CopyTree( Form("(%s<%f && %s && %s)",category_.c_str(),cutValue_,bin_.c_str(),additionalCutHiMt_.c_str()) );


  fullTreeDataHiMtCutF->Draw("event_met_pfmet>>hMetF");
  float NumDataHiMtF = hMetF->Integral();
  cout<<"NumDataHiMtF "<<NumDataHiMtF<<" fullTreeDataHiMtCutF "<<fullTreeDataHiMtCutF->GetEntries()<<endl;
  hMetF->Reset();
  fullTreeSgnHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumSgnHiMtF = hMetF->Integral();
  cout<<"NumSgnHiMtF "<<NumSgnHiMtF<<" fullTreeSgnHiMtCutF "<<fullTreeSgnHiMtCutF->GetEntries()<<endl;
  hMetF->Reset();
  fullTreeWJetsHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWJetsHiMtF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZttHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZttHiMtF = hMetF->Integral();
  hMetF->Reset();
  fullTreeTTJetsHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumTTJetsHiMtF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWWHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWWHiMtF = hMetF->Integral();
  hMetF->Reset();
  fullTreeWZHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumWZHiMtF = hMetF->Integral();
  hMetF->Reset();
  fullTreeZZHiMtCutF->Draw("event_met_pfmet>>hMetF","tag_puMCWeightRun2012");
  float NumZZHiMtF = hMetF->Integral();
  hMetF->Reset();;

  delete hMet;

  //Normalization failing

  McP2->cd();

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  //RooDataSet DataDataSetF("DataDataSetF", "dataset for Data fail", RooArgSet(mass), Import( *fullTreeDataCutF ) );
  //std::cout << "data dataset Fail " << DataDataSetF.numEntries() << "  " << std::endl;
  //RooDataHist DataDataHistF("DataDataHistF", "", RooArgSet(mass), DataDataSetF, 1.0);
  RooDataHist DataDataHistF = dataSetProducer(fullTreeDataCutF, mass, 1.0, 1.0, 1.0);
  //float nFail = DataDataHistF.sum(false);

  //RooDataSet dataDataSetSSF("dataDataSetSSF","dataset for signal-pass template SS", RooArgSet(mass), Import( *fullTreeDataSSCutF ) );
  //RooDataHist dataDataHistSSF("dataDataHistSSF", "", RooArgSet(mass), dataDataSetSSF, 1.0);
  RooDataHist dataDataHistSSF = dataSetProducer(fullTreeDataSSCutF, mass, 1.0, 1.0, 1.0);

  //RooDataSet dataDataSetHiMtF("dataDataSetHiMtF","dataset for signal-pass template HiMt", RooArgSet(mass), Import( *fullTreeDataHiMtCutF ) );
  //RooDataHist dataDataHistHiMtF("dataDataHistHiMtF", "", RooArgSet(mass), dataDataSetHiMtF, 1.0);
  RooDataHist dataDataHistHiMtF = dataSetProducer(fullTreeDataHiMtCutF, mass, 1.0, 1.0, 1.0);

  delete fullTreeDataCutF;
  delete fullTreeDataSSCutF;
  delete fullTreeDataHiMtCutF;

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist sgnDataHistF = dataSetProducer(fullTreeSgnCutF, mass, NumSgnF, ScaleFactorSgn, 1.0);
  RooHistPdf sgnTemplatePdfF("sgnTemplatePdfF", "", RooArgSet(mass), sgnDataHistF, 4);

  RooDataHist sgnDataHistSSF = dataSetProducer(fullTreeSgnSSCutF, mass, NumSgnSSF, ScaleFactorSgn, -1.0);

  RooDataHist sgnDataHistHiMtF = dataSetProducer(fullTreeSgnHiMtCutF, mass, NumSgnHiMtF, ScaleFactorSgn, -1.0);

  delete fullTreeSgnCutF;
  delete fullTreeSgnSSCutF;
  delete fullTreeSgnHiMtCutF;

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wjetsDataHistF = dataSetProducer(fullTreeWJetsCutF, mass, NumWJetsF, ScaleFactorWJets, 1.0);
  RooHistPdf wjetsHistPdfF("wjetsHistPdfF", "", RooArgSet(mass), wjetsDataHistF, 4);

  RooDataHist wjetsDataHistSSF = dataSetProducer(fullTreeWJetsSSCutF, mass, NumWJetsSSF, ScaleFactorWJets, -1.0);

  RooDataHist wjetsDataHistHiMtF = dataSetProducer(fullTreeWJetsHiMtCutF, mass, NumWJetsHiMtF, ScaleFactorWJets, -1.0);

  delete fullTreeWJetsCutF;
  delete fullTreeWJetsSSCutF;
  delete fullTreeWJetsHiMtCutF;

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zttDataHistF = dataSetProducer(fullTreeZttCutF, mass, NumZttF, ScaleFactorZtt, 1.0);
  RooHistPdf zttPdfF("zttPdfF", "", RooArgSet(mass), zttDataHistF, 4);  

  RooDataHist zttDataHistSSF = dataSetProducer(fullTreeZttSSCutF, mass, NumZttSSF, ScaleFactorZtt, -1.0);

  RooDataHist zttDataHistHiMtF = dataSetProducer(fullTreeZttHiMtCutF, mass, NumZttHiMtF, ScaleFactorZtt, -1.0);

  delete fullTreeZttCutF;
  delete fullTreeZttSSCutF;
  delete fullTreeZttHiMtCutF;

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist ttjetsDataHistF = dataSetProducer(fullTreeTTJetsCutF, mass, NumTTJetsF, ScaleFactorTTJets, 1.0);
  RooHistPdf ttjetsPdfF("ttjetsPdfF", "", RooArgSet(mass), ttjetsDataHistF, 4);

  RooDataHist ttjetsDataHistSSF = dataSetProducer(fullTreeTTJetsSSCutF, mass, NumTTJetsSSF, ScaleFactorTTJets, -1.0);

  RooDataHist ttjetsDataHistHiMtF = dataSetProducer(fullTreeTTJetsHiMtCutF, mass, NumTTJetsHiMtF, ScaleFactorTTJets, -1.0);

  delete fullTreeTTJetsCutF;
  delete fullTreeTTJetsSSCutF;
  delete fullTreeTTJetsHiMtCutF;

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wwDataHistF = dataSetProducer(fullTreeWWCutF, mass, NumWWF, ScaleFactorWW, 1.0);
  RooHistPdf wwPdfF("wwPdfF", "", RooArgSet(mass), wwDataHistF, 4);

  mass.setBins( 50 );
  RooDataHist wwDataHistSSF = dataSetProducer(fullTreeWWSSCutF, mass, NumWWSSF, ScaleFactorWW, -1.0);

  mass.setBins( 50 );
  RooDataHist wwDataHistHiMtF = dataSetProducer(fullTreeWWHiMtCutF, mass, NumWWHiMtF, ScaleFactorWW, -1.0);

  delete fullTreeWWCutF;
  delete fullTreeWWSSCutF;
  delete fullTreeWWHiMtCutF;

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist wzDataHistF = dataSetProducer(fullTreeWZCutF, mass, NumWZF, ScaleFactorWZ, 1.0);
  RooHistPdf wzPdfF("wzPdfF", "", RooArgSet(mass), wzDataHistF, 4);

  mass.setBins( 50 );
  RooDataHist wzDataHistSSF = dataSetProducer(fullTreeWZSSCutF, mass, NumWZSSF, ScaleFactorWZ, -1.0);

  mass.setBins( 50 );
  RooDataHist wzDataHistHiMtF = dataSetProducer(fullTreeWZHiMtCutF, mass, NumWZHiMtF, ScaleFactorWZ, -1.0);

  delete fullTreeWZCutF;
  delete fullTreeWZSSCutF;
  delete fullTreeWZHiMtCutF;

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  mass.setBins( 50 );
  RooDataHist zzDataHistF = dataSetProducer(fullTreeZZCutF, mass, NumZZF, ScaleFactorZZ, 1.0);
  RooHistPdf zzPdfF("zzPdfF", "", RooArgSet(mass), zzDataHistF, 4);

  mass.setBins( 50 );
  RooDataHist zzDataHistSSF = dataSetProducer(fullTreeZZSSCutF, mass, NumZZSSF, ScaleFactorZZ, -1.0);

  mass.setBins( 50 );
  RooDataHist zzDataHistHiMtF = dataSetProducer(fullTreeZZHiMtCutF, mass, NumZZHiMtF, ScaleFactorZZ, -1.0);

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
  float exWJetsF = NumWJetsF/NumWJetsHiMtF;
  cout<<"exWJetsF "<<exWJetsF<<endl;
  RooDataHist dataDataHistHiMtScaledF("dataDataHistHiMtScaledF", "", RooArgSet(mass), dataDataHistHiMtF, exWJetsF);
  //dataDataHistHiMtScaledF.plotOn(massFrame2,MarkerColor(kViolet));
  RooHistPdf WJetsDataDrivenHistPdfF("WJetsDataDrivenHistPdfF", "", RooArgSet(mass), dataDataHistHiMtScaledF, 4);

  // QCD
  dataDataHistSSF.add(sgnDataHistSSF);
  dataDataHistSSF.add(wjetsDataHistSSF);
  dataDataHistSSF.add(ttjetsDataHistSSF);
  dataDataHistSSF.add(zttDataHistSSF);
  dataDataHistSSF.add(wwDataHistSSF);
  dataDataHistSSF.add(wzDataHistSSF);
  dataDataHistSSF.add(zzDataHistSSF);

  RooHistPdf QCDHistPdfF("QCDHistPdfF", "", RooArgSet(mass), dataDataHistSSF, 4);

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

  DataDataHistP.reset();
  DataDataHistP.add(((RooDataHist)dataDataHistSS)); //QCD
  DataDataHistP.add(((RooDataHist)zttDataHist));
  DataDataHistP.add(((RooDataHist)sgnDataHist));
  DataDataHistP.add(((RooDataHist)dataDataHistHiMtScaled)); //WJets
  DataDataHistP.add(((RooDataHist)ttjetsDataHist));
  DataDataHistP.add(((RooDataHist)wwDataHist));
  DataDataHistP.add(((RooDataHist)wzDataHist));
  DataDataHistP.add(((RooDataHist)zzDataHist));
  //float nPassTemplate =  DataDataHistP.sum(false);

  DataDataHistF.reset();
  DataDataHistF.add(((RooDataHist)dataDataHistSSF)); //QCD
  DataDataHistF.add(((RooDataHist)zttDataHistF));
  DataDataHistF.add(((RooDataHist)sgnDataHistF));
  DataDataHistF.add(((RooDataHist)wjetsDataHistF)); //WJets
  DataDataHistF.add(((RooDataHist)ttjetsDataHistF));
  DataDataHistF.add(((RooDataHist)wwDataHistF));
  DataDataHistF.add(((RooDataHist)wzDataHistF));
  DataDataHistF.add(((RooDataHist)zzDataHistF));
  //float nFailTemplate =  DataDataHistF.sum(false);

  //////////////////////////////////////////////////////////////// SimFit ////////////////////////////////////////////////////////////////

  //Create pdfs for passing
  RooRealVar CoeffSgnP("CoeffSgnP","",0,10000000);
  RooRealVar CoeffZttP("CoeffZttP","",0,10000000);
  RooRealVar CoeffWJetsP("CoeffWJetsP","",0,10000000);
  RooRealVar CoeffTTJetsP("CoeffTTJetsP","",0,10000000);
  RooRealVar CoeffQCDP("CoeffQCDP","",0,10000000);
  RooRealVar CoeffWWP("CoeffWWP","",0,10000000);
  RooRealVar CoeffWZP("CoeffWZP","",0,10000000);
  RooRealVar CoeffZZP("CoeffZZP","",0,10000000);

  RooAddPdf DataModelP("DataModelP", "", RooArgList(sgnTemplatePdf,zttPdf,WJetsHistPdfP_C,ttjetsPdf,QCDHistPdfP_C,wwPdf,wzPdf,zzPdf), RooArgList(CoeffSgnP/*DataNumSgnP*/,CoeffZttP,CoeffWJetsP,CoeffTTJetsP,CoeffQCDP,CoeffWWP,CoeffWZP,CoeffZZP));

  //Create pdfs for failing
  RooRealVar CoeffSgnF("CoeffSgnF","",0,10000000);
  RooRealVar CoeffZttF("CoeffZttF","",0,10000000);
  RooRealVar CoeffWJetsF("CoeffWJetsF","",0,10000000);
  RooRealVar CoeffTTJetsF("CoeffTTJetsF","",0,10000000);
  RooRealVar CoeffQCDF("CoeffQCDF","",0,10000000);
  RooRealVar CoeffWWF("CoeffWWF","",0,10000000);
  RooRealVar CoeffWZF("CoeffWZF","",0,10000000);
  RooRealVar CoeffZZF("CoeffZZF","",0,10000000);

  RooAddPdf DataModelF("DataModelF", "", RooArgList(sgnTemplatePdfF,zttPdfF,WJetsHistPdfF_C,ttjetsPdfF,QCDHistPdfF_C,wwPdfF,wzPdfF,zzPdfF), RooArgList(CoeffSgnF/*DataNumSgnP*/,CoeffZttF,CoeffWJetsF,CoeffTTJetsF,CoeffQCDF,CoeffWWF,CoeffWZF,CoeffZZF));

  mass.setBins(nBins_);

  RooCategory category("category","category") ;
  category.defineType("pass") ;
  category.defineType("fail") ;

  // binned combined dataset
  RooDataHist DataCombData("DataCombData","combined data",mass,Index(category),Import("pass", *(DataDataHistP.createHistogram("histoDataP",mass)) ) ,Import("fail", *(DataDataHistF.createHistogram("histoDataF",mass))), Weight(0.5) ) ;
  // unbinned combined dataset
  //RooDataSet DataCombDataUnBinned("DataCombDataUnBinned","combined data",mass,Index(category),Import("pass", DataDataHistP) ,Import("fail",DataDataHistF), Weight(0.5) ) ;

  RooSimultaneous DataSimPdf("DataSimPdf","simultaneous pdf",category) ;
  DataSimPdf.addPdf(DataModelP,"pass") ;
  DataSimPdf.addPdf(DataModelF,"fail") ;

  //mass.setBins( 10000, "fft" );
  mass.setBins( nBins_ );
  RooFitResult* ResDataCombinedFit =  0;
  if(doBinned_)  ResDataCombinedFit = DataSimPdf.fitTo(DataCombData , Extended(1), Minos(1), Save(1), NumCPU(4), /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/  SumW2Error(1));
  //else ResDataCombinedFit = DataSimPdf.fitTo(DataCombDataUnBinned , Extended(1), Minos(1), Save(1), NumCPU(4),  /*ExternalConstraints( RooArgSet(alfaSgn_CPdf,nSgn_CPdf) )*/ SumW2Error(1));

  RooArgSet DataFitParam(ResDataCombinedFit->floatParsFinal());

  string theSample = "MC";

  RooPlot* DataFrameP = mass.frame(Bins(40),Title(Form("CMS Preliminary 2012 #sqrt{s}=8 TeV %s  L=%.0f pb^{-1}: passing probe",theSample.c_str(),Lumi_)));
  DataCombData.plotOn(DataFrameP,Cut("category==category::pass"),Name("dataP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), LineColor(kBlack),Name("modelP"));
  DataSimPdf.plotOn(DataFrameP,Slice(category,"pass"), ProjWData(category,DataCombData), Components("sgnTemplatePdf"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyP"));
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
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("sgnTemplatePdfF"), LineColor(kBlue), LineStyle(kSolid),Name("signal onlyF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("zttPdfF"), LineColor(kRed), LineStyle(kSolid),Name("ZttF"));
  DataSimPdf.plotOn(DataFrameF,Slice(category,"fail"), ProjWData(category,DataCombData), Components("WJetsDataDrivenHistPdfF"), LineColor(kGreen), LineStyle(kSolid),Name("WJetsF"));
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
  leg1->AddEntry("signal onlyP","Signal only", "L");
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
  leg2->AddEntry("signal onlyF","Signal only", "L");
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

  RooFormulaVar DataEffFit("DataEffFit","DataEffFit","CoeffSgnP/(CoeffSgnP+CoeffSgnF)", RooArgSet(CoeffSgnP,CoeffSgnP,CoeffSgnF));
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

void calcutateFitClosure(){

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

