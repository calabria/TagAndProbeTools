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

TH1F * massPlotProducer(TTree * fitter_tree){

	Float_t mass1;
	fitter_tree->SetBranchAddress("mass", &mass1);
	TH1F *hMass = new TH1F("mass","mass distribution",200,0,200);
	Int_t nentries = (Int_t)fitter_tree->GetEntries();
	for (Int_t i=0; i<nentries; i++) {
		fitter_tree->GetEntry(i);
		hMass->Fill(mass1);
	}

	return hMass;

}

RooDataHist dataSetProducer(TTree * fullTreeSgnCut, RooRealVar mass, float NumSgnP, float ScaleFactorSgn, float weight){

  	RooDataSet sgnDataSet("sgnDataSet", "dataset for signal-pass template", RooArgSet(mass), Import( *fullTreeSgnCut ) );
  	//if(NumSgnP*ScaleFactorSgn != 1) sgnDataSet.reduce(EventRange(1,(int)( NumSgnP*ScaleFactorSgn)));
        cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
  	RooDataHist sgnDataHist("sgnDataHist", "", RooArgSet(mass), sgnDataSet, weight*ScaleFactorSgn);

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
			const string additionalCut_      = "pt > 0 && pair_charge == 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_muPFIsolation < 0.1",
			const string additionalCutSS_    = "pt > 0 && pair_charge != 0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_muPFIsolation < 0.1",
			const string additionalCutHiMt_  = "pt > 0 && pair_charge == 0 && tag_Mt > 60 && tag_IsoMu24_eta2p1 && tag_muPFIsolation < 0.1",
			float scale_                     = 0.0
			){

  string path = "./InputFile/";
  string pathData = "./InputFile/";

  // signal
  TFile fsgn((path + "testTagAndProbe_DYToMuMu_TNP_TNP.root").c_str());
  TTree *fullTreeSgn = (TTree*)fsgn.Get((tnp_+"/fitter_tree").c_str());

  // MC truth Efficiency

  TH1F* hS           = new TH1F("hS","",1,0,1500);
  TH1F* hSP          = new TH1F("hSP","",1,0,1500);

  fullTreeSgn->Draw("mass>>hS",Form("1.0*(%s && mass>%f && mass<%f && pair_charge==0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_muPFIsolation < 0.1 )",bin_.c_str(),xLow_,xHigh_));
  double SGNtrue = hS->Integral();
  fullTreeSgn->Draw("mass>>hSP",Form("1.0*(%s && %s>=%f && mass>%f && mass<%f && pair_charge==0 && tag_Mt < 40 && tag_IsoMu24_eta2p1 && tag_muPFIsolation < 0.1 )",bin_.c_str(),category_.c_str(),cutValue_,xLow_,xHigh_));
  double SGNtruePass = hSP->Integral();

  double McTruthEff    = SGNtruePass/SGNtrue;
  double BinomialError = TMath::Sqrt(SGNtruePass/SGNtrue*(1-SGNtruePass/SGNtrue)/SGNtrue);
  
  cout << bin_.c_str() << " ==> MCTRUTH: " << McTruthEff << " +/- " << BinomialError << endl;

  delete hS; delete hSP;

}

void calcutateFitMc(){

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
