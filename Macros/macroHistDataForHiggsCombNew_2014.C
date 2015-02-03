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
#include "TLorentzVector.h"
#include <vector>
#include "THStack.h"
using namespace std;

std::vector<TH1F *> CreateESMassShape(TTree * fullTreeSgn, TString processName, double scaleFactor){


  std::vector<TH1F *> hTemp;
  TH1F *h1 = new TH1F ("hTemp0","hTemp distribution",50,70,120);
  TH1F *h2 = new TH1F ("hTemp1","hTemp distribution",50,70,120);
  TH1F *h3 = new TH1F ("hTemp2","hTemp distribution",50,70,120);
  TH1F *h4 = new TH1F ("hTemp3","hTemp distribution",50,70,120);
  hTemp.push_back (h1);
  hTemp.push_back (h2);
  hTemp.push_back (h3);
  hTemp.push_back (h4);


  Float_t pt, eta, phi, En;
  Float_t tag_pt, tag_eta, tag_phi, tag_En;
 
  fullTreeSgn->SetBranchAddress("pt",&pt);
  fullTreeSgn->SetBranchAddress("tag_pt",&tag_pt);
  fullTreeSgn->SetBranchAddress("eta",&eta);
  fullTreeSgn->SetBranchAddress("tag_eta",&tag_eta);
  fullTreeSgn->SetBranchAddress("phi",&phi);
  fullTreeSgn->SetBranchAddress("tag_phi",&tag_phi);
  fullTreeSgn->SetBranchAddress("En",&En);
  fullTreeSgn->SetBranchAddress("tag_En",&tag_En);
  
 
  Int_t nentries = (Int_t)fullTreeSgn->GetEntries();

  TLorentzVector Tag; Tag.SetPtEtaPhiE(0.,0.,0.,0.);
  TLorentzVector Probe; Probe.SetPtEtaPhiE(0.,0.,0.,0.);

  TLorentzVector Tag_Up; Tag_Up.SetPtEtaPhiE(0.,0.,0.,0.);
  TLorentzVector Probe_Up; Probe_Up.SetPtEtaPhiE(0.,0.,0.,0.);

  TLorentzVector Tag_Down; Tag_Down.SetPtEtaPhiE(0.,0.,0.,0.);
  TLorentzVector Probe_Down; Probe_Down.SetPtEtaPhiE(0.,0.,0.,0.);

  TLorentzVector Pair_TauUp, Pair_TauDown, Pair_MuUp, Pair_MuDown, Pair;
 
  for (Int_t i=0; i<nentries; i++) {
    fullTreeSgn->GetEntry(i);

    // std::cout<<" ********************************px************************ tau "<<pt<<" mu  "<<tag_pt<<std::endl;
    Tag.SetPtEtaPhiE(tag_pt, tag_eta,tag_phi, tag_En);
    Probe.SetPtEtaPhiE(pt, eta, phi, En);

    Tag_Up.SetPtEtaPhiE(tag_pt*(1.03), tag_eta,tag_phi, (1.03)*tag_En);
    Probe_Up.SetPtEtaPhiE(pt*1.01, eta, phi, En*1.01);

    Tag_Down.SetPtEtaPhiE(tag_pt*0.97, tag_eta,tag_phi, tag_En*0.97);
    Probe_Down.SetPtEtaPhiE(pt*0.99, eta, phi, En*0.99);

    Pair = Tag + Probe;
   
    Pair_TauUp = Probe_Up + Tag;
    Pair_TauDown = Probe_Down + Tag;

    Pair_MuUp = Tag_Up + Probe;
    Pair_MuDown = Tag_Down + Probe;

    // cout<<" Mass   "<<Pair.M()<<" mass Up  "<<Pair_TauUp.M()<<endl;
    hTemp[0]->Fill(Pair_TauUp.M());
    hTemp[1]->Fill(Pair_TauDown.M());
    hTemp[2]->Fill(Pair_MuUp.M());
    hTemp[3]->Fill(Pair_MuDown.M());
   
  }

  hTemp[0]->Scale(scaleFactor);
  hTemp[1]->Scale(scaleFactor);
  hTemp[2]->Scale(scaleFactor);
  hTemp[3]->Scale(scaleFactor);
  

  hTemp[0]->Write(processName+"_CMS_scale_t_mutau_8TeVUp");
  hTemp[1]->Write(processName+"_CMS_scale_t_mutau_8TeVDown");
  hTemp[2]->Write(processName+"_CMS_scale_m_mutau_8TeVUp");
  hTemp[3]->Write(processName+"_CMS_scale_m_mutau_8TeVDown");

  return hTemp;
  //  delete hTemp;

  
}

TH1F * histoProducer(TTree * fullTreeSgnCut, float ScaleFactorSgn, float weight, std::string var, int min, int max, int num){

	TH1F* hMass = new TH1F("hMass","",num,min,max);
	std::string varTmp = var + ">>hMass"; 
  	if(ScaleFactorSgn != 1) {
		fullTreeSgnCut->Draw(varTmp.c_str(),"tag_puMCWeightRun2012");
		hMass->Scale(ScaleFactorSgn*weight*0.987883333*0.985533333*0.9548436313);
	}
	else fullTreeSgnCut->Draw(varTmp.c_str());
        //cout<<"sgnDataSet "<<sgnDataSet.sumEntries()<<endl;
        //cout<<"ScaleFactorSgn "<<ScaleFactorSgn<<" sgnDataSet "<<fullTreeSgnCut->GetEntries()<<endl;
        //cout<<"scaled "<<hMass->Integral(0,51)<<endl;

	return hMass;

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
  	const std::string additionalCut_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5 && tag_muMVAMet < 25)",
  	const std::string additionalCutSS_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5 && tag_muMVAMet < 25)",
  	const std::string additionalCutHiMt_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5 && tag_muMVAMet < 25)"
	){

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

  //Low MT No Discriminator
  TTree* fullTreeSgnCut2 = fullTreeSgn->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsCut2 = fullTreeWJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttCut2 = fullTreeZtt->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsCut2 = fullTreeTTJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWCut2 = fullTreeWW->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZCut2 = fullTreeWZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZCut2 = fullTreeZZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeDataCut2 = fullTreeData->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutHiMt_.c_str()) );

  //SS No Cut
  TTree* fullTreeDataSSCut2 = fullTreeData->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut2 = fullTreeSgn->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut2 = fullTreeWJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut2 = fullTreeZtt->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut2 = fullTreeTTJets->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut2 = fullTreeWW->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut2 = fullTreeWZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut2 = fullTreeZZ->CopyTree( Form("(%s && %s)",bin_.c_str(),additionalCutSS_.c_str()) );

  McP->cd();

  //Normalization passing

  //float Lumi_ = 19484.55;
  //  float Lumi_ = 15829.435;
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

  std::map<std::string, TTree*> reducedTrees;
  std::map<std::string, TTree*> reducedTreesSS;
  std::map<std::string, TTree*> reducedTreesHiMt;

  reducedTrees.insert ( std::pair<std::string, TTree*>("data", fullTreeDataCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZCut) );

  reducedTrees.insert ( std::pair<std::string, TTree*>("zmmMatch", fullTreeSgnCutTemp) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zmmJets", fullTreeSgnCutTempJet) );

  reducedTreesSS.insert ( std::pair<std::string, TTree*>("data", fullTreeDataSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZSSCut) );
  reducedTreesSS.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZSSCut) );

  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("data", fullTreeDataHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZHiMtCut) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZHiMtCut) );

  std::map<std::string, float> scaleFactors;

  scaleFactors.insert ( std::pair<std::string, float>("zmm", ScaleFactorSgn) );
  scaleFactors.insert ( std::pair<std::string, float>("ztt", ScaleFactorZtt) );
  scaleFactors.insert ( std::pair<std::string, float>("wjets", ScaleFactorWJets) );
  scaleFactors.insert ( std::pair<std::string, float>("ttjets", ScaleFactorTTJets) );
  scaleFactors.insert ( std::pair<std::string, float>("ww", ScaleFactorWW) );
  scaleFactors.insert ( std::pair<std::string, float>("wz", ScaleFactorWZ) );
  scaleFactors.insert ( std::pair<std::string, float>("zz", ScaleFactorZZ) );

  ///////////////////////////////////////

  string TreeNameP = "./HistosForHiggsCombine_2014/histForSaving_"+tnp_+"_"+category_;
  string TreeNamePass = "./HistosForHiggsCombine_2014/histForSaving_"+tnp_+"_"+category_;
  TFile *SaveP = new TFile(Form("%s%s.root",TreeNamePass.c_str(), binCenter_.c_str()),"RECREATE");

  SaveP->mkdir("pass");
  SaveP->cd("pass");

  TH1F* hDummy1 = new TH1F("hDummy","",50,70,120); 
  hDummy1->SetBinContent(1,1);
  hDummy1->Write("ggH");

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoProducer(reducedTrees["data"], 1, 1, "mass", 70, 120, 50);
  TH1F * hDataDataHistPSS = histoProducer(reducedTreesSS["data"], 1, 1, "mass", 70, 120, 50);
  TH1F * hDataDataHistPHiMt = histoProducer(reducedTreesHiMt["data"], 1, 1, "mass", 70, 120, 50);
 
  hDataDataHistP->Write("data_obs");
  std::vector<TH1F*> dataShift = CreateESMassShape(fullTreeDataCut,"data_obs",1.0);

  TH1F * hDataDataHistP2 = histoProducer(fullTreeDataCut2, 1, 1, "mass", 70, 120, 50);
  TH1F * hDataDataHistPSS2 = histoProducer(fullTreeDataSSCut2, 1, 1, "mass", 70, 120, 50);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoProducer(reducedTrees["zmm"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistTemp = histoProducer(reducedTrees["zmmMatch"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistTempNo = histoProducer(reducedTrees["zmmJets"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistSS = histoProducer(reducedTreesSS["zmm"], scaleFactors["zmm"], -1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistHiMt = histoProducer(reducedTreesHiMt["zmm"], scaleFactors["zmm"], -1, "mass", 70, 120, 50);

  hsgnDataHistTemp->Write("ZL"); //signal matched
  hsgnDataHist->Write("zmumuAll");
  hsgnDataHistTempNo->Write("zmumuJets");

  std::vector<TH1F*> zmumuShift = CreateESMassShape(fullTreeSgnCutTemp,"ZL",ScaleFactorSgn);

  TH1F * hsgnDataHist2 = histoProducer(fullTreeSgnCut2, ScaleFactorSgn, 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistSS2 = histoProducer(fullTreeSgnSSCut2, ScaleFactorSgn, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoProducer(reducedTrees["wjets"], scaleFactors["wjets"], 1, "mass", 70, 120, 50);
  TH1F * hwjetsDataHistSS = histoProducer(reducedTreesSS["wjets"], scaleFactors["wjets"], -1, "mass", 70, 120, 50);
  TH1F * hwjetsDataHistHiMt = histoProducer(reducedTreesHiMt["wjets"], scaleFactors["wjets"], 1, "mass", 70, 120, 50);

  hwjetsDataHist->Write("wjets");

  std::vector<TH1F*> WjetsShift = CreateESMassShape(fullTreeWJetsCut,"wjets",ScaleFactorWJets);

  TH1F * hwjetsDataHist2 = histoProducer(fullTreeWJetsCut2, ScaleFactorWJets, 1, "mass", 70, 120, 50);
  TH1F * hwjetsDataHistSS2 = histoProducer(fullTreeWJetsSSCut2, ScaleFactorWJets, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoProducer(reducedTrees["ztt"], scaleFactors["ztt"], 1, "mass", 70, 120, 50);
  TH1F * hzttDataHistSS = histoProducer(reducedTreesSS["ztt"], scaleFactors["ztt"], -1, "mass", 70, 120, 50);
  TH1F * hzttDataHistHiMt = histoProducer(reducedTreesHiMt["ztt"], scaleFactors["ztt"], -1, "mass", 70, 120, 50);

  hzttDataHist->Write("ZTT");

  std::vector<TH1F*> ztautauhift = CreateESMassShape(fullTreeZttCut,"ZTT",ScaleFactorZtt);

  TH1F * hzttDataHist2 = histoProducer(fullTreeZttCut2, ScaleFactorZtt, 1, "mass", 70, 120, 50);
  TH1F * hzttDataHistSS2 = histoProducer(fullTreeZttSSCut2, ScaleFactorZtt, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoProducer(reducedTrees["ttjets"], scaleFactors["ttjets"], 1, "mass", 70, 120, 50);
  TH1F * httjetsDataHistSS = histoProducer(reducedTreesSS["ttjets"], scaleFactors["ttjets"], -1, "mass", 70, 120, 50);
  TH1F * httjetsDataHistHiMt = histoProducer(reducedTreesHiMt["ttjets"], scaleFactors["ttjets"], -1, "mass", 70, 120, 50);

  httjetsDataHist->Write("TT");

  std::vector<TH1F*> ttjetsShift = CreateESMassShape(fullTreeTTJetsCut,"TT",ScaleFactorTTJets);

  TH1F * httjetsDataHist2 = histoProducer(fullTreeTTJetsCut2, ScaleFactorTTJets, 1, "mass", 70, 120, 50);
  TH1F * httjetsDataHistSS2 = histoProducer(fullTreeTTJetsSSCut2, ScaleFactorTTJets, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoProducer(reducedTrees["ww"], scaleFactors["ww"], 1, "mass", 70, 120, 50);
  TH1F * hwwDataHistSS = histoProducer(reducedTreesSS["ww"], scaleFactors["ww"], -1, "mass", 70, 120, 50);
  TH1F * hwwDataHistHiMt = histoProducer(reducedTreesHiMt["ww"], scaleFactors["ww"], -1, "mass", 70, 120, 50);

  hwwDataHist->Write("ww");

  std::vector<TH1F*> wwShift = CreateESMassShape(fullTreeWWCut,"ww",ScaleFactorWW);

  TH1F * hwwDataHist2 = histoProducer(fullTreeWWCut2, ScaleFactorWW, 1, "mass", 70, 120, 50);
  TH1F * hwwDataHistSS2 = histoProducer(fullTreeWWSSCut2, ScaleFactorWW, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoProducer(reducedTrees["wz"], scaleFactors["wz"], 1, "mass", 70, 120, 50);
  TH1F * hwzDataHistSS = histoProducer(reducedTreesSS["wz"], scaleFactors["wz"], -1, "mass", 70, 120, 50);
  TH1F * hwzDataHistHiMt = histoProducer(reducedTreesHiMt["wz"], scaleFactors["wz"], -1, "mass", 70, 120, 50);

  hwzDataHist->Write("wz");

  std::vector<TH1F*> wzShift = CreateESMassShape(fullTreeWZCut,"wz",ScaleFactorWZ);

  TH1F * hwzDataHist2 = histoProducer(fullTreeWZCut2, ScaleFactorWZ, 1, "mass", 70, 120, 50);
  TH1F * hwzDataHistSS2 = histoProducer(fullTreeWZSSCut2, ScaleFactorWZ, -1, "mass", 70, 120, 50);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoProducer(reducedTrees["zz"], scaleFactors["zz"], 1, "mass", 70, 120, 50);
  TH1F * hzzDataHistSS = histoProducer(reducedTreesSS["zz"], scaleFactors["zz"], -1, "mass", 70, 120, 50);
  TH1F * hzzDataHistHiMt = histoProducer(reducedTreesHiMt["zz"], scaleFactors["zz"], -1, "mass", 70, 120, 50);

  hzzDataHist->Write("zz");

  std::vector<TH1F*> zzShift = CreateESMassShape(fullTreeZZCut,"zz",ScaleFactorZZ);

  TH1F * hzzDataHist2 = histoProducer(fullTreeZZCut2, ScaleFactorZZ, 1, "mass", 70, 120, 50);
  TH1F * hzzDataHistSS2 = histoProducer(fullTreeZZSSCut2, ScaleFactorZZ, -1, "mass", 70, 120, 50);


  hDataDataHistP->Sumw2(); 
  hDataDataHistPSS->Sumw2(); 
  hDataDataHistPHiMt->Sumw2(); 

  hsgnDataHist->Sumw2(); 
  hsgnDataHistSS->Sumw2(); 
  hsgnDataHistHiMt->Sumw2(); 
  hsgnDataHistTemp->Sumw2(); 
  hsgnDataHistTempNo->Sumw2(); 

  hwjetsDataHist->Sumw2(); 
  hwjetsDataHistSS->Sumw2();
  hwjetsDataHistHiMt->Sumw2();

  hzttDataHist->Sumw2(); 
  hzttDataHistSS->Sumw2();
  hzttDataHistHiMt->Sumw2();

  httjetsDataHist->Sumw2();
  httjetsDataHistSS->Sumw2();
  httjetsDataHistHiMt->Sumw2();

  hwwDataHist->Sumw2();
  hwwDataHistSS->Sumw2(); 
  hwwDataHistHiMt->Sumw2(); 

  hwzDataHist->Sumw2();
  hwzDataHistSS->Sumw2(); 
  hwzDataHistHiMt->Sumw2(); 

  hzzDataHist->Sumw2(); 
  hzzDataHistSS->Sumw2(); 
  hzzDataHistHiMt->Sumw2();

  /////////////////////////////////////////////////////////////////////////////////////

  hDataDataHistP2->Sumw2(); 
  hDataDataHistPSS2->Sumw2(); 

  hsgnDataHist2->Sumw2(); 
  hsgnDataHistSS2->Sumw2(); 

  hwjetsDataHist2->Sumw2(); 
  hwjetsDataHistSS2->Sumw2(); 

  hzttDataHist2->Sumw2(); 
  hzttDataHistSS2->Sumw2();

  httjetsDataHist2->Sumw2(); 
  httjetsDataHistSS2->Sumw2();

  hwwDataHist2->Sumw2(); 
  hwwDataHistSS2->Sumw2(); 

  hwzDataHist2->Sumw2(); 
  hwzDataHistSS2->Sumw2();

  hzzDataHist2->Sumw2(); 
  hzzDataHistSS2->Sumw2();

  /////////////////////////////////////////////////////////////////////////////////////

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
  hDataDataHistP2->Divide(hsgnDataHist2); //at this point histogram with scale factors
  hDataDataHistP2->Write("scaleFactorsHiMt");

  /////////////////////////////////////////////////////////////////////////////////////

  TH1F *h1VV = new TH1F("h1VV", "hVV", 50, 70, 120);
  h1VV->Add(hzzDataHist);
  h1VV->Add(hwzDataHist);
  h1VV->Add(hzzDataHist);
  h1VV->Write("VV");
  
  TH1F *hVVTES_Up = new TH1F("hVVTES_Up", "hVVTES_Up", 50, 70, 120);
  hVVTES_Up->Add(zzShift[0]);
  hVVTES_Up->Add(wzShift[0]);
  hVVTES_Up->Add(wwShift[0]);
  hVVTES_Up->Write("VV_CMS_scale_t_mutau_8TeVUp");

  TH1F *hVVTES_Down = new TH1F("hVVTES_Down", "hVVTES_Down", 50, 70, 120);
  hVVTES_Down->Add(zzShift[1]);
  hVVTES_Down->Add(wzShift[1]);
  hVVTES_Down->Add(wwShift[1]);
  hVVTES_Down->Write("VV_CMS_scale_t_mutau_8TeVDown");

  TH1F *hVVMMS_Up = new TH1F("hVVMMS_Up", "hVVMMS_Up", 50, 70, 120);
  hVVMMS_Up->Add(zzShift[2]);
  hVVMMS_Up->Add(wzShift[2]);
  hVVMMS_Up->Add(wwShift[2]);
  hVVMMS_Up->Write("VV_CMS_scale_m_mutau_8TeVUp");

  TH1F *hVVMMS_Down = new TH1F("hVVMMS_Down", "hVVMMS_Down", 50, 70, 120);
  hVVMMS_Down->Add(zzShift[3]);
  hVVMMS_Down->Add(wzShift[3]);
  hVVMMS_Down->Add(wwShift[3]);
  hVVMMS_Down->Write("VV_CMS_scale_m_mutau_8TeVDown");

  //WJets from SideBand
  hsgnDataHistHiMt->Sumw2();
  hDataDataHistP2->Sumw2();
  hsgnDataHistHiMt->Multiply(hDataDataHistP2); //rescaling signal

  hDataDataHistPHiMt->Add(hsgnDataHistHiMt);
  hDataDataHistPHiMt->Add(hzttDataHistHiMt);
  hDataDataHistPHiMt->Add(httjetsDataHistHiMt);
  hDataDataHistPHiMt->Add(hwwDataHistHiMt);
  hDataDataHistPHiMt->Add(hwzDataHistHiMt);
  hDataDataHistPHiMt->Add(hzzDataHistHiMt);
  float expW = hwjetsDataHist->Integral(0,hwjetsDataHist->GetSize()) / hwjetsDataHistHiMt->Integral(0,hwjetsDataHistHiMt->GetSize());
  cout<<"expW "<<expW<<endl;
  hDataDataHistPHiMt->Scale(expW);
  hDataDataHistPHiMt->Write("wjetsDataDriven");

  float numWJetsDataDriven = hDataDataHistPHiMt->Integral(0,hDataDataHistPHiMt->GetSize());
  float normWJetsMC = numWJetsDataDriven / hwjetsDataHist->Integral(0,hwjetsDataHist->GetSize());
  TH1F * wjetsNorm = (TH1F*)hwjetsDataHist->Clone();
  wjetsNorm->Scale(normWJetsMC);
  wjetsNorm->Write("W");

  //QCD datadriven
  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  for(int i = 1; i < hDataDataHistPSS->GetSize(); i++){

      /*double err_zz  =  (hzz_ss->GetBinError(i))*(hzz_ss->GetBinError(i));
      double err_wz  =  (hwz_ss->GetBinError(i))*(hwz_ss->GetBinError(i));
      double err_ww  =  (hww_ss->GetBinError(i))*(hww_ss->GetBinError(i));
     
      double err_zmm  =  (hzmm_ss->GetBinError(i))*(hzmm_ss->GetBinError(i));
      double err_ztt  =  (hztt_ss->GetBinError(i))*(hztt_ss->GetBinError(i));

      double err_wjets  =  (hwjets_ss->GetBinError(i))*(hwjets_ss->GetBinError(i));
      double err_ttjets  =  (httjets_ss->GetBinError(i))*(httjets_ss->GetBinError(i));

      double err_data =(dataP_SS->GetBinError(i))*(dataP_SS->GetBinError(i));

      double err_QCD =TMath::Sqrt(err_zz+err_wz+err_ww+err_ztt+err_zmm+err_wjets+err_ttjets+err_data);*/

      if(hDataDataHistPSS->GetBinContent(i) < 0){

 	hDataDataHistPSS->SetBinContent(i,0);
        hDataDataHistPSS->SetBinError(i,0);

      }

  }

  hDataDataHistPSS->Write("QCD");

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

  std::map<std::string, TTree*> reducedTreesF;
  std::map<std::string, TTree*> reducedTreesSSF;
  std::map<std::string, TTree*> reducedTreesHiMtF;

  reducedTreesF.insert ( std::pair<std::string, TTree*>("data", fullTreeDataCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZCutF) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZCutF) );

  reducedTreesF.insert ( std::pair<std::string, TTree*>("zmmMatch", fullTreeSgnCutFTemp) );
  reducedTreesF.insert ( std::pair<std::string, TTree*>("zmmJets", fullTreeSgnCutFTempJet) );

  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("data", fullTreeDataSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZSSCutF) );
  reducedTreesSSF.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZSSCutF) );

  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("data", fullTreeDataHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZHiMtCutF) );
  reducedTreesHiMtF.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZHiMtCutF) );

  //////////////////////////////////////

  // string TreeNameF = "./Prova2506/histForSaving_"+tnp_+"_"+category_ + "_fail";
  TFile *SaveF = new TFile(Form("%s%s.root",TreeNameP.c_str(), binCenter_.c_str()),"UPDATE");

  SaveF->mkdir("fail");
  SaveF->cd("fail");

  TH1F* hDummy2 = new TH1F("hDummy","",50,70,120); 
  hDummy2->SetBinContent(1,1);
  hDummy2->Write("ggH");

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistF = histoProducer(reducedTreesF["data"], 1, 1, "mass", 70, 120, 50);
  TH1F * hDataDataHistFSS = histoProducer(reducedTreesSSF["data"], 1, 1, "mass", 70, 120, 50);
  TH1F * hDataDataHistFHiMt = histoProducer(reducedTreesHiMtF["data"], 1, 1, "mass", 70, 120, 50);
 
  hDataDataHistF->Write("data_obs");
  std::vector<TH1F*> dataFailShift = CreateESMassShape(fullTreeDataCutF,"data_obs",1.0);

  delete fullTreeDataCutF;
  delete fullTreeDataSSCutF;
  delete fullTreeDataHiMtCutF;

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHistF = histoProducer(reducedTreesF["zmm"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistFTemp = histoProducer(reducedTreesF["zmmMatch"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistFTempNo = histoProducer(reducedTreesF["zmmJets"], scaleFactors["zmm"], 1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistSSF = histoProducer(reducedTreesSSF["zmm"], scaleFactors["zmm"], -1, "mass", 70, 120, 50);
  TH1F * hsgnDataHistHiMtF = histoProducer(reducedTreesHiMtF["zmm"], scaleFactors["zmm"], -1, "mass", 70, 120, 50);

  hsgnDataHistFTemp->Write("ZL");
  hsgnDataHistF->Write("zmumuAll");
  hsgnDataHistFTempNo->Write("zmumuJets");

  std::vector<TH1F*> zmmFailShift = CreateESMassShape(fullTreeSgnCutFTemp,"ZL",ScaleFactorSgn);

  delete fullTreeSgnCutF;
  delete fullTreeSgnSSCutF;
  delete fullTreeSgnHiMtCutF;

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHistF = histoProducer(reducedTreesF["wjets"], scaleFactors["wjets"], 1, "mass", 70, 120, 50);
  TH1F * hwjetsDataHistSSF = histoProducer(reducedTreesSSF["wjets"], scaleFactors["wjets"], -1, "mass", 70, 120, 50);
  TH1F * hwjetsDataHistHiMtF = histoProducer(reducedTreesHiMtF["wjets"], scaleFactors["wjets"], 1, "mass", 70, 120, 50);

  hwjetsDataHistF->Write("wjets");

  std::vector<TH1F*> wjetsFailShift = CreateESMassShape(fullTreeWJetsCutF,"wjets",ScaleFactorWJets);

  delete fullTreeWJetsCutF;
  delete fullTreeWJetsSSCutF;
  delete fullTreeWJetsHiMtCutF;

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHistF = histoProducer(reducedTreesF["ztt"], scaleFactors["ztt"], 1, "mass", 70, 120, 50);
  TH1F * hzttDataHistSSF = histoProducer(reducedTreesSSF["ztt"], scaleFactors["ztt"], -1, "mass", 70, 120, 50);
  TH1F * hzttDataHistHiMtF = histoProducer(reducedTreesHiMtF["ztt"], scaleFactors["ztt"], -1, "mass", 70, 120, 50);

  hzttDataHistF->Write("ZTT");

  std::vector<TH1F*> zttFailShift = CreateESMassShape(fullTreeZttCutF,"ZTT",ScaleFactorZtt);

  delete fullTreeZttCutF;
  delete fullTreeZttSSCutF;
  delete fullTreeZttHiMtCutF;

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHistF = histoProducer(reducedTreesF["ttjets"], scaleFactors["ttjets"], 1, "mass", 70, 120, 50);
  TH1F * httjetsDataHistSSF = histoProducer(reducedTreesSSF["ttjets"], scaleFactors["ttjets"], -1, "mass", 70, 120, 50);
  TH1F * httjetsDataHistHiMtF = histoProducer(reducedTreesHiMtF["ttjets"], scaleFactors["ttjets"], -1, "mass", 70, 120, 50);

  httjetsDataHistF->Write("TT");

  std::vector<TH1F*> ttjetsFailShift = CreateESMassShape(fullTreeTTJetsCutF,"TT",ScaleFactorTTJets);

  delete fullTreeTTJetsCutF;
  delete fullTreeTTJetsSSCutF;
  delete fullTreeTTJetsHiMtCutF;

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHistF = histoProducer(reducedTreesF["ww"], scaleFactors["ww"], 1, "mass", 70, 120, 50);
  TH1F * hwwDataHistSSF = histoProducer(reducedTreesSSF["ww"], scaleFactors["ww"], -1, "mass", 70, 120, 50);
  TH1F * hwwDataHistHiMtF = histoProducer(reducedTreesHiMtF["ww"], scaleFactors["ww"], -1, "mass", 70, 120, 50);

  hwwDataHistF->Write("ww");

  std::vector<TH1F*> wwFailShift = CreateESMassShape(fullTreeWWCutF,"ww",ScaleFactorWW);
  
  delete fullTreeWWCutF;
  delete fullTreeWWSSCutF;
  delete fullTreeWWHiMtCutF;

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHistF = histoProducer(reducedTreesF["wz"], scaleFactors["wz"], 1, "mass", 70, 120, 50);
  TH1F * hwzDataHistSSF = histoProducer(reducedTreesSSF["wz"], scaleFactors["wz"], -1, "mass", 70, 120, 50);
  TH1F * hwzDataHistHiMtF = histoProducer(reducedTreesHiMtF["wz"], scaleFactors["wz"], -1, "mass", 70, 120, 50);

  hwzDataHistF->Write("wz");

  std::vector<TH1F*> wzFailShift = CreateESMassShape(fullTreeWZCutF,"wz",ScaleFactorWZ);

  delete fullTreeWZCutF;
  delete fullTreeWZSSCutF;
  delete fullTreeWZHiMtCutF;

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHistF = histoProducer(reducedTreesF["zz"], scaleFactors["zz"], 1, "mass", 70, 120, 50);
  TH1F * hzzDataHistSSF = histoProducer(reducedTreesSSF["zz"], scaleFactors["zz"], -1, "mass", 70, 120, 50);
  TH1F * hzzDataHistHiMtF = histoProducer(reducedTreesHiMtF["zz"], scaleFactors["zz"], -1, "mass", 70, 120, 50);

  hzzDataHistF->Write("zz");

  std::vector<TH1F*> zzFailShift = CreateESMassShape(fullTreeZZCutF,"zz",ScaleFactorZZ);


  hDataDataHistF->Sumw2(); 
  hDataDataHistFSS->Sumw2(); 
  hDataDataHistFHiMt->Sumw2(); 

  hsgnDataHistF->Sumw2(); 
  hsgnDataHistSSF->Sumw2(); 
  hsgnDataHistHiMtF->Sumw2(); 
  hsgnDataHistFTemp->Sumw2(); 
  hsgnDataHistFTempNo->Sumw2(); 

  hwjetsDataHistF->Sumw2(); 
  hwjetsDataHistSSF->Sumw2();
  hwjetsDataHistHiMtF->Sumw2();

  hzttDataHistF->Sumw2(); 
  hzttDataHistSSF->Sumw2();
  hzttDataHistHiMtF->Sumw2();

  httjetsDataHistF->Sumw2();
  httjetsDataHistSSF->Sumw2();
  httjetsDataHistHiMtF->Sumw2();

  hwwDataHistF->Sumw2();
  hwwDataHistSSF->Sumw2(); 
  hwwDataHistHiMtF->Sumw2(); 

  hwzDataHistF->Sumw2();
  hwzDataHistSSF->Sumw2(); 
  hwzDataHistHiMtF->Sumw2(); 

  hzzDataHistF->Sumw2(); 
  hzzDataHistSSF->Sumw2(); 
  hzzDataHistHiMtF->Sumw2(); 

  /////////////////////////////////////////////////////////////////////////////////////

  TH1F *h2VV = new TH1F("h2VV", "h2VV", 50, 70, 120);
  h2VV->Add(hzzDataHistF);
  h2VV->Add(hwzDataHistF);
  h2VV->Add(hzzDataHistF);
  h2VV->Write("VV");

  TH1F *hVVTES_UpF = new TH1F("hVVTES_UpF", "hVVTES_UpF", 50, 70, 120);
  hVVTES_UpF->Add(zzFailShift[0]);
  hVVTES_UpF->Add(wzFailShift[0]);
  hVVTES_UpF->Add(wwFailShift[0]);
  hVVTES_UpF->Write("VV_CMS_scale_t_mutau_8TeVUp");
  
  TH1F *hVVTES_DownF = new TH1F("hVVTES_DownF", "hVVTES_DownF", 50, 70, 120);
  hVVTES_DownF->Add(zzFailShift[1]);
  hVVTES_DownF->Add(wzFailShift[1]);
  hVVTES_DownF->Add(wwFailShift[1]);
  hVVTES_DownF->Write("VV_CMS_scale_t_mutau_8TeVDown");
  
  TH1F *hVVMMS_UpF = new TH1F("hVVMMS_UpF", "hVVMMS_UpF", 50, 70, 120);
  hVVMMS_UpF->Add(zzFailShift[2]);
  hVVMMS_UpF->Add(wzFailShift[2]);
  hVVMMS_UpF->Add(wwFailShift[2]);
  hVVMMS_UpF->Write("VV_CMS_scale_m_mutau_8TeVUp");

  TH1F *hVVMMS_DownF = new TH1F("hVVMMS_DownF", "hVVMMS_DownF", 50, 70, 120);
  hVVMMS_DownF->Add(zzFailShift[3]);
  hVVMMS_DownF->Add(wzFailShift[3]);
  hVVMMS_DownF->Add(wwFailShift[3]);
  hVVMMS_DownF->Write("VV_CMS_scale_m_mutau_8TeVDown");

  delete fullTreeZZCutF;
  delete fullTreeZZSSCutF;
  delete fullTreeZZHiMtCutF;

  //WJets from SideBand
  hsgnDataHistHiMtF->Sumw2();
  hDataDataHistP2->Sumw2();
  hsgnDataHistHiMtF->Multiply(hDataDataHistP2); //rescaling signal

  hDataDataHistFHiMt->Add(hsgnDataHistHiMtF);
  hDataDataHistFHiMt->Add(hzttDataHistHiMtF);
  hDataDataHistFHiMt->Add(httjetsDataHistHiMtF);
  hDataDataHistFHiMt->Add(hwwDataHistHiMtF);
  hDataDataHistFHiMt->Add(hwzDataHistHiMtF);
  hDataDataHistFHiMt->Add(hzzDataHistHiMtF);
  float expWF = hwjetsDataHistF->Integral(0,hwjetsDataHistF->GetSize()) / hwjetsDataHistHiMtF->Integral(0,hwjetsDataHistHiMtF->GetSize());
  cout<<"expWF "<<expWF<<endl;
  hDataDataHistFHiMt->Scale(expWF);
  hDataDataHistFHiMt->Write("wjetsDataDriven");

  float numWJetsDataDrivenF = hDataDataHistFHiMt->Integral(0,hDataDataHistFHiMt->GetSize());
  float normWJetsMCF = numWJetsDataDrivenF / hwjetsDataHistF->Integral(0,hwjetsDataHistF->GetSize());
  TH1F * wjetsNormF = (TH1F*)hwjetsDataHistF->Clone();
  wjetsNormF->Scale(normWJetsMCF);
  wjetsNormF->Write("W");

  //QCD datadriven
  hDataDataHistFSS->Add(hsgnDataHistSSF);
  hDataDataHistFSS->Add(hwjetsDataHistSSF);
  hDataDataHistFSS->Add(hzttDataHistSSF);
  hDataDataHistFSS->Add(httjetsDataHistSSF);
  hDataDataHistFSS->Add(hwwDataHistSSF);
  hDataDataHistFSS->Add(hwzDataHistSSF);
  hDataDataHistFSS->Add(hzzDataHistSSF);

  for(int i = 1; i < hDataDataHistFSS->GetSize(); i++){

      /*double err_zz  =  (hzz_ss->GetBinError(i))*(hzz_ss->GetBinError(i));
      double err_wz  =  (hwz_ss->GetBinError(i))*(hwz_ss->GetBinError(i));
      double err_ww  =  (hww_ss->GetBinError(i))*(hww_ss->GetBinError(i));
     
      double err_zmm  =  (hzmm_ss->GetBinError(i))*(hzmm_ss->GetBinError(i));
      double err_ztt  =  (hztt_ss->GetBinError(i))*(hztt_ss->GetBinError(i));

      double err_wjets  =  (hwjets_ss->GetBinError(i))*(hwjets_ss->GetBinError(i));
      double err_ttjets  =  (httjets_ss->GetBinError(i))*(httjets_ss->GetBinError(i));

      double err_data =(dataP_SS->GetBinError(i))*(dataP_SS->GetBinError(i));

      double err_QCD =TMath::Sqrt(err_zz+err_wz+err_ww+err_ztt+err_zmm+err_wjets+err_ttjets+err_data);*/

      if(hDataDataHistFSS->GetBinContent(i) < 0){

 	hDataDataHistFSS->SetBinContent(i,0);
        hDataDataHistFSS->SetBinError(i,0);

      }

  }

  hDataDataHistFSS->Write("QCD");

  SaveF->Close();

}

void calculateFitV1(){

	std::string looseV1 = "HpsAntiMuLoose";
	std::string mediumV1 = "HpsAntiMuMedium";
	std::string tightV1 = "HpsAntiMuTight";

  	cout<<"passingIsoLooseMuonVeto_LooseV1"<<endl;
	fitStudyTemplatesFromMC("muToTau",looseV1,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",looseV1,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",looseV1,"abseta>1.7",50,"Endcap");
  
  	cout<<"passingIsoLooseMuonVeto_MediumV1"<<endl;
	fitStudyTemplatesFromMC("muToTau",mediumV1,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",mediumV1,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",mediumV1,"abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_TightV1"<<endl;
  	fitStudyTemplatesFromMC("muToTau",tightV1,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",tightV1,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",tightV1,"abseta>1.7",50,"Endcap");

}

void calculateFitV2(){

	std::string looseV2 = "HpsAntiMuLoose2";
	std::string mediumV2 = "HpsAntiMuMedium2";
	std::string tightV2 = "HpsAntiMuTight2";

	/*cout<<"passingIsoLooseMuonVeto_LooseV2"<<endl;
	fitStudyTemplatesFromMC("muToTau",looseV2,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",looseV2,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",looseV2,"abseta>1.7",50,"Endcap");

       	cout<<"passingIsoLooseMuonVeto_MediumV2"<<endl;
       	fitStudyTemplatesFromMC("muToTau",mediumV2,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",mediumV2,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",mediumV2,"abseta>1.7",50,"Endcap");*/

	cout<<"passingIsoLooseMuonVeto_TightV2"<<endl;
	fitStudyTemplatesFromMC("muToTau",tightV2,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",tightV2,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",tightV2,"abseta>1.7",50,"Endcap");

}

void calculateFitV3(){

	std::string looseV3 = "HpsAntiMuLoose3";
	std::string tightV3 = "HpsAntiMuTight3";

	cout<<"passingIsoLooseMuonVeto_LooseV3"<<endl;
	fitStudyTemplatesFromMC("muToTau",looseV3,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",looseV3,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",looseV3,"abseta>1.7",50,"Endcap");
  
	cout<<"passingIsoLooseMuonVeto_TightV3"<<endl;
	fitStudyTemplatesFromMC("muToTau",tightV3,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",tightV3,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",tightV3,"abseta>1.7",50,"Endcap");

}

void calculateFitMVA(){

	std::string looseMVA = "HpsAntiMuLooseMVA";
	std::string mediumMVA = "HpsAntiMuMediumMVA";
	std::string tightMVA = "HpsAntiMuTightMVA";

    	cout<<"passingIsoLooseMuonVeto_LooseMVA"<<endl;
	fitStudyTemplatesFromMC("muToTau",looseMVA,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",looseMVA,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",looseMVA,"abseta>1.7",50,"Endcap");

	cout<<"passingIsoLooseMuonVeto_MediumMVA"<<endl;
	fitStudyTemplatesFromMC("muToTau",mediumMVA,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",mediumMVA,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",mediumMVA,"abseta>1.7",50,"Endcap");
	
	cout<<"passingIsoLooseMuonVeto_TightMVA"<<endl;
	fitStudyTemplatesFromMC("muToTau",tightMVA,"abseta<1.2",50,"Barrel");
	fitStudyTemplatesFromMC("muToTau",tightMVA,"abseta>1.2 && abseta<1.7",50,"Overlap");
	fitStudyTemplatesFromMC("muToTau",tightMVA,"abseta>1.7",50,"Endcap");

}

