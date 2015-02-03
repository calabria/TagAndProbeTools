
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

// tdrGrid: Turns the grid lines on (true) or off (false)

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.03);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.035, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.5);
  tdrStyle->SetTitleYOffset(1.5);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.035, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

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

typedef std::map<std::string, TTree*> mymap;
typedef std::map<std::string, float> mymap2;
typedef std::map<std::string, TH1F*> mymap3;

void makeHistos(mymap3 map, std::string nameSave){

  TCanvas *c2 = new TCanvas("canvas","canvas",10,30,700,700);
  c2->SetGrid(0,0);
  c2->SetFillStyle(4000);
  c2->SetFillColor(10);
  c2->SetTicky();
  c2->SetObjectStat(0);

  setTDRStyle();

  THStack * hs = new THStack("VisMass","VisMass");
  THStack * hs2 = new THStack("VisMass","VisMass");  
  TLegend *leg2 = new TLegend(0.7,0.6,0.95,0.9);
  leg2->SetFillColor(kWhite);
  leg2->SetLineColor(kWhite);

  map["ztt"]->SetLineColor(1);
  map["ztt"]->SetFillColor(95);
  hs->Add(map["ztt"]);
  hs2->Add(map["ztt"]);
  map["ttjets"]->SetLineColor(1);
  map["ttjets"]->SetFillColor(4);
  hs->Add(map["ttjets"]);
  hs2->Add(map["ttjets"]);
  map["ww"]->SetLineColor(1);
  map["ww"]->SetFillColor(5);
  hs->Add(map["ww"]);
  hs2->Add(map["ww"]);
  map["wz"]->SetLineColor(1);
  map["wz"]->SetFillColor(50);
  hs->Add(map["wz"]);
  hs2->Add(map["wz"]);
  map["zz"]->SetLineColor(1);
  map["zz"]->SetFillColor(7);
  hs->Add(map["zz"]);
  hs2->Add(map["zz"]);

  map["wjets"]->SetLineColor(1);
  map["wjets"]->SetFillColor(8);
  map["wjetsDD"]->SetLineColor(1);
  map["wjetsDD"]->SetFillColor(8);
  if(nameSave.find("failing") != std::string::npos){

	  hs->Add(map["wjets"]);
	  hs2->Add(map["wjets"]);

  }
  else if(nameSave.find("Met") != std::string::npos){

	  hs->Add(map["wjets"]);
	  hs2->Add(map["wjets"]);

  }
  else if(nameSave.find("Mt") != std::string::npos){

	  hs->Add(map["wjets"]);
	  hs2->Add(map["wjets"]);

  }
  else{

	  hs->Add(map["wjetsDD"]);
	  hs2->Add(map["wjetsDD"]);

  }

  map["qcd"]->SetLineColor(1);
  map["qcd"]->SetFillColor(6);
  hs->Add(map["qcd"]);
  hs2->Add(map["qcd"]);
  map["zmmMatch"]->SetLineColor(1);
  map["zmmMatch"]->SetFillColor(2);
  hs->Add(map["zmmMatch"]);
  hs2->Add(map["zmmMatch"]);
  map["zmmJets"]->SetLineColor(1);
  map["zmmJets"]->SetFillColor(40);
  hs->Add(map["zmmJets"]);
  hs2->Add(map["zmmJets"]);

  //gPad->SetLogy();
  //hs->SetMinimum(1);

  TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0, 0.2, 1.0, 1.0, 10);
  TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0, 0.05, 1.0, 0.2, 10);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  pad1->SetTickx();
  pad1->SetTicky();

  TPaveText * cmsPreliminaryLabel = new TPaveText(0.13, 0.92, 0.99, 1.01, "NDC");
  cmsPreliminaryLabel->AddText("CMS Preliminary Simulation #int L = 19.9 fb^{-1}");
  cmsPreliminaryLabel->SetTextAlign(12);
  cmsPreliminaryLabel->SetTextFont(42);
  cmsPreliminaryLabel->SetTextSize(0.03);
  cmsPreliminaryLabel->SetFillStyle(0);
  cmsPreliminaryLabel->SetBorderSize(0);

  hs->Draw("HIST");
  cmsPreliminaryLabel->Draw();
  TH1F * stack1 = (TH1F*)hs->GetStack()->Last()->Clone();
  stack1->SetMarkerSize(0);
  stack1->SetFillColor(1);
  stack1->SetFillStyle(3013);
  stack1->SetLineWidth(1);
  stack1->Draw("E2SAME");
  hs->SetTitle("");
  hs->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  if(nameSave.find("nPV") != std::string::npos) hs->GetXaxis()->SetTitle("# PV");
  if(nameSave.find("muMVAMet") != std::string::npos) hs->GetXaxis()->SetTitle("MET [GeV]");
  if(nameSave.find("tag_Mt") != std::string::npos) hs->GetXaxis()->SetTitle("M_{T} [GeV/c^{2}]");
  if(nameSave.find("pt") != std::string::npos) hs->GetXaxis()->SetTitle("p^{Probe}_{T} [GeV/c]");
  if(nameSave.find("tag_pt") != std::string::npos) hs->GetXaxis()->SetTitle("p^{Tag}_{T} [GeV/c]");
  if(nameSave.find("eta") != std::string::npos) hs->GetXaxis()->SetTitle("#eta^{Probe}");
  if(nameSave.find("eta_pt") != std::string::npos) hs->GetXaxis()->SetTitle("#eta^{Tag}");
  hs->GetYaxis()->SetTitle("Entries");
  map["data"]->SetMarkerStyle(20);
  map["data"]->SetMarkerColor(1);
  map["data"]->SetMarkerSize(1.0);
  map["data"]->Draw("EPSAME");

  leg2->AddEntry(map["data"], "Data", "p");
  leg2->AddEntry(map["zmmMatch"], "Signal (#mu#rightarrow#tau)", "f");
  leg2->AddEntry(map["zmmJets"], "Signal (jet#rightarrow#tau)", "f");
  leg2->AddEntry(map["ztt"], "Z#rightarrow#tau#tau", "f");
  if(nameSave.find("failing") != std::string::npos) leg2->AddEntry(map["wjets"], "WJets", "f");
  else leg2->AddEntry(map["wjetsDD"], "WJets", "f");
  leg2->AddEntry(map["ttjets"], "TTJets", "f");
  leg2->AddEntry(map["qcd"], "QCD", "f");
  leg2->AddEntry(map["ww"], "WW", "f");
  leg2->AddEntry(map["wz"], "WZ", "f");
  leg2->AddEntry(map["zz"], "ZZ", "f");
  leg2->Draw();

  pad2->cd();

  pad2->SetTickx();
  pad2->SetTicky();

  TH1F * stack = (TH1F*)hs2->GetStack()->Last()->Clone();
  TH1F * hDataDataHistPClone = (TH1F*)map["data"]->Clone();
  stack->Sumw2();
  hDataDataHistPClone->Sumw2();
  hDataDataHistPClone->Divide(stack);
  hDataDataHistPClone->SetStats(kFALSE);
  hDataDataHistPClone->SetMarkerStyle(20);
  hDataDataHistPClone->SetMarkerColor(1);
  hDataDataHistPClone->SetMinimum(0.5);
  hDataDataHistPClone->SetMaximum(1.5);
  hDataDataHistPClone->GetXaxis()->SetTitle("m_{tp} [GeV/c^{2}]");
  if(nameSave.find("nPV") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("# PV");
  if(nameSave.find("muMVAMet") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("MET [GeV]");
  if(nameSave.find("tag_Mt") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("M_{T} [GeV/c^{2}]");
  if(nameSave.find("pt") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("p^{Probe}_{T} [GeV/c]");
  if(nameSave.find("tag_pt") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("p^{Tag}_{T} [GeV/c]");
  if(nameSave.find("eta") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("#eta^{Probe}");
  if(nameSave.find("eta_pt") != std::string::npos) hDataDataHistPClone->GetXaxis()->SetTitle("#eta^{Tag}");
  hDataDataHistPClone->GetXaxis()->SetTitleSize(0.08);
  hDataDataHistPClone->GetXaxis()->SetTitleOffset(1.5);
  hDataDataHistPClone->GetXaxis()->SetLabelSize(0.08);
  hDataDataHistPClone->GetYaxis()->SetTitle("obs/exp");
  hDataDataHistPClone->GetYaxis()->SetTitleSize(0.08);
  hDataDataHistPClone->GetYaxis()->SetTitleOffset(0.5);
  hDataDataHistPClone->GetYaxis()->SetLabelSize(0.08);
  //hDataDataHistPClone->Fit("pol0");
  hDataDataHistPClone->Draw("ep");

  c2->SaveAs((nameSave+".png").c_str());
  c2->SaveAs((nameSave+".pdf").c_str());

}

void checkMass(mymap map, mymap mapSS, mymap mapHiMt, mymap2 mapSF, std::string var, int min, int max, int num, std::string name){

  //////////////////////////////////////
  //    Data
  //////////////////////////////////////

  TH1F * hDataDataHistP = histoProducer(map["data"], 1, 1, var.c_str(), min, max, num);
  TH1F * hDataDataHistPSS = histoProducer(mapSS["data"], 1, 1, var.c_str(), min, max, num);
  TH1F * hDataDataHistPHiMt = histoProducer(mapHiMt["data"], 1, 1, var.c_str(), min, max, num);

  //////////////////////////////////////
  //    sgn (Zmumu)
  //////////////////////////////////////

  TH1F * hsgnDataHist = histoProducer(map["zmm"], mapSF["zmm"], 1, var.c_str(), min, max, num);
  TH1F * hsgnDataHistTemp = histoProducer(map["zmmMatch"], mapSF["zmm"], 1, var.c_str(), min, max, num);
  TH1F * hsgnDataHistTempNo = histoProducer(map["zmmJets"], mapSF["zmm"], 1, var.c_str(), min, max, num);
  TH1F * hsgnDataHistSS = histoProducer(mapSS["zmm"], mapSF["zmm"], -1, var.c_str(), min, max, num);
  TH1F * hsgnDataHistHiMt = histoProducer(mapHiMt["zmm"], mapSF["zmm"], -1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    Wjets
  /////////////////////////////////////////

  TH1F * hwjetsDataHist = histoProducer(map["wjets"], mapSF["wjets"], 1, var.c_str(), min, max, num);
  TH1F * hwjetsDataHistSS = histoProducer(mapSS["wjets"], mapSF["wjets"], -1, var.c_str(), min, max, num);
  TH1F * hwjetsDataHistHiMt = histoProducer(mapHiMt["wjets"], mapSF["wjets"], 1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    Ztautau
  /////////////////////////////////////////

  TH1F * hzttDataHist = histoProducer(map["ztt"], mapSF["ztt"], 1, var.c_str(), min, max, num);
  TH1F * hzttDataHistSS = histoProducer(mapSS["ztt"], mapSF["ztt"], -1, var.c_str(), min, max, num);
  TH1F * hzttDataHistHiMt = histoProducer(mapHiMt["ztt"], mapSF["ztt"], -1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    TTJets
  /////////////////////////////////////////

  TH1F * httjetsDataHist = histoProducer(map["ttjets"], mapSF["ttjets"], 1, var.c_str(), min, max, num);
  TH1F * httjetsDataHistSS = histoProducer(mapSS["ttjets"], mapSF["ttjets"], -1, var.c_str(), min, max, num);
  TH1F * httjetsDataHistHiMt = histoProducer(mapHiMt["ttjets"], mapSF["ttjets"], -1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    WW
  /////////////////////////////////////////

  TH1F * hwwDataHist = histoProducer(map["ww"], mapSF["ww"], 1, var.c_str(), min, max, num);
  TH1F * hwwDataHistSS = histoProducer(mapSS["ww"], mapSF["ww"], -1, var.c_str(), min, max, num);
  TH1F * hwwDataHistHiMt = histoProducer(mapHiMt["ww"], mapSF["ww"], -1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    WZ
  /////////////////////////////////////////

  TH1F * hwzDataHist = histoProducer(map["wz"], mapSF["wz"], 1, var.c_str(), min, max, num);
  TH1F * hwzDataHistSS = histoProducer(mapSS["wz"], mapSF["wz"], -1, var.c_str(), min, max, num);
  TH1F * hwzDataHistHiMt = histoProducer(mapHiMt["wz"], mapSF["wz"], -1, var.c_str(), min, max, num);

  /////////////////////////////////////////
  //    ZZ
  /////////////////////////////////////////

  TH1F * hzzDataHist = histoProducer(map["zz"], mapSF["zz"], 1, var.c_str(), min, max, num);
  TH1F * hzzDataHistSS = histoProducer(mapSS["zz"], mapSF["zz"], -1, var.c_str(), min, max, num);
  TH1F * hzzDataHistHiMt = histoProducer(mapHiMt["zz"], mapSF["zz"], -1, var.c_str(), min, max, num);

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


  //QCD datadriven
  hDataDataHistPSS->Add(hsgnDataHistSS);
  hDataDataHistPSS->Add(hwjetsDataHistSS);
  hDataDataHistPSS->Add(hzttDataHistSS);
  hDataDataHistPSS->Add(httjetsDataHistSS);
  hDataDataHistPSS->Add(hwwDataHistSS);
  hDataDataHistPSS->Add(hwzDataHistSS);
  hDataDataHistPSS->Add(hzzDataHistSS);

  //Wjets datadriven
  hDataDataHistPHiMt->Add(hsgnDataHistHiMt);
  hDataDataHistPHiMt->Add(hzttDataHistHiMt);
  hDataDataHistPHiMt->Add(httjetsDataHistHiMt);
  hDataDataHistPHiMt->Add(hwwDataHistHiMt);
  hDataDataHistPHiMt->Add(hwzDataHistHiMt);
  hDataDataHistPHiMt->Add(hzzDataHistHiMt);
  float expW = hwjetsDataHist->Integral(0,hwjetsDataHist->GetSize()) / hwjetsDataHistHiMt->Integral(0,hwjetsDataHistHiMt->GetSize());
  cout<<"expW "<<expW<<endl;
  hDataDataHistPHiMt->Scale(expW);

  std::map<std::string, TH1F*> histos;
  histos.insert ( std::pair<std::string, TH1F*>("data", hDataDataHistP) );
  histos.insert ( std::pair<std::string, TH1F*>("zmm", hsgnDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("zmmMatch", hsgnDataHistTemp) );
  histos.insert ( std::pair<std::string, TH1F*>("zmmJets", hsgnDataHistTempNo) );
  histos.insert ( std::pair<std::string, TH1F*>("ztt", hzttDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("wjets", hwjetsDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("wjetsDD", hDataDataHistPHiMt) );
  histos.insert ( std::pair<std::string, TH1F*>("ttjets", httjetsDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("ww", hwwDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("wz", hwzDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("zz", hzzDataHist) );
  histos.insert ( std::pair<std::string, TH1F*>("qcd", hDataDataHistPSS) );

  makeHistos(histos, var+"_"+name);

}

void makePlots(
	const string tnp_                = "muToTau",
	const string category_           = "HpsAntiMuLoose > 0.5",
	const string bin_                = "abseta<1.2"
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

  TFile *McP = new TFile("dummy.root","RECREATE");

  std::string additionalCut_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";
  std::string additionalCutSS_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_Mt < 40 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";
  std::string additionalCutHiMt_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";
  std::string additionalCutHiMtSS_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_Mt > 60 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";
  std::string additionalCutNoMt_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge == 0 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";
  std::string additionalCutNoMtSS_ = "(pt > 20 && tag_pt > 25 && tag_muPFIsolation < 0.1 && abseta < 2.3 && DecayMode > 0.5 && HpsLooseCombIsoDBCorr3Hits > 0.5 && pair_charge != 0 && tag_muTriggerMatching > 0.5 && tag_triggerBitSingleMu > 0.5)";

  //Low MT
  TTree* fullTreeDataCut = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCut = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWJetsCut = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZttCut = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeTTJetsCut = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWWCut = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeWZCut = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeZZCut = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );

  TTree* fullTreeSgnCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );
  TTree* fullTreeSgnCutTempNo = fullTreeSgn->CopyTree( Form("(!mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCut_.c_str()) );

  //SS
  TTree* fullTreeDataSSCut = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeSgnSSCut = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWJetsSSCut = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZttSSCut = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeTTJetsSSCut = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWWSSCut = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeWZSSCut = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );
  TTree* fullTreeZZSSCut = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutSS_.c_str()) );

  //High MT
  TTree* fullTreeDataHiMtCut = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeSgnHiMtCut = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWJetsHiMtCut = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZttHiMtCut = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeTTJetsHiMtCut = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWWHiMtCut = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeWZHiMtCut = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeZZHiMtCut = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );

  TTree* fullTreeSgnHiMtCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );
  TTree* fullTreeSgnHiMtCutTempNo = fullTreeSgn->CopyTree( Form("(!mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMt_.c_str()) );

  //High MT SS
  TTree* fullTreeDataHiMtCutSS = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeSgnHiMtCutSS = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeWJetsHiMtCutSS = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeZttHiMtCutSS = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeTTJetsHiMtCutSS = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeWWHiMtCutSS = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeWZHiMtCutSS = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );
  TTree* fullTreeZZHiMtCutSS = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutHiMtSS_.c_str()) );

  //No MT
  TTree* fullTreeDataNoMtCut = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeSgnNoMtCut = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeWJetsNoMtCut = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeZttNoMtCut = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeTTJetsNoMtCut = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeWWNoMtCut = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeWZNoMtCut = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeZZNoMtCut = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );

  TTree* fullTreeSgnNoMtCutTemp = fullTreeSgn->CopyTree( Form("(mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );
  TTree* fullTreeSgnNoMtCutTempNo = fullTreeSgn->CopyTree( Form("(!mcTrue && %s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMt_.c_str()) );

  //No MT SS
  TTree* fullTreeDataNoMtCutSS = fullTreeData->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeSgnNoMtCutSS = fullTreeSgn->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeWJetsNoMtCutSS = fullTreeWJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeZttNoMtCutSS = fullTreeZtt->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeTTJetsNoMtCutSS = fullTreeTTJets->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeWWNoMtCutSS = fullTreeWW->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeWZNoMtCutSS = fullTreeWZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );
  TTree* fullTreeZZNoMtCutSS = fullTreeZZ->CopyTree( Form("(%s && %s && %s)",category_.c_str(),bin_.c_str(),additionalCutNoMtSS_.c_str()) );

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

  std::map<std::string, TTree*> reducedTrees;
  std::map<std::string, TTree*> reducedTreesSS;
  std::map<std::string, TTree*> reducedTreesHiMt;
  std::map<std::string, TTree*> reducedTreesHiMtSS;
  std::map<std::string, TTree*> reducedTreesNoMt;
  std::map<std::string, TTree*> reducedTreesNoMtSS;

  reducedTrees.insert ( std::pair<std::string, TTree*>("data", fullTreeDataCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZCut) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zmmMatch", fullTreeSgnCutTemp) );
  reducedTrees.insert ( std::pair<std::string, TTree*>("zmmJets", fullTreeSgnCutTempNo) );

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
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("zmmMatch", fullTreeSgnHiMtCutTemp) );
  reducedTreesHiMt.insert ( std::pair<std::string, TTree*>("zmmJets", fullTreeSgnHiMtCutTempNo) );

  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("data", fullTreeDataHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZHiMtCutSS) );
  reducedTreesHiMtSS.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZHiMtCutSS) );

  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("data", fullTreeDataNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZNoMtCut) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("zmmMatch", fullTreeSgnNoMtCutTemp) );
  reducedTreesNoMt.insert ( std::pair<std::string, TTree*>("zmmJets", fullTreeSgnNoMtCutTempNo) );

  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("data", fullTreeDataNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("zmm", fullTreeSgnNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("ztt", fullTreeZttNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("wjets", fullTreeWJetsNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("ttjets", fullTreeTTJetsNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("ww", fullTreeWWNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("wz", fullTreeWZNoMtCutSS) );
  reducedTreesNoMtSS.insert ( std::pair<std::string, TTree*>("zz", fullTreeZZNoMtCutSS) );

  std::map<std::string, float> scaleFactors;

  scaleFactors.insert ( std::pair<std::string, float>("zmm", ScaleFactorSgn) );
  scaleFactors.insert ( std::pair<std::string, float>("ztt", ScaleFactorZtt) );
  scaleFactors.insert ( std::pair<std::string, float>("wjets", ScaleFactorWJets) );
  scaleFactors.insert ( std::pair<std::string, float>("ttjets", ScaleFactorTTJets) );
  scaleFactors.insert ( std::pair<std::string, float>("ww", ScaleFactorWW) );
  scaleFactors.insert ( std::pair<std::string, float>("wz", ScaleFactorWZ) );
  scaleFactors.insert ( std::pair<std::string, float>("zz", ScaleFactorZZ) );

  std::string antiMu;
  std::string passing;
  std::string region;

  if(category_.find(">") != std::string::npos) passing = "passing";
  else passing = "failing";

  if(category_.find("HpsAntiMuLoose") != std::string::npos) antiMu = "loose";
  if(category_.find("HpsAntiMuMedium") != std::string::npos) antiMu = "medium";
  if(category_.find("HpsAntiMuTight") != std::string::npos) antiMu = "tight";
  if(category_.find("HpsAntiMuLoose2") != std::string::npos) antiMu = "loose2";
  if(category_.find("HpsAntiMuMedium2") != std::string::npos) antiMu = "medium2";
  if(category_.find("HpsAntiMuTight2") != std::string::npos) antiMu = "tight2";
  if(category_.find("HpsAntiMuLoose3") != std::string::npos) antiMu = "loose3";
  if(category_.find("HpsAntiMuTight3") != std::string::npos) antiMu = "tight3";
  if(category_.find("HpsAntiMuLooseMVA") != std::string::npos) antiMu = "looseMVA";
  if(category_.find("HpsAntiMuMediumMVA") != std::string::npos) antiMu = "mediumMVA";
  if(category_.find("HpsAntiMuTightMVA") != std::string::npos) antiMu = "tightMVA";

  if(bin_ == "abseta<1.2") region = "barrel";
  else if(bin_ == "abseta>1.2 && abseta<1.7") region = "overlap";
  else if(bin_ == "abseta>1.7") region = "endcap";

  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "mass", 70, 120, 50, antiMu+"_"+passing+"_"+region);
  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "event_nPV", 0, 40, 40, antiMu+"_"+passing+"_"+region);
  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "tag_muMVAMet", 0, 100, 50, antiMu+"_"+passing+"_"+region);

  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "pt", 0, 100, 50, antiMu+"_"+passing+"_"+region);
  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "tag_pt", 0, 100, 50, antiMu+"_"+passing+"_"+region);

  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region);
  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "tag_eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region);
  checkMass(reducedTrees, reducedTreesSS, reducedTreesHiMt, scaleFactors, "tag_Mt", 0, 160, 32, antiMu+"_"+passing+"_"+region);

  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "tag_Mt", 0, 160, 32, antiMu+"_"+passing+"_"+region+"_noMt");
  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "mass", 70, 120, 50, antiMu+"_"+passing+"_"+region+"_noMt");
  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "event_nPV", 0, 40, 40, antiMu+"_"+passing+"_"+region+"_noMt");
  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "tag_muMVAMet", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_noMt");

  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "pt", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_noMt");
  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "tag_pt", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_noMt");

  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region+"_noMt");
  checkMass(reducedTreesNoMt, reducedTreesNoMtSS, reducedTreesHiMt, scaleFactors, "tag_eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region+"_noMt");

  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "mass", 70, 120, 50, antiMu+"_"+passing+"_"+region+"_HiMt");
  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "event_nPV", 0, 40, 40, antiMu+"_"+passing+"_"+region+"_HiMt");
  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "tag_muMVAMet", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_HiMt");

  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "pt", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_HiMt");
  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "tag_pt", 0, 100, 50, antiMu+"_"+passing+"_"+region+"_HiMt");

  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region+"_HiMt");
  checkMass(reducedTreesHiMt, reducedTreesHiMtSS, reducedTreesHiMt, scaleFactors, "tag_eta", -2.5, 2.5, 50, antiMu+"_"+passing+"_"+region+"_HiMt");

}

void calcutateMass(){

	std::string looseV1P = "HpsAntiMuLoose > 0.5";
	std::string mediumV1P = "HpsAntiMuMedium > 0.5";
	std::string tightV1P = "HpsAntiMuTight > 0.5";

	std::string looseV2P = "HpsAntiMuLoose2 > 0.5";
	std::string mediumV2P = "HpsAntiMuMedium2 > 0.5";
	std::string tightV2P = "HpsAntiMuTight2 > 0.5";

	std::string looseV3P = "HpsAntiMuLoose3 > 0.5";
	std::string tightV3P = "HpsAntiMuTight3 > 0.5";

	std::string looseMVAP = "HpsAntiMuLooseMVA > 0.5";
	std::string mediumMVAP = "HpsAntiMuMediumMVA > 0.5";
	std::string tightMVAP = "HpsAntiMuTightMVA > 0.5";

	std::string looseV1F = "HpsAntiMuLoose < 0.5";
	std::string mediumV1F = "HpsAntiMuMedium < 0.5";
	std::string tightV1F = "HpsAntiMuTight < 0.5";

	std::string looseV2F = "HpsAntiMuLoose2 < 0.5";
	std::string mediumV2F = "HpsAntiMuMedium2 < 0.5";
	std::string tightV2F = "HpsAntiMuTight2 < 0.5";

	std::string looseV3F = "HpsAntiMuLoose3 < 0.5";
	std::string tightV3F = "HpsAntiMuTight3 < 0.5";

	std::string looseMVAF = "HpsAntiMuLooseMVA < 0.5";
	std::string mediumMVAF = "HpsAntiMuMediumMVA < 0.5";
	std::string tightMVAF = "HpsAntiMuTightMVA < 0.5";

	/////////////////////////////////////////////

	/*cout<<"passingIsoLooseMuonVetoLoose"<<endl;
	makePlots("muToTau",looseV1P,"abseta<1.2");
	makePlots("muToTau",looseV1F,"abseta<1.2");
	
	makePlots("muToTau",looseV1P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV1F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV1P,"abseta>1.7");
	makePlots("muToTau",looseV1F,"abseta>1.7");

	cout<<"passingIsoLooseMuonVetoMedium"<<endl;
	makePlots("muToTau",mediumV1P,"abseta<1.2");
	makePlots("muToTau",mediumV1F,"abseta<1.2");
	
	makePlots("muToTau",mediumV1P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumV1F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumV1P,"abseta>1.7");
	makePlots("muToTau",mediumV1F,"abseta>1.7");
	
	cout<<"passingIsoLooseMuonVetoTight"<<endl;
	makePlots("muToTau",tightV1P,"abseta<1.2");
	makePlots("muToTau",tightV1F,"abseta<1.2");
	makePlots("muToTau",tightV1P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV1F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV1P,"abseta>1.7");
	makePlots("muToTau",tightV1F,"abseta>1.7");*/

	/////////////////////////////////////////////

	/*cout<<"passingIsoLooseMuonVetoLooseV2"<<endl;
	makePlots("muToTau",looseV2P,"abseta<1.2");
	makePlots("muToTau",looseV2F,"abseta<1.2");
	
	makePlots("muToTau",looseV2P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV2F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV2P,"abseta>1.7");
	makePlots("muToTau",looseV2F,"abseta>1.7");

	cout<<"passingIsoLooseMuonVetoMediumV2"<<endl;
	makePlots("muToTau",mediumV2P,"abseta<1.2");
	makePlots("muToTau",mediumV2F,"abseta<1.2");
	
	makePlots("muToTau",mediumV2P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumV2F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumV2P,"abseta>1.7");
	makePlots("muToTau",mediumV2F,"abseta>1.7");
	
	cout<<"passingIsoLooseMuonVetoTightV2"<<endl;
	makePlots("muToTau",tightV2P,"abseta<1.2");
	makePlots("muToTau",tightV2F,"abseta<1.2");
	makePlots("muToTau",tightV2P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV2F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV2P,"abseta>1.7");
	makePlots("muToTau",tightV2F,"abseta>1.7");*/

	/////////////////////////////////////////////

	/*cout<<"passingIsoLooseMuonVetoLooseV3"<<endl;
	makePlots("muToTau",looseV3P,"abseta<1.2");
	makePlots("muToTau",looseV3F,"abseta<1.2");
	
	makePlots("muToTau",looseV3P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV3F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseV3P,"abseta>1.7");
	makePlots("muToTau",looseV3F,"abseta>1.7");

	cout<<"passingIsoLooseMuonVetoTightV3"<<endl;
	makePlots("muToTau",tightV3P,"abseta<1.2");
	makePlots("muToTau",tightV3F,"abseta<1.2");
	makePlots("muToTau",tightV3P,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV3F,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightV3P,"abseta>1.7");
	makePlots("muToTau",tightV3F,"abseta>1.7");*/

	/////////////////////////////////////////////

	/*cout<<"passingIsoLooseMuonVetoLooseMVA"<<endl;
	makePlots("muToTau",looseMVAP,"abseta<1.2");
	makePlots("muToTau",looseMVAF,"abseta<1.2");
	
	makePlots("muToTau",looseMVAP,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseMVAF,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",looseMVAP,"abseta>1.7");
	makePlots("muToTau",looseMVAF,"abseta>1.7");

	cout<<"passingIsoLooseMuonVetoMediumMVA"<<endl;
	makePlots("muToTau",mediumMVAP,"abseta<1.2");
	makePlots("muToTau",mediumMVAF,"abseta<1.2");
	
	makePlots("muToTau",mediumMVAP,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumMVAF,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",mediumMVAP,"abseta>1.7");
	makePlots("muToTau",mediumMVAF,"abseta>1.7");
	*/
	cout<<"passingIsoLooseMuonVetoTightMVA"<<endl;
	makePlots("muToTau",tightMVAP,"abseta<1.2");
	makePlots("muToTau",tightMVAF,"abseta<1.2");
	/*makePlots("muToTau",tightMVAP,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightMVAF,"abseta>1.2 && abseta<1.7");
	makePlots("muToTau",tightMVAP,"abseta>1.7");
	makePlots("muToTau",tightMVAF,"abseta>1.7");*/

}

