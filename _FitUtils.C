#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include <TROOT.h>
#include <TMath.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"


using namespace RooFit;

//std::vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMuPi", "kIncohPsi2sToMu", "kTwoGammaToMuLow", "kTwoGammaToMuMedium"};
//std::vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMu"};

std::list<TCut> DefineCuts(std::string period) {
	TCut ZNcut[2], ADcut[2], V0cut[2];
	
	int i = -1;
	if (period == "LHC16r") i = 0;
	else if (period == "LHC16s") i = 1;
	if (i == -1) {std::cout << "Wrong period in DefineCuts" << std::endl; return {};}
	
	ZNcut[0] = "!fZNAgoodTiming";	// no event in beam-beam time window (because UPC only)
	ZNcut[1] = "!fZNCgoodTiming";
	//ZNcut = "fZNAfired==0";
	
	/*
	 ZNcut[0] = "fZNATDC[0]>-2 && fZNATDC[0]<2";
	 ZNcut[1] = "fZNCTDC[0]>-2 && fZNCTDC[0]<2";
	 
	 for (int j = 1; j<4; j++) {
	 ZNcut[0] = ZNcut[0] || Form("fZNATDC[%d]>-2 && fZNATDC[%d]<2", j, j);
	 ZNcut[1] = ZNcut[1] || Form("fZNCTDC[%d]>-2 && fZNCTDC[%d]<2", j, j);
	 }
	 ZNcut[0] = !ZNcut[0];
	 ZNcut[1] = !ZNcut[1];
	 */
	V0cut[0] = "fV0ADecision == 0";		// Pb does not dissociate
	V0cut[1] = "";
	
	ADcut[0] = "!(fADADecision==1) && fADABBNHits==0";	// Pb does not dissociate
	ADcut[1] = "!(fADCDecision==1) && fADCBBNHits==0";
	
	bool exclusiveOnly = false;	// if true, p is not allowed to break
	bool inclusiveOnly = false;

	
	TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
	std::list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i], "fTracklets == 0"};
	//std::list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ADcut[i], V0cut[i], "fV0CBBNHits < 3"};
	
	if (exclusiveOnly) {
		std::list<TCut> cutListExclusive = {"fTracklets == 0", "fV0CBBNHits < 3", "fV0CBGNHits == 0", "fTrkTrkPt<4", "(fADCDecision==0) && (fADADecision==0)"};
		for (std::list<TCut>::iterator it=cutListExclusive.begin(); it!=cutListExclusive.end(); it++)
		mCutList.push_back(*it);
	}
	if (inclusiveOnly) {
		//std::list<TCut> cutListInclusive = {"fV0CBBNHits > 2"};
		std::list<TCut> cutListInclusive = {"fTracklets > 0"};
		for (std::list<TCut>::iterator it=cutListInclusive.begin(); it!=cutListInclusive.end(); it++)
		mCutList.push_back(*it);
	}

	//std::list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
	//std::list<TCut> mCutList = {"", "fAnaType==0", ZNcut[i], ADcut[i], V0cut[i]};
	//std::list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i]};
	return mCutList;
}


void ImportDataSet(RooWorkspace* ws, TTree* tree, TCut cut = "", Double_t mMin = 1.5, Double_t mMax = 5, Double_t ptMin = 0, Double_t ptMax = 8) {
	// import data
	
	RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)", mMin, mMax);
	RooRealVar pt("fTrkTrkPt","Dimuon p_{T} (GeV/c)", ptMin, ptMax);
	
	RooArgSet variables(m, pt);
	
	//RooRealVar pt1("fTrkPt1","pt1",0,12);
	RooRealVar anaType("fAnaType","fAnaType",-2,4);
	RooRealVar fZNAgoodTiming("fZNAgoodTiming","fZNAgoodTiming",-1,2);
	RooRealVar fZNCgoodTiming("fZNCgoodTiming","fZNCgoodTiming",-1,2);
	
	RooRealVar fZNAEnergy("fZNAEnergy","fZNAEnergy",-1.e5,1.e5);
	RooRealVar fZNCEnergy("fZNCEnergy","fZNCEnergy",-1.e5,1.e5);
	RooRealVar fZPAEnergy("fZPAEnergy","fZPAEnergy",-1.e5,1.e5);
	RooRealVar fZPCEnergy("fZPCEnergy","fZPCEnergy",-1.e5,1.e5);
	
	// V0
	RooRealVar fV0ADecision("fV0ADecision","fV0ADecision",-1,5);
	RooRealVar fV0CDecision("fV0CDecision","fV0CDecision",-1,5);
	RooRealVar fV0ACounts("fV0ACounts","fV0ACounts",-1,1000);
	RooRealVar fV0CCounts("fV0CCounts","fV0CCounts",-1,1000);
	RooRealVar fV0CBBNHits("fV0CBBNHits","fV0CBBNHits", 0, 33);
	RooRealVar fV0CBGNHits("fV0CBGNHits","fV0CBGNHits", 0, 33);
	
	RooRealVar fADADecision("fADADecision","fADADecision",-1,5);
	RooRealVar fADCDecision("fADCDecision","fADCDecision",-1,5);
	
	RooRealVar fADABBNHits("fADABBNHits","fADABBNHits",-1,30);
	RooRealVar fADCBBNHits("fADCBBNHits","fADCBBNHits",-1,30);
	
	//RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-20,20);
	RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-4.0,-2.5);
	//RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-3.2,-2.7);
	
	RooRealVar fTracklets("fTracklets","fTracklets",-2.0,3.e5);
	
	RooRealVar fTrkQ1("fTrkQ1","fTrkQ1",-5,5);
	RooRealVar fTrkQ2("fTrkQ2","fTrkQ2",-5,5);
	
	variables.add(anaType);
	variables.add(fZNAgoodTiming);
	variables.add(fZNCgoodTiming);
	
	/*
	variables.add(fZNAEnergy);
	variables.add(fZNCEnergy);
	variables.add(fZPAEnergy);
	variables.add(fZPCEnergy);
	 */
	
	variables.add(fV0ADecision);
	variables.add(fV0CDecision);

	variables.add(fV0ACounts);
	variables.add(fV0CCounts);
	variables.add(fV0CBBNHits);
	variables.add(fV0CBGNHits);
	
	variables.add(fADADecision);
	variables.add(fADCDecision);
	
	variables.add(fADABBNHits);
	variables.add(fADCBBNHits);
	
	variables.add(fTrkTrkY);
	
	variables.add(fTrkQ1);
	variables.add(fTrkQ2);
	variables.add(fTracklets);
	
	RooDataSet* data = new RooDataSet("data","data",variables,Import(*tree),Cut(cut));
	
	data->Print();
	//return data;
	ws->import(*data, Rename("data"));
	
	Int_t nData = data->numEntries();
	Int_t nEntries = tree->GetEntries();
	std::cout << "\n\nNumber of entries = " << nData << "\nAnd entries in TTree = "<< nEntries << "\n\n" << std::endl;
}


//RooDataHist* GetPtHist(std::string period, TString process){
void GetPtHistMC(RooWorkspace* ws, std::string rootfilePathMC, std::string period, TString process, Double_t ptMin = 0, Double_t ptMax = 4.){
	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	fAnaTree->Draw("fTrkTrkPt");
	
	RooRealVar pt("fTrkTrkPt","Dimuon p_{T} (GeV/c)",ptMin,ptMax);
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(pt),Import(*fAnaTree));
	
	Int_t nBins = int((ptMax-ptMin)*100);
	TH1F* histPt = new TH1F("hPt"+process, "hPt"+process, nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hPt"+process);
	RooDataHist* dTemplatePt = new RooDataHist("ptHist"+process,"ptHist"+process, RooArgList(pt),histPt);
	RooHistPdf* ptPdf = new RooHistPdf("pt"+process, "pt"+process, pt, *dTemplatePt);
	ws->import(*ptPdf);
}

