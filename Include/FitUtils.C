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
#include "RooFormulaVar.h"
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

#include "ExtendedCrystalBall.h"
#include "Utils.C"


using namespace RooFit;
using namespace std;

//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMuPi", "kIncohPsi2sToMu", "kTwoGammaToMuLow", "kTwoGammaToMuMedium"};
//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMu"};

list<TCut> DefineCuts(string period, bool excOnly) {
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
	V0cut[0] = "fV0ADecision != 1";		// Pb does not dissociate
	V0cut[1] = "";
	
	ADcut[0] = "!(fADADecision==1) && fADABBNHits==0";	// Pb does not dissociate
	ADcut[1] = "!(fADCDecision==1) && fADCBBNHits==0";
	
	//bool exclusiveOnly = false;	// if true, p is not allowed to break
	bool inclusiveOnly = false;

	
	TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
	list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
	//list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ADcut[i], V0cut[i], "fV0CBBNHits < 3"};
	
	if (excOnly) {
		list<TCut> cutListExclusive = {"fTracklets == 0", "fV0CBBNHits < 3", "fV0CBGNHits == 0", "(fADCDecision==0) && (fADADecision==0)"};
		for (list<TCut>::iterator it=cutListExclusive.begin(); it!=cutListExclusive.end(); it++)
		mCutList.push_back(*it);
	}
	if (inclusiveOnly) {
		//list<TCut> cutListInclusive = {"fV0CBBNHits > 2"};
		list<TCut> cutListInclusive = {"fTracklets > 0"};
		for (list<TCut>::iterator it=cutListInclusive.begin(); it!=cutListInclusive.end(); it++)
		mCutList.push_back(*it);
	}

	//list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
	//list<TCut> mCutList = {"", "fAnaType==0", ZNcut[i], ADcut[i], V0cut[i]};
	//list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i]};
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

void ImportDataSet(RooWorkspace* ws, TTree* tree, Double_t mMin = 1.5, Double_t mMax = 5, Double_t ptMin = 0, Double_t ptMax = 8) {
	ImportDataSet(ws, tree, "fTrkTrkM>0", mMin, mMax, ptMin, ptMax);
	
}


//RooDataHist* GetPtHist(string period, TString process){
void GetPtHistMC(RooWorkspace* ws, string rootfilePathMC, string period, TString process, Double_t ptMin = 0, Double_t ptMax = 4., TCut cutMc = ""){
	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	Int_t nBins = int((ptMax-ptMin)*100);
	TH1F* histPt = new TH1F("hPt"+process, "hPt"+process, nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hPt"+process, cutMc);
	RooDataHist* dTemplatePt = new RooDataHist("ptHist"+process,"ptHist"+process, RooArgList(pt),histPt);
	RooHistPdf* ptPdf = new RooHistPdf("pt"+process, "pt"+process, pt, *dTemplatePt);
	ws->import(*ptPdf);
}

void GetMHistMC(RooWorkspace* ws, string rootfilePathMC, string period, TString process, Double_t mMin = 2.5, Double_t mMax = 3.5, TCut cutMc = ""){
	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	
	RooRealVar m = *ws->var("fTrkTrkM");
	
	Int_t nBins = int((mMax-mMin)*100);
	TH1F* histM = new TH1F("hM"+process, "hM"+process, nBins, mMin, mMax);
	fAnaTree->Draw("fTrkTrkM>>hM"+process, cutMc);
	RooDataHist* dTemplateM = new RooDataHist("mHist"+process,"mHist"+process, RooArgList(m),histM);
	RooHistPdf* mPdf = new RooHistPdf("m"+process, "m"+process, m, *dTemplateM);
	ws->import(*mPdf);
}


void GetV0Template(RooWorkspace* ws, string rootfilePath, string period) {
	
	// Open the file
	TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	
	// Connect to the tree
	TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
	
	TCut newCut = "fV0CBBNHits>7 && fTrkTrkM<3";
	if (period == "LHC16s") newCut = "fADABBNHits>1 && fTrkTrkM<3";
	Double_t ptMin = 0, ptMax = 8;
	
	Int_t nBins = int((ptMax-ptMin)*10);
	TH1F* hist = new TH1F("hV0C7", "hV0C7", nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hV0C7", newCut);
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataHist* dTemplatePt = new RooDataHist("hV0C7","hV0C7", RooArgList(pt),hist);
	RooHistPdf* ptPdf = new RooHistPdf("ptV0C7", "ptV0C7", pt, *dTemplatePt);
	
	ws->import(*ptPdf);
	//gStyle->SetOptStat(0);
}


void LoadMassFitFunctions(RooWorkspace* ws, string period) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	// First mass PDFs
	// J/Psi peak
	// Take tails parameters from TailParameters.C
	double alphaL = 0.961839, nL = 7.521515, alphaR = 2.641260, nR = 3.325886;
	if (period == "LHC16s") {alphaL = 0.993482; nL = 6.845735; alphaR = 2.669157; nR = 3.078395;}
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3,3.0,3.3);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.0811359, 0.05, 0.15);
	//RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.081);
	RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L", alphaL);
	RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L", nL);
	RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R", alphaR);
	RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R", nR);
	//RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
	ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", m,
														mean_jpsi, sigma_jpsi, alpha_jpsi_L,
														n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
	
	ws->import(*jpsi);
	
	// Then Psi(2s)
	// Most of the time not needed
	// That is kIncohPsi2sToMu in MC data
	double alphaL2 = 1.200001, nL2 = 3.017759, alphaR2 = 2.928444, nR2 = 2.256593;
	if (period == "LHC16s") {alphaL2 = 1.192599; nL2 = 3.118819; alphaR2 = 2.927051; nR2 = 2.055075;}
	
	Double_t factorMean = 1.1902; // maybe to adjust
								  //RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07,0,0.1);
								  // scaling factor Psi(2S) / JPsi
	Double_t factorSigma = 1.05; // maybe to adjust
	RooRealVar mean_scaling("mean_scaling", "", factorMean);
	RooRealVar sigma_scaling("sigma_scaling", "", factorSigma);
	RooFormulaVar mean_psi("mean_psi","mean_jpsi*mean_scaling",RooArgSet(mean_jpsi, mean_scaling));
	RooFormulaVar sigma_psi("sigma_psi", "sigma_jpsi*sigma_scaling", RooArgSet(sigma_jpsi, sigma_scaling));
	//RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07, 0, 0.15);
	RooRealVar alpha_psi_L("alpha_psi_L","alpha_psi_L",alphaL2);
	RooRealVar n_psi_L("n_psi_L","n_psi_L",nL2);
	RooRealVar alpha_psi_R("alpha_psi_R","alpha_psi_R",alphaR2);
	RooRealVar n_psi_R("n_psi_R","n_psi_R",nR2);
	//RooCBShape *psi = new RooCBShape("psi","crystal ball PDF",m,mean_psi,sigma_psi,alpha_psi,n_psi);
	ExtendedCrystalBall *psi2s = new ExtendedCrystalBall("psi2s","crystal ball PDF", m,
													   mean_psi, sigma_psi, alpha_psi_L,
													   n_psi_L, alpha_psi_R, n_psi_R);
	ws->import(*psi2s);
	
	// Finally background
	
	//Exponential background for mass
	RooRealVar a1("a1","a1",-1,-6,-0.5);
	RooExponential* bkg = new RooExponential("bkg","bkg",m,a1);
	/*
	 // using mass template from MC data
	 GetMHistMC(ws, rootfilePathMC, period, "kTwoGammaToMuLow", mMin, mMax);
	 RooAbsPdf* bkg = ws->pdf("mkTwoGammaToMuLow");
	 bkg->SetName("exp");
	 */
	ws->import(*bkg);
}

void LoadJpsiPtFitFunctions(RooWorkspace* ws, string rootfilePathMC, string period, bool exp, bool excOnly) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");

	// Contribution from exclusive J/Psis
	
	// Using H1 formula (b is free or not)
	double bExcValue = 4;
	double unused;
	if (! GetParameters(period, bExcValue, unused) ) {cout << "Parameters not loaded" << endl; return;}
	
	RooRealVar bExc("bExc","bExc", bExcValue, 3., 8);
	if (!excOnly) {bExc.setVal(bExcValue); bExc.setConstant();}
	
	// H1 formula
	RooGenericPdf *ptJpsiExclusive = new RooGenericPdf("ptJpsiExclusive","exclusive jPsi PDF","(2*fTrkTrkPt*exp(-bExc*(fTrkTrkPt**2)))",RooArgSet(pt,bExc)) ;
	
	/*
	 // using template pre-defined in sPlot
	 TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	 TH1F* hPtExclusive = (TH1F*)fTemplates->Get("hPtExclusive");
	 RooDataHist* ptHistExclusive = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hPtExclusive);
	 RooHistPdf* ptJpsiExclusive = new RooHistPdf("ptJpsiExclusive", "ptExclusive", pt, *ptHistExclusive);
	 */
	/*
	 // using pt template from MC data
	 GetPtHistMC(ws, rootfilePathMC, period, "kIncohJpsiToMu");
	 RooAbsPdf* ptJpsiExclusive = ws->pdf("ptkIncohJpsiToMu");
	 ptJpsiExclusive->SetName("ptJpsiExclusive");
	 */
	ws->import(*ptJpsiExclusive);
	
	// JPsi Dissociative
	/*
	 // using template pre-defined in sPlot
	 TH1F* hptJpsiDissociative = (TH1F*)fTemplates->Get("hptJpsiDissociative");
	 RooDataHist* ptHistDissociative = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hptJpsiDissociative);
	 RooHistPdf* ptJpsiDissociative = new RooHistPdf("ptJpsiDissociative", "ptJpsiDissociative", pt, *ptHistDissociative);
	 */
	RooGenericPdf *ptJpsiDissociative = nullptr;
	RooRealVar* bDiss = nullptr;
	RooRealVar* nDiss = nullptr;
	if (exp) {
		// H1 formula (the first one, with the exp)
		//RooRealVar bDiss("bDiss","bDiss", 0.323027, 0, 2);
		bDiss = new RooRealVar("bDiss","bDiss", 0.323027, 0.23, 1.8);
		ptJpsiDissociative = new RooGenericPdf("ptJpsiDissociative","Dissociative jPsi PDF","(2*fTrkTrkPt*exp(-bDiss*(fTrkTrkPt**2)))", RooArgSet(pt,*bDiss)) ;
	}
	else {
		// H1 formula (the second one, with the power law)
		bDiss = new RooRealVar("bDiss","bDiss", 1.6, 0.8, 4);
		nDiss = new RooRealVar("nDiss","nDiss", 3.6, 1, 10);
		ptJpsiDissociative = new RooGenericPdf("ptJpsiDissociative","Dissociative jPsi PDF","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bDiss/nDiss))**(-nDiss))", RooArgSet(pt, *nDiss, *bDiss)) ;
	}
	ws->import(*ptJpsiDissociative);
	
	// Contribution from gamma Pb events
	// using pt template from MC data
	GetPtHistMC(ws, rootfilePathMC, period, "kCohJpsiToMu");
	RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptkCohJpsiToMu");
	ptJpsiGammaPb->SetName("ptJpsiGammaPb");
	ws->import(*ptJpsiGammaPb);
	
	// Inclusive events
	/*
	 // From data with additional cuts
	 GetInclusiveTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
	 ptJpsiInclusive->SetName("ptJpsiInclusive");
	 */
	// With a new formula

	RooRealVar *pt0 = new RooRealVar("pt0","pt0", 4.5, 2.5, 10);
	RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.5, 2, 10);
/*
	 RooRealVar *pt0 = new RooRealVar("pt0","pt0", 4.5);
	 RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.5);
 */
	RooGenericPdf* ptJpsiInclusive = new RooGenericPdf("ptJpsiInclusive","Inclusive jPsi PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(pt, *pt0, *nInc)) ;
	
	ws->import(*ptJpsiInclusive);
	
	// Psi(2s)
	GetPtHistMC(ws, rootfilePathMC, period, "kIncohPsi2sToMu");
	RooAbsPdf* ptPsi2s = ws->pdf("ptkIncohPsi2sToMu");
	ptPsi2s->SetName("ptPsi2s");
	ws->import(*ptPsi2s);
}

void LoadBkgPtFitFunctions(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool excOnly) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");

	// Background
	// Using MC: gamma gamma to Mu Mu
	TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
	GetPtHistMC(ws, rootfilePathMC, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
	RooAbsPdf* ptTwoGamma = ws->pdf("ptkTwoGammaToMuLow");
	ptTwoGamma->SetName("ptTwoGamma");
	ws->import(*ptTwoGamma);
	
	// Extra background: pt distribution from background obtained with sPlot
	string suffix = "";
	if (excOnly) suffix = "-exclusive-only";
	TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s%s.root", rootfilePath.c_str(), period.c_str(), suffix.c_str()),"READ");
	//TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtExtraBackground__fTrkTrkPt");
	TH1F* hPtBackground = (TH1F*)fTemplates->Get("hSubSmooth");
	RooDataHist* ptHistBackground = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), hPtBackground);
	RooHistPdf* ptBackground = new RooHistPdf("ptBackground", "ptBackground", pt, *ptHistBackground);
	
	ws->import(*ptBackground);
	
	/*
	 // pt distribution from background is obtained with sPlot
	 TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	 //TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBackground");
	 TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBkgSmooth");
	 RooDataHist* ptHistBackground = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), hPtBackground);
	 RooHistPdf* ptBackground = new RooHistPdf("ptBackground", "ptBackground", pt, *ptHistBackground);
	 */
	/*
	 // using sidebands
	 GetSidebandsTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptBackground = ws->pdf("ptSidebands");
	 ptBackground->SetName("ptBackground");
	 */
	
}


void LoadYields(RooWorkspace* ws, string period, bool excOnly) {
	
	RooRealVar yieldJpsiExclusive("yieldJpsiExclusive","yieldJpsiExclusive",1200,400.,1700);
	RooRealVar yieldJpsiDissociative("yieldJpsiDissociative","yieldJpsiDissociative",1400,0.,1900);
	RooRealVar yieldJpsiGammaPb("yieldJpsiGammaPb","yieldJpsiGammaPb",30,0.,140);
	RooRealVar yieldJpsiInclusive("yieldJpsiInclusive","yieldJpsiInclusive",1000,0.,1600);
	
	RooRealVar yieldPsi2s("yieldPsi2s","yieldPsi2s",0,0.,1000);
	RooRealVar yieldTwoGamma("yieldTwoGamma","yieldTwoGamma",100,0.,240);
	RooRealVar yieldBkg("yieldBkg","yieldBkg",100,0.,800);
	if (period == "LHC16r") {
		yieldJpsiExclusive.setRange(900, 1400); yieldJpsiExclusive.setVal(1100);
		yieldJpsiInclusive.setRange(0,2000); yieldJpsiInclusive.setVal(1500);
	}
	else {
		yieldJpsiExclusive.setRange(400, 1100); yieldJpsiExclusive.setVal(900);
		yieldJpsiDissociative.setRange(0,400); yieldJpsiDissociative.setVal(200);
		//yieldJpsiInclusive.setRange(0,400);
		yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
	}
	 if (excOnly) {
		 yieldJpsiExclusive.setRange(500,1500); yieldJpsiExclusive.setVal(1000);
		 yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
		 //if (period == "LHC16s") {yieldJpsiDissociative.setVal(100);}
		 if (period == "LHC16s") {
			 yieldJpsiDissociative.setVal(0); yieldJpsiDissociative.setConstant();
		 }
		 else {
			 yieldJpsiDissociative.setRange(0,600); yieldJpsiDissociative.setVal(80);
		 }
	 }
	double yieldGammaPbValue = 40;
	double unused;
	if (! GetParameters(period, unused, yieldGammaPbValue) ) {cout << "Parameters not loaded" << endl; return;}
	yieldJpsiGammaPb.setVal(yieldGammaPbValue); yieldJpsiGammaPb.setConstant();
	
	ws->import(yieldJpsiExclusive);
	ws->import(yieldJpsiDissociative);
	ws->import(yieldJpsiGammaPb);
	ws->import(yieldJpsiInclusive);
	ws->import(yieldPsi2s);
	ws->import(yieldTwoGamma);
	ws->import(yieldBkg);
}

