#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include <TROOT.h>
#include <TMath.h>
#include "TString.h"
#include "TDatime.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include "_FitUtils.C"


using namespace RooFit;
using namespace std;

TH1F* smoothHist(TH1F* h, Double_t ptMin, Double_t ptMax) {
	int nBins = h->GetNbinsX();
	
	cout << "\n\n\n\n\n" << h->GetBinContent(1) << "\n\n\n\n\n" << endl;
	
	TH1F* sHist = new TH1F("hSmooth", "smooth pt hist", nBins, ptMin, ptMax);
	
	for (Int_t k = 0; k< nBins; k++) {
		int binNum = k+1;
		double size = binNum*0.8;	// number of bins to smooth
		if (size > 10) size = 10.;
		int val = h->GetBinContent(binNum);
		int nSmoothBins = 1;
		for (int j = 1; j < (int)size; j++ ) {
			if (binNum-j > ptMin) {
				nSmoothBins++;
				val += h->GetBinContent(binNum-j);
			}
			if (binNum+j < ptMax) {
				nSmoothBins++;
				val += h->GetBinContent(binNum+j);
			}
		}
		double mean = (double)val/nSmoothBins;
		sHist->SetBinContent(binNum, mean);
	}
	
	return sHist;
}

Double_t FitInclusive(Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*x[0]/pow(1+ (pow(x[0]/par[0],2)), par[1]); }

TF1* GetFitCurve(TH1F* h) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = 0;
	Int_t fitRangeMax = xMax + 5 * h->GetRMS();
	
	TF1* f = new TF1( "FitFunction", FitInclusive, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("pt0", "n", "amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}

void GetSidebandsTemplate(RooWorkspace* ws, std::string rootfilePath, std::string period,  Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax) {
	
	// Open the file
	TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	
	// Connect to the tree
	TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
	
	TCut newCut = "fTrkTrkM<2.5 && fTrkTrkM>1.5";
	
	Int_t nBins = int((ptMax-ptMin)*10);
	TH1F* hist = new TH1F("hSidebands", "hSidebands", nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hSidebands", newCut);
	//gStyle->SetOptStat(1111);
	TCanvas* cv = new TCanvas();
	hist->Draw();
	cv->SaveAs(Form("Plots/Sidebands-%s.pdf", period.c_str()));
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataHist* dTemplatePt = new RooDataHist("hSidebands","hSidebands", RooArgList(pt),hist);
	RooHistPdf* ptPdf = new RooHistPdf("ptSidebands", "ptSidebands", pt, *dTemplatePt);
	
	ws->import(*ptPdf);
	//gStyle->SetOptStat(0);
}

void GetInclusiveTemplate(RooWorkspace* ws, std::string rootfilePath, std::string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax) {
	
	// Open the file
	TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	
	// Connect to the tree
	TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
	
	//TCut newCut = mCut + Form("fTrkTrkM>%f && fTrkTrkM<%f && fTrkTrkPt>%f && fTrkTrkPt<%f", mMin, mMax, ptMin, ptMax);
	//TCut newCut = Form("fTrkTrkM>%f && fTrkTrkM<%f && (fTracklets>0 || fV0ARingMultiplicity[3]>10)", mMin, mMax);
	TCut newCut = Form("fTrkTrkM>%f && fTrkTrkM<%f && (fTracklets>1)", mMin, mMax);
	//TCut newCut = Form("fTrkTrkM>%f && fTrkTrkM<%f && (fTracklets>0)", mMin, mMax);
	//TCut newCut = Form("fTrkTrkM>%f && fTrkTrkM<%f && (fTracklets>1) && fV0CBBNHits>3", mMin, mMax);
	//newCut += "(fV0CCounts > 12 || fV0ACounts > 12) && fTrkTrkM>2.9 && fTrkTrkM<3.2";
	
	Int_t nBins = int((ptMax-ptMin)*10);
	TH1F* hist = new TH1F("hInclusive", "hInclusive", nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hInclusive", newCut);
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataHist* dTemplatePt = new RooDataHist("hInclusive","hInclusive", RooArgList(pt), hist);
	RooHistPdf* ptPdf = new RooHistPdf("ptInclusive", "ptInclusive", pt, *dTemplatePt);
	
	//gStyle->SetOptStat(1111);
	TCanvas* cv = new TCanvas();
	hist->SetLineColor(kBlack);
	hist->Draw();
	TH1F* sHist = smoothHist(hist, ptMin, ptMax);
	sHist->SetLineColor(kBlue);
	sHist->Draw("hist same");
	TF1* fitCurve = GetFitCurve(sHist);
	fitCurve->Draw("same");
	
	double pt0 = fitCurve->GetParameter(0);
	double pt0Error = fitCurve->GetParError(0);
	double slope = fitCurve->GetParameter(1);
	double slopeError = fitCurve->GetParError(1);
	double xPos = (ptMax-ptMin)*2./3;
	double yPos1 = hist->GetMaximum()*0.8;
	double yPos2 = hist->GetMaximum()*0.7;
	TLatex* t1 = new TLatex(xPos,yPos1, Form("p_{T0} = %.2f #pm %.2f", pt0, pt0Error));
	TLatex* t2 = new TLatex(xPos,yPos2, Form("n = %.2f #pm %.2f", slope, slopeError));
	t1->Draw();
	t2->Draw();
	cv->SaveAs(Form("Plots/ptTemplateInclusive-%s.pdf", period.c_str()));
	
	ws->import(*ptPdf);
	//gStyle->SetOptStat(0);
	
	/*
	 double pt0Val, nIncVal;
	 //if (period == "LHC16r") {nIncVal = 2.22749; pt0Val = 1.44421;}
	 //else {nIncVal = 2.97885; pt0Val = 0.813416;}
	 nIncVal = 2.22749; pt0Val = 1.44421;
	 RooRealVar *pt0 = new RooRealVar("pt0","pt0", 4, 0.7, 6);
	 //RooRealVar *pt0 = new RooRealVar("pt0","pt0", pt0Val);
	 RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.7, 2, 8);
	 //RooRealVar *nInc = new RooRealVar("nInc","nInc", nIncVal);
	 RooGenericPdf* ptInclusive = new RooGenericPdf("ptInclusive","Inclusive jPsi PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(pt, *pt0, *nInc)) ;
	 
	 RooRealVar yieldJpsiInclusive("yieldJpsiInclusive","yieldJpsiInclusive",100,1.,3.e3);
	 RooArgList yieldList = RooArgList(yieldJpsiInclusive);
	 RooArgList* pdfList = new RooArgList(*ptInclusive);
	 // Create fit model
	 RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
	 
	 RooDataSet* data = new RooDataSet("data2","data2",pt);
	 RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
	 
	 ws->import(*ptInclusive);
	 */
}



void GetTemplates(std::string rootfilePath, std::vector<std::string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2., Double_t mMax = 4.2, Double_t ptMin = 0., Double_t ptMax = 4) {
	
	//gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.05);
	gStyle->SetTextSize(.06);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
	const int nPeriod = periods.size();
	
	for (int k = 0; k<nPeriod; k++) {
		std::string period = periods[k];
		
		// Define cuts
		std::list<TCut> mCutList = DefineCuts(period);
		TCut mCut = "";
		
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
		
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax);
		
		GetInclusiveTemplate(wspace, rootfilePath, period, mMin, mMax, ptMin, ptMax);
		GetSidebandsTemplate(wspace, rootfilePath, period, mMin, mMax, ptMin, ptMax);
		
		
	}
}


