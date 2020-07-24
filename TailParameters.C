#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>
#include <TROOT.h>

#include "TSystem.h"
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
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "_FitUtils.C"
#include "ExtendedCrystalBall.h"

using namespace RooFit;

Double_t mMin = 2.5;
Double_t mMax = 4.0;

RooPlot* DrawCB(std::string rootfilePath, std::string period, TString process, bool log) {
	

	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePath.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fSimuTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	
	RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)",mMin, mMax);
	
	RooArgSet variables(m);

	// Import binned data (otherwise it never ends)
	TH1F* histM = new TH1F("hM"+process, "hM"+process, 300, 2.5, 4.1);
	fSimuTree->Draw("fTrkTrkM>>hM"+process);
	RooDataHist* dTemplateM = new RooDataHist("mHist"+process,"mHist"+process, RooArgList(m),histM);
	
	
	double x = 3., xmin = 2.9, xmax = 3.3;
	if (process == "kIncohPsi2sToMu") {x = 3.6, xmin = 3.5, xmax = 3.7;}
	
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",x,xmin,xmax);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.05,0,0.5);
	RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L",1, 0, 30);
	RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L",3, 0, 10);
	RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R",1, 0, 30);
	RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R",3, 0, 10);
	
	ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", m,
														mean_jpsi, sigma_jpsi, alpha_jpsi_L,
														n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
	RooRealVar yield1("yield"+process,"yield"+process,1.e6, 10, 1.e7);
	RooArgList* pdfList = new RooArgList(*jpsi);
	RooArgList yieldList = RooArgList(yield1);
	// Create fit model
	RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
	//fitModel->fitTo(*data,Minos(true),Strategy(2));
	RooFitResult* r = fitModel->fitTo(*dTemplateM ,Minos(true),Strategy(2), Save());
	Double_t alphaL = alpha_jpsi_L.getVal();
	Double_t nL = n_jpsi_L.getVal();
	Double_t alphaR = alpha_jpsi_R.getVal();
	Double_t nR = n_jpsi_R.getVal();
	
	std::cout << "alphaL = " << alphaL << std::endl;
	std::cout << "nL = " << nL << std::endl;
	std::cout << "alphaR = " << alphaR << std::endl;
	std::cout << "nR = " << nR << std::endl;
	
	RooPlot* mframe = m.frame(Title("Fit of invariant mass ("+process+")"));
	dTemplateM->plotOn(mframe);
	fitModel->plotOn(mframe, Name("sum"), LineColor(kRed));
	
	double yMax = mframe->GetMaximum();
	double y1 = 0.8*yMax, y2 = 0.7*yMax, y3 = 0.6*yMax, y4 = 0.5*yMax, y5 = 0.4*yMax;
	if (log) {y1 = yMax/pow(2,1); y2 = yMax/pow(2,2); y3 = yMax/pow(2,3); y4 = yMax/pow(2,4); y5 = yMax/pow(2,5);}
	double xPos = mMin + (mMax-mMin)*2./3;
	if (process == "kIncohPsi2sToMu") xPos = mMin + (mMax-mMin)*0.1;
	TLatex* txt1 = new TLatex(xPos,y1,Form("alpha_L = %f", alphaL));
	TLatex* txt2 = new TLatex(xPos,y2,Form("n_{L} = %f", nL));
	TLatex* txt3 = new TLatex(xPos,y3,Form("alpha_{R} = %f", alphaR));
	TLatex* txt4 = new TLatex(xPos,y4,Form("n_{R} = %f", nR));
	TLatex* txt5 = new TLatex(xPos,y5,Form("# candidates = %.1f", yield1.getVal()));
	mframe->addObject(txt1) ;
	mframe->addObject(txt2) ;
	mframe->addObject(txt3) ;
	mframe->addObject(txt4) ;
	mframe->addObject(txt5) ;
	
	return mframe;
}


void TailParameters(std::string rootfilePath="", std::vector<std::string> periods = {"LHC16r", "LHC16s"}, bool logScale = false) {
	
	gStyle->SetOptStat(0);
	
	gROOT->ProcessLine(".L ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./ExtendedCrystalBall_cxx.so") ;
	
	const int nPeriod = periods.size();
	for (int k = 0; k<nPeriod;k++) {
		std::string period = periods[k];
		
		RooPlot* mframeJpsi = DrawCB(rootfilePath, period, "kIncohJpsiToMu", logScale);
		RooPlot* mframePsi2s = DrawCB(rootfilePath, period, "kIncohPsi2sToMu", logScale);
		
		TCanvas* cv = new TCanvas("cv","cv",800,800) ;
		cv->Divide(1,2);
		
		cv->cd(1);
		if (logScale) gPad->SetLogy();
		gPad->SetLeftMargin(0.15) ;
		mframeJpsi->Draw();
		cv->cd(2);
		if (logScale) gPad->SetLogy();
		gPad->SetLeftMargin(0.15) ;
		mframePsi2s->Draw();
		
		cv->SaveAs(Form("Plots/2sidedCB-%s.pdf", period.c_str()));
	}
	
}
