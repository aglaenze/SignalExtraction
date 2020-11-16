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


void Background(string rootfilePath, string rootfilePathMc, vector<string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2., Double_t mMax = 4.2, Double_t ptMin = 0., Double_t ptMax = 4) {
	
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
		string period = periods[k];
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
		
		TCut cut = Form("fTrkTrkM < %f && fTrkTrkM > %f", mMax, mMin);
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)", mMin, mMax);
		RooRealVar pt("fTrkTrkPt","Dimuon p_{T} (GeV/c)", ptMin, ptMax);
		RooArgSet variables(m, pt);
		RooDataSet* data = new RooDataSet("data","data",variables,Import(*fAnaTree),Cut(cut));
		wspace->import(*data, Rename("data"));
		
		data->Print();
		wspace->Print();
		RooPlot* ptframe = pt.frame(Title("Fit of pt"));
		//data->plotOn(ptframe, DataError(RooAbsData::SumW2));
		data->plotOn(ptframe, Binning(50));
		Double_t ptRangeFitMax = 6;
		ptframe->GetXaxis()->SetRangeUser(0, ptRangeFitMax);
		
		TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
		GetPtHistMC(wspace, rootfilePathMc, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
		//GetV0Template(wspace, rootfilePath, period);
		
		RooAbsPdf* pdfMC = wspace->pdf("ptkTwoGammaToMuLow");
		pdfMC->SetName("ptMC");
		RooAbsPdf* pdfV0 = wspace->pdf("ptV0C7");
		
		RooRealVar* yieldMC = new RooRealVar("yieldMC","yieldMC",100,0,10000);
		RooRealVar* yieldV0 = new RooRealVar("yieldV0","yieldV0",100,0,10000);
		
		RooArgList* pdfList = new RooArgList(*pdfMC, *pdfV0);
		RooArgList yieldList = RooArgList(*yieldMC, *yieldV0);
		
		RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
		
		RooFitResult* r = fitModel->fitTo(*data, Extended(), Range(0,ptRangeFitMax, kTRUE), Minos(true), Strategy(1), Save());
		
		fitModel->plotOn(ptframe, Name("sum"), LineColor(kRed), LineWidth(1));
		fitModel->plotOn(ptframe,Name("ptMC"),Components(*pdfMC),LineStyle(kDashed), LineColor(2));
		//fitModel->plotOn(ptframe,Name("ptV0"),Components(*pdfV0),LineStyle(kDashed), LineColor(3));
		
		
		TCanvas* cv = new TCanvas("", "", 1200, 400);
		cv->Divide(3);
		cv->cd(2);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetLogy();
		ptframe->Draw();
		cv->cd(1);
		// Write number of gamma gamma candidates
		//Double_t xText = ptMin+(ptMax-ptMin)*2./3;
		Double_t xText = 1;
		Double_t yText = ptframe->GetMaximum()*2./3;
		TLatex* txtGg = new TLatex(xText, yText, Form("N_{#gamma#gamma} = %.2f #pm %.2f", yieldMC->getVal(), yieldMC->getError()));
		ptframe->addObject(txtGg) ;
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		ptframe->Draw();
		cv->cd(3);
		TLegend* lgd = new TLegend(0.1, 0.5, 0.9, 0.9);
		lgd->AddEntry(ptframe->findObject("ptMC"), "Monte Carlo", "L");
		lgd->AddEntry(ptframe->findObject("ptV0"), "Data with V0C cells > 7", "L");
		lgd->Draw();
		cv->SaveAs(Form("Plots/BackgroundTest-%s.pdf", period.c_str()));
		
	}
}


