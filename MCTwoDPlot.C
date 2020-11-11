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

#include "ExtendedCrystalBall.h"

using namespace RooFit;
using namespace std;

Double_t mMin = 2.6, mMax = 3.5;
Double_t ptMin = 0., ptMax = 8;
int ptBinNumber = int(10*ptMax);

void ImportDataSet(RooWorkspace* ws, TTree* tree) {
	// import data
	
	RooRealVar m("m","M_{#mu#mu} (GeV/c2)", mMin, mMax);
	RooRealVar pt("pt","Dimuon p_{T} (GeV/c)", ptMin, ptMax);
	
	RooArgSet variables(m, pt);
	RooDataSet* data = new RooDataSet("data","data",variables,Import(*tree));
	
	data->Print();
	ws->import(*data, Rename("data"));
	
	Int_t nData = data->numEntries();
	Int_t nEntries = tree->GetEntries();
	cout << "\n\n Number of entries = " << nData << "\nAnd entries in TTree = "<< nEntries << "\n\n" << endl;
}

void AddModel(RooWorkspace* ws, string rootfilePathMC, string period, bool useOriginalPt) {
	// Define 2D model
	// First define fits in mass
	// Second define fits in pt
	// Then make the product to get 2D PDFs
	
	RooRealVar m = *ws->var("m");
	RooRealVar pt = *ws->var("pt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	// First mass PDFs
	// J/Psi peak
	// Take tails parameters from TailParameters.C
	double alphaL = 0.961839, nL = 7.521515, alphaR = 2.641260, nR = 3.325886;
	if (period == "LHC16s") {alphaL = 0.993482; nL = 6.845735; alphaR = 2.669157; nR = 3.078395;}
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3,2.9,3.3);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.0811359, 0.07, 1);
	//RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.081);
	RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L", alphaL);
	RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L", nL);
	RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R", alphaR);
	RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R", nR);
	//RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
	ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", m,
														mean_jpsi, sigma_jpsi, alpha_jpsi_L,
														n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
	// Finally background
	//Exponential background for mass
	//RooRealVar a1("a1","a1",-2,-10,-0.5);
	//RooRealVar a1("a1","a1",0);
	RooRealVar a1("a1","a1", 1, 0, 3);
	RooExponential bkg("exp","exp",m,a1);
	
	// Now pt PDFs
	
	// Contribution from exclusive J/Psis
	
	// Using H1 formula (b is free)
	RooRealVar bExc("bExc","bExc", 2.79116, 2, 8);
	// H1 formula
	RooGenericPdf *ptPdfExclusive = new RooGenericPdf("jpsiExc","exclusive jPsi PDF","(pt*exp(-bExc*(pt**2)))",RooArgSet(pt,bExc)) ;
	
	// Background
	// pt distribution from background is obtained with sPlot
	TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePathMC.c_str(), period.c_str()),"READ");
	//TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBackground__pt");
	TH1F* hPtBackground = nullptr;
	if (useOriginalPt) hPtBackground = (TH1F*)fTemplates->Get("hPtBkgOriginal");
	else hPtBackground = (TH1F*)fTemplates->Get("hPtBkgSmooth");
	RooDataHist* ptHistBackground = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), hPtBackground);
	RooHistPdf* ptBackground = new RooHistPdf("ptBackground", "ptBackground", pt, *ptHistBackground);
	//RooExponential* ptBackground = new RooExponential("exp","exp",pt,a1);
	
	// Last step:
	// Product fit(m) x fit(pt)
	RooProdPdf* pdfJpsiExclusive = new RooProdPdf("pdfJpsiExclusive","jpsi*ptkIncohJpsiToMu",RooArgList(*jpsi,*ptPdfExclusive));
	RooProdPdf* pdfBackground = new RooProdPdf("pdfBackground","bkg*ptBackground",RooArgList(bkg,*ptBackground));
	
	// All yields
	RooRealVar yieldJpsiExclusive("yieldJpsiExclusive","yieldJpsiExclusive",2000, 0.1, 1.e4);
	RooRealVar yieldBkg("yieldBkg","yieldBkg",500,0.,2.e4);
	
	// Assemble all components in sets
	
	RooArgList* pdfList = new RooArgList(*pdfJpsiExclusive, *pdfBackground);
	RooArgList yieldList = RooArgList(yieldJpsiExclusive, yieldBkg);
	
	// Create fit model
	RooAbsPdf* fitModel;
	fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
	
	ws->import(*fitModel);
}

void MakePlots(RooWorkspace* ws, string period, bool drawPulls, bool logScale, bool useOriginalPt, bool write) {
	
	string suffix = "";
	if (useOriginalPt) suffix = "original";
	else suffix = "extracted";
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("m");
	RooRealVar* pt = ws->var("pt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooAbsPdf* fitModel = ws->pdf("model");
	
	// Fit data
	RooFitResult* r = fitModel->fitTo(*data, Extended(),Minos(true),Strategy(1), Save());
	
	RooAbsPdf* pdfJpsiExclusive = ws->pdf("pdfJpsiExclusive");
	RooAbsPdf* pdfBackground = ws->pdf("pdfBackground");
	
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldBkg = ws->var("yieldBkg");
	
	RooRealVar* bExc = ws->var("bExc");
	
	
	if (!write) {
		TCanvas* c1 = new TCanvas("2Dplot","2D fit",1800,1200) ;
		//TCanvas* c1 = new TCanvas("2Dplot","2D fit",800,300) ;
		if (drawPulls) c1->Divide(3,3) ;
		else {c1->Divide(3,2) ; c1->SetCanvasSize(1800, 800);}
		
		// Define mass frame
		RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
		data->plotOn(mframe, Binning(50));
		fitModel->plotOn(mframe, Name("sum"), LineColor(kRed), LineWidth(1));
		fitModel->plotOn(mframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(2), LineWidth(1));
		
		fitModel->plotOn(mframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(7), LineWidth(1));
		
		// Define pt frame
		RooPlot* ptframe = pt->frame(Title("Fit of pt"));
		data->plotOn(ptframe, Binning(ptBinNumber));
		fitModel->plotOn(ptframe, Name("sum"), LineColor(kRed), LineWidth(1));
		fitModel->plotOn(ptframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(2), LineWidth(1));
		fitModel->plotOn(ptframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(7), LineWidth(1));
		
		// Mass Plot
		c1->cd(2) ;
		gPad->SetLeftMargin(0.15) ;
		gPad->SetBottomMargin(0.15) ;
		if (logScale) gPad->SetLogy() ;
		mframe->Draw();
		
		// pt plot
		c1->cd(3) ;
		gPad->SetLeftMargin(0.15) ;
		gPad->SetBottomMargin(0.15) ;
		if (logScale) gPad->SetLogy() ;
		double yMax2 = ptframe->GetMaximum();
		double y1 = 0.75*yMax2, y2 = 0.65*yMax2, y3 = 0.57*yMax2;
		if (logScale) {y1 = yMax2/pow(2.,1), y2 = yMax2/pow(2.,2), y3 = yMax2/pow(2.,2.8);}
		TLatex* txtExc = new TLatex(ptMin+(ptMax-ptMin)*1/2, y1,Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
		ptframe->addObject(txtExc) ;
		ptframe->GetXaxis()->SetRangeUser(0,2.5);
		ptframe->Draw();
		
		// Quality plots: (data-fit)/sigma
		c1->cd(5);
		gPad->SetLeftMargin(0.15) ;
		//gPad->SetTopMargin(0.15) ;
		data->plotOn(mframe, Binning(50));
		fitModel->plotOn(mframe, Binning(50));
		RooHist* hpullM = mframe->pullHist();
		hpullM->GetXaxis()->SetRangeUser(mMin, mMax);
		hpullM->SetTitle("(data - fit)/#sigma for m distribution");
		hpullM->Draw("");
		
		c1->cd(6);
		gPad->SetLeftMargin(0.15) ;
		data->plotOn(ptframe, Binning(ptBinNumber));
		fitModel->plotOn(ptframe, Binning(ptBinNumber));
		RooHist* hpullPt = ptframe->pullHist();
		hpullPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
		hpullPt->SetTitle("(data - fit)/#sigma for p_{T} distribution");
		hpullPt->GetXaxis()->SetRangeUser(0,2.5);
		hpullPt->Draw("");
		
		// Quality plots: data-fit
		if (drawPulls) {
			c1->cd(8);
			gPad->SetLeftMargin(0.15) ;
			RooHist* hresidM = mframe->residHist();
			hresidM->GetXaxis()->SetRangeUser(mMin, mMax);
			hresidM->SetTitle("Residuals (data - fit) for m distribution");
			hresidM->Draw("");
			
			c1->cd(9);
			gPad->SetLeftMargin(0.15) ;
			RooHist* hresidPt = ptframe->residHist();
			hresidPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
			hresidPt->SetTitle("Residuals (data - fit) for p_{T} distribution");
			hresidPt->Draw("");
		}
		
		// Legend in a subcanvas
		c1->cd(1);
		TLegend* legend = new TLegend(0.1, 0.3, 0.9, 0.9);
		legend->SetFillColor(kWhite);
		legend->SetLineColor(kWhite);
		//legend->AddEntry(ptframe->findObject("pdfkCohJpsiToMu"), "kCohJpsiToMu","L");
		legend->AddEntry(ptframe->findObject("pdfJpsiExclusive"), "Exclusive J/Psi","L");
		legend->AddEntry(ptframe->findObject("pdfBackground"), "#gamma#gamma #rightarrow #mu^{+} #mu^{-} + other background","L");
		legend->AddEntry(ptframe->findObject("sum"),"sum","L");
		legend->Draw();
		
		// Number of different contributions in a subcanvas
		c1->cd(4);
		
		// Write number of candidates
		//TLatex* txt1 = new TLatex(3.3,0.9*yMax,Form("Coherent Jpsi : %.1f", yieldCohJpsi.getVal()));
		TLatex* txt2 = new TLatex(0.2,0.9,Form("Exclusive J/#Psi : %.1f #pm %.1f", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
		
		TLatex* txt6 = new TLatex(0.2,0.8,Form("#gamma#gamma #rightarrow #mu^{+} #mu^{-} : %.1f #pm %.1f", yieldBkg->getVal(), yieldBkg->getError()));
		txt2->Draw(); txt6->Draw();
		c1->SaveAs(Form("Plots/MC-Fit-2D-%s-%s-%.1f-%.1f.pdf", suffix.c_str(), period.c_str(), mMin, mMax));
	}
	
	// Ecrire les valeurs dans un fichier
	if (write) {
		string const nomFichier(Form("Candidate-Numbers/2Dfit-%s-Numbers-%s.txt", suffix.c_str(), period.c_str()));
		ofstream monFlux(nomFichier.c_str(), ios::app);
		
		if(monFlux) {
			monFlux << yieldJpsiExclusive->getVal() << " ";
			monFlux << yieldBkg->getVal() << endl;
		}
		else { cout << "File not opened" << endl;}
	}
	
}

void MCTwoDPlot(string rootfilePathMC, bool logScale = false, bool drawPulls = false, bool useOriginalPt = true, bool write = false) {
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.05);
	gStyle->SetTextSize(.07);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
	gROOT->ProcessLine(".L ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./ExtendedCrystalBall_cxx.so") ;
	
	vector<string> periods = {"LHC16r", "LHC16s"};
	const int nPeriod = periods.size();
	
	for (int k = 0; k<nPeriod; k++) {
		string period = periods[k];
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"READ");
		
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("tMixed");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		ImportDataSet(wspace, fAnaTree);
		AddModel(wspace, rootfilePathMC, period, useOriginalPt);
		wspace->Print();
		MakePlots(wspace, period, drawPulls, logScale, useOriginalPt, write);
		
		
	}
}


