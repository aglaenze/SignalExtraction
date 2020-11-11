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

#include "GetTemplates.C"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;


void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax) {
	
	int ptBinNumber = int(10*ptMax);
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooAbsPdf* fitModel = ws->pdf("model");

	
	// Define mass frame
	RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
	data->plotOn(mframe, Binning(50));
	
	
	// Define pt frame
	vector<RooPlot*> ptframes= {pt->frame(Title("Fit of p_{T} (log scale and full p_{T} range)")), pt->frame(Title("Fit of pt (zoom)"))};
	for (int k = 0; k<(int)ptframes.size(); k++) {
		data->plotOn(ptframes[k], Binning(ptBinNumber));
	}
	
	TCanvas* c1 = new TCanvas("2Dplot","2D fit",1800,1000) ;
	//TCanvas* c1 = new TCanvas("2Dplot","2D fit",800,300) ;
	c1->Divide(3,1) ; c1->SetCanvasSize(1350, 300);
	
	// Mass Plot
	c1->cd(1) ;
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	mframe->Draw();

	// pt plot
	c1->cd(2) ;
	gPad->SetLogy();
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	double yMax2 = ptframes[0]->GetMaximum();
	double y1 = 0.75*yMax2, y2 = 0.65*yMax2, y3 = 0.57*yMax2, y4 = 0.5*yMax2;
	double ptRangeMax=ptMax;
	ptframes[0]->GetXaxis()->SetRangeUser(0, ptRangeMax);
	ptframes[0]->Draw();
	ptframes[0]->Print();
	
	// Zoom pt
	c1->cd(3);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	ptframes[1]->GetXaxis()->SetRangeUser(0, ptRangeMax);
	ptframes[1]->Draw();

	
	// save plot
	string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	c1->SaveAs(Form("Plots/%s/NakedPlot%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax));
	
}

void NakedPlot(string rootfilePath, string rootfilePathMC, vector<string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2., Double_t mMax = 4.2, Double_t ptMin = 0., Double_t ptMax = 4, bool useCuts = true) {
	
	bool excOnly = false;
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.07);
	gStyle->SetTextSize(.07);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
	const int nPeriod = periods.size();
	
	for (int k = 0; k<nPeriod; k++) {
		string period = periods[k];
		
		if (excOnly) {
			if (period == "LHC16r") ptMax = 1.8;
			else ptMax = 0.8;
		}
		else {
			if (period == "LHC16s") ptMax = 3;
		}
		
		// Define cuts
		std::list<TCut> mCutList = DefineCuts(period, excOnly);
		TCut mCut = "";
		
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		if (!useCuts) mCut = "";
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
		
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax);

		wspace->Print();
		MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax);
		
	}
}


