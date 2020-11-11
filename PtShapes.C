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
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

using namespace std;


void PtShapes() {
	gStyle->SetOptStat(0);
	
	double ptMax = 4;
	int nBins = 400;
	double step = (double)ptMax/(nBins+1);
	
	// Get exclusive contribution
	double bH1 = 4;
	TH1F* histH1 = new TH1F("hexp", "exp function", nBins, 0, ptMax);
	
	// Get dissociative contribution
	double n = 3.2;
	double b = 1.2;
	TH1F* histDiss = new TH1F("hpowerlaw", "power law function", nBins, 0, ptMax);
	double x, y, yH1, y2;
	for (int k = 0; k < nBins; k++) {
		x = k*step;
		y = x*pow(1+pow(x,2) * (b/n),-n);
		yH1 = x*exp(-bH1*pow(x,2));
		histDiss->Fill(x,y);
		histH1->Fill(x,yH1);
	}
	histH1->Scale(1./histH1->GetMaximum());
	histDiss->Scale(1./histDiss->GetMaximum());
	
	// Get gamma-Pb
	TFile *fAnaSimu = new TFile("../p-Pb-2016/rootFiles/MC-std/AnalysisResults_LHC16r_MC_kCohJpsiToMu.root","READ");
	TTree* t = (TTree*)fAnaSimu->Get("MyTask/fAnaTree");
	TH1F* hist = new TH1F("hist", "hist", nBins, 0, ptMax);
	t->Draw("fTrkTrkPt>>hist");
	hist->Scale(1./hist->GetMaximum());
	TCanvas* cv = new TCanvas();
	histH1->SetTitle("p_{T} distributions");
	histH1->GetXaxis()->SetTitle("p_{T}");
	histH1->SetLineColor(kRed);
	histH1->Draw("hist");
	histDiss->SetLineColor(kGreen);
	histDiss->Draw("hist same");
	hist->SetLineColor(kBlue);
	hist->Draw("hist same");
	
	TLegend* lgd = new TLegend(0.4, 0.7, 0.9, 0.9);
	lgd->AddEntry(histH1, "Exclusive J/#Psi: p_{T} exp(-b p_{T}^{2})", "l");
	lgd->AddEntry(histDiss, "Inclusive J/#Psi: p_{T} (1-(b/n)*p_{T}^{2})^{-n}   ", "l");
	lgd->AddEntry(hist, "#gamma-Pb: template from MC", "l");
	lgd->Draw();
	cv->SaveAs("Plots/allpt.pdf");
	
}
