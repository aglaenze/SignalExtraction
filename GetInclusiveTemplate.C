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

Double_t FitInclusive(Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*x[0]/(1+ (x[0]/par[0])**2)**par[1]); }

TF1* GetFitCurve(TH1F* h) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = 0;
	Int_t fitRangeMax = xMax + 1.1 * h->GetRMS();
	
	TF1* f = new TF1( "FitFunction", FitInclusive, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("pt0", "n", "amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}

void GetInclusiveTemplate() {

	vector<string> periods = {"LHC16r", "LHC16s"};

	int nPeriods = (int)periods.size();
	for (int i = 0; i<nPeriods; i++) {
		string period = periods[i];
		TFile *f = new TFile(Form("../p-Pb-2016/rootFiles/std/AnalysisResults_%s.root", period.c_str()),"READ");
		TTree* t = (TTree*)f->Get("MyTask/fAnaTree");
		
		TCanvas* cv = new TCanvas();

		TH1F* hPt = new TH1F("hPtInclusive", "hPtInclusive", 200, 0, 20);
		t->Draw("fTrkTrkPt>>hPtInclusive", "fTracklets>1 && fTrkTrkM > 2.8 && fTrkTrkM < 3.3");
		hPt->Draw("hist");
		
		TF1* fitCurve = GetFitCurve(hPt);
		fitCurve->Draw("same");

		cv->SaveAs(Form("Plots/ptInclusive-%s.pdf", period.c_str()));
	}
}
