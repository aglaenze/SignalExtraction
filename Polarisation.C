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

#include "_FitUtils.C"

using namespace std;


void Polarisation() {
	gStyle->SetOptStat(0);
	
	Double_t mMin = 2.6;
	Double_t mMax = 3.5;
	vector<string> periods = {"LHC16r", "LHC16s"};
	
	for (int i = 0; i<(int)periods.size(); i++) {
		string period = periods[i];
		
		std::list<TCut> mCutList = DefineCuts(period);
		TCut mCut = "";
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		//mCut += "fTrkTrkM>1.7 && fTrkTrkM<2.2";
		mCut += Form("fTrkTrkM>%f && fTrkTrkM<%f", mMin, mMax);
		mCut += "fTrkTrkPt>0 && fTrkTrkPt<8";
		
		// Get data file
		TFile *fAna = new TFile(Form("../p-Pb-2016/rootFiles/std/AnalysisResults_%s.root", period.c_str()),"READ");
		TTree* t = (TTree*)fAna->Get("MyTask/fAnaTree");
		
		TCanvas* cv = new TCanvas();
		cv->Divide(3,2);
		cv->cd(1);
		t->Draw("fTrkTrkPt:fTrkTrkPhi", mCut, "colz");
		cv->cd(2);
		t->Draw("fTrkTrkPt:fTrkPhi1", mCut, "colz");
		cv->cd(3);
		t->Draw("fTrkTrkPt:fTrkPhi2", mCut, "colz");
		
		cv->cd(4);
		t->Draw("fTrkTrkPt:fTrkTrkY", mCut, "colz");
		cv->cd(5);
		t->Draw("fTrkTrkPt:fTrkEta1", mCut, "colz");
		cv->cd(6);
		t->Draw("fTrkTrkPt:fTrkEta2", mCut, "colz");
		
		cv->SaveAs(Form("Plots/polarisation-%s.pdf", period.c_str()));
	}
	
}
