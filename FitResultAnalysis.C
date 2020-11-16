#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <sstream>
#include <cmath>

#include <TROOT.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLine.h"


using namespace std;

//____________________________________________
TGraphErrors* CreateTGraph( Int_t size, const Double_t* x, const Double_t* y, const Double_t* xErr, const Double_t* yErr )
{
	TGraphErrors* tg = new TGraphErrors();
	for( Int_t i = 0; i < size; ++i )
	{
		tg->SetPoint( i, x[i], y[i] );
		tg->SetPointError( i, xErr[i], yErr[i] );
	}
	
	return tg;
}

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, const std::vector<Double_t>& xErr, const std::vector<Double_t>& yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr[0], &yErr[0] ); }

void MakeTree() {
	
	vector<string> periods = {"LHC16r", "LHC16s"};
	const int nPeriods = periods.size();
	
	Double_t a1, a1Err, bExc, bExcErr, bDiss, bDissErr, nDiss, nDissErr, nInc, nIncErr, pt0, pt0Err;
	Double_t jpsiExc, jpsiExcErr, jpsiDiss, jpsiDissErr, jpsiGammaPb, jpsiGammaPbErr, jpsiInc, jpsiIncErr, twoGamma, twoGammaErr, bkg, bkgErr, ratio, ratioErr;
	
	for (int k = 0; k<nPeriods; k++) {
		string period = periods[k];
		
		// Here build a tree
		TTree* t = new TTree("tFitResults", "Parameters resulting of the fit");
		t->Branch("a1", &a1);
		t->Branch("a1Err", &a1Err);
		t->Branch("bExc", &bExc);
		t->Branch("bExcErr", &bExcErr);
		t->Branch("bDiss", &bDiss);
		t->Branch("bDissErr", &bDissErr);
		t->Branch("nInc", &nInc);
		t->Branch("nIncErr", &nIncErr);
		t->Branch("pt0", &pt0);
		t->Branch("pt0Err", &pt0Err);
		t->Branch("jpsiExc", &jpsiExc);
		t->Branch("jpsiExcErr", &jpsiExcErr);
		t->Branch("jpsiDiss", &jpsiDiss);
		t->Branch("jpsiDissErr", &jpsiDissErr);
		t->Branch("jpsiGammaPb", &jpsiGammaPb);
		t->Branch("jpsiGammaPbErr", &jpsiGammaPbErr);
		t->Branch("jpsiInc", &jpsiInc);
		t->Branch("jpsiIncErr", &jpsiIncErr);
		t->Branch("twoGamma", &twoGamma);
		t->Branch("twoGammaErr", &twoGammaErr);
		t->Branch("bkg", &bkg);
		t->Branch("bkgErr", &bkgErr);
		t->Branch("ratio", &ratio);
		t->Branch("ratioErr", &ratioErr);
		
		ifstream file("output-" + period + ".txt", ios::in);
		if (file) {
			string line;
			while(getline(file, line)) {
				stringstream stream(line);
				stream >> a1 >> a1Err >> bExc >> bExcErr >> bDiss >> bDissErr >> nDiss >> nDissErr >> nInc >> nIncErr >> pt0 >> pt0Err;
				getline(file, line);
				stringstream stream2(line);
				stream2 >> jpsiExc >> jpsiExcErr >> jpsiDiss >> jpsiDissErr >> jpsiGammaPb >> jpsiGammaPbErr >> jpsiInc >> jpsiIncErr >> twoGamma >> twoGammaErr >> bkg >> bkgErr >> ratio >> ratioErr;
				t->Fill();
			}
			t->SaveAs(Form("FitResults-%s.root", period.c_str()));
		}
		
		else cout << "Error: not possible to open output-" << period << ".txt file in reading mode" << endl;
		
	}
	
}


void PlotResults() {
	
	vector<string> periods = {"LHC16r", "LHC16s"};
	const int nPeriods = periods.size();
	
	Double_t a1, a1Err, bExc, bExcErr, bDiss, bDissErr, nDiss, nDissErr, nInc, nIncErr, pt0, pt0Err;
	Double_t jpsiExc, jpsiExcErr, jpsiDiss, jpsiDissErr, jpsiGammaPb, jpsiGammaPbErr, jpsiInc, jpsiIncErr, twoGamma, twoGammaErr, bkg, bkgErr, ratio, ratioErr;
	
	vector<Double_t> a1Vec, a1ErrVec, bExcVec, bExcErrVec, bDissVec, bDissErrVec, nDissVec, nDissErrVec, nIncVec, nIncErrVec, pt0Vec, pt0ErrVec;
	vector<Double_t> jpsiExcVec, jpsiExcErrVec, jpsiDissVec, jpsiDissErrVec, jpsiGammaPbVec, jpsiGammaPbErrVec, jpsiIncVec, jpsiIncErrVec, twoGammaVec, twoGammaErrVec, bkgVec, bkgErrVec, ratioVec, ratioErrVec;
	
	for (int k = 0; k<nPeriods; k++) {
		string period = periods[k];
		TFile* f = new TFile(Form("FitResults-%s.root", period.c_str()), "READ");
		
		// Here build a tree
		TTree* t = (TTree*)f->Get("tFitResults");
		t->SetBranchAddress("a1", &a1);
		t->SetBranchAddress("a1Err", &a1Err);
		t->SetBranchAddress("bExc", &bExc);
		t->SetBranchAddress("bExcErr", &bExcErr);
		t->SetBranchAddress("bDiss", &bDiss);
		t->SetBranchAddress("bDissErr", &bDissErr);
		t->SetBranchAddress("nInc", &nInc);
		t->SetBranchAddress("nIncErr", &nIncErr);
		t->SetBranchAddress("pt0", &pt0);
		t->SetBranchAddress("pt0Err", &pt0Err);
		t->SetBranchAddress("jpsiExc", &jpsiExc);
		t->SetBranchAddress("jpsiExcErr", &jpsiExcErr);
		t->SetBranchAddress("jpsiDiss", &jpsiDiss);
		t->SetBranchAddress("jpsiDissErr", &jpsiDissErr);
		t->SetBranchAddress("jpsiGammaPb", &jpsiGammaPb);
		t->SetBranchAddress("jpsiGammaPbErr", &jpsiGammaPbErr);
		t->SetBranchAddress("jpsiInc", &jpsiInc);
		t->SetBranchAddress("jpsiIncErr", &jpsiIncErr);
		t->SetBranchAddress("twoGamma", &twoGamma);
		t->SetBranchAddress("twoGammaErr", &twoGammaErr);
		t->SetBranchAddress("bkg", &bkg);
		t->SetBranchAddress("bkgErr", &bkgErr);
		t->SetBranchAddress("ratio", &ratio);
		t->SetBranchAddress("ratioErr", &ratioErr);
		
		bool exp = false;
		
		const int n = t->GetEntries();
		for (int i = 0; i<n; i++) {
			t->GetEntry(i);
			a1Vec.push_back(a1);
			a1ErrVec.push_back(a1Err);
			bExcVec.push_back(bExc);
			bExcErrVec.push_back(bExcErr);
			bDissVec.push_back(bDiss);
			bDissErrVec.push_back(bDissErr);
			nDissVec.push_back(nDiss);
			nDissErrVec.push_back(nDissErr);
			nIncVec.push_back(nInc);
			nIncErrVec.push_back(nIncErr);
			pt0Vec.push_back(pt0);
			pt0ErrVec.push_back(pt0Err);
			jpsiExcVec.push_back(jpsiExc);
			jpsiExcErrVec.push_back(jpsiExcErr);
			jpsiDissVec.push_back(jpsiDiss);
			jpsiDissErrVec.push_back(jpsiDissErr);
			jpsiGammaPbVec.push_back(jpsiGammaPb);
			jpsiGammaPbErrVec.push_back(jpsiGammaPbErr);
			jpsiIncVec.push_back(jpsiInc);
			jpsiIncErrVec.push_back(jpsiIncErr);
			twoGammaVec.push_back(twoGamma);
			twoGammaErrVec.push_back(twoGammaErr);
			bkgVec.push_back(bkg);
			bkgErrVec.push_back(bkgErr);
			ratioVec.push_back(ratio);
			ratioErrVec.push_back(ratioErr);
			
			if (nDiss == 0) exp = true;
		}
		
		// Draw in a canvas
		// First parameters
		TCanvas* cv = new TCanvas();
		cv->Divide(2, 3);
		TGraphErrors* gr1 = CreateTGraph(bExcVec, a1Vec, bExcErrVec, a1ErrVec);
		gr1->SetTitle("a1 = f(b_{exc})");
		TGraphErrors* gr2 = CreateTGraph(bExcVec, bDissVec, bExcErrVec, bDissErrVec);
		gr2->SetTitle("b_{diss} = f(b_{exc})");
		TGraphErrors* gr3 = CreateTGraph(bExcVec, nDissVec, bExcErrVec, nDissErrVec);
		gr3->SetTitle("n_{diss} = f(b_{exc})");
		TGraphErrors* gr4 = CreateTGraph(bExcVec, pt0Vec, bExcErrVec, pt0ErrVec);
		gr4->SetTitle("p_{T0} = f(b_{exc})");
		TGraphErrors* gr5 = CreateTGraph(bExcVec, nIncVec, bExcErrVec, nIncErrVec);
		gr5->SetTitle("n_{inc} = f(b_{exc})");
		
		
		cv->cd(1);
		gr1->Draw("ALP");
		cv->cd(2);
		gr2->Draw("ALP");
		cv->cd(3);
		if (!exp) gr3->Draw("ALP");
		cv->cd(4);
		gr4->Draw("ALP");
		cv->cd(5);
		gr5->Draw("ALP");
		
		cv->SaveAs(Form("Plots/Param-%s.pdf", period.c_str()));
		
		//Then yields
		TCanvas* cv2 = new TCanvas();
		cv2->Divide(2, 4);
		TGraphErrors* grr1 = CreateTGraph(bExcVec, jpsiExcVec, bExcErrVec, jpsiExcErrVec);
		grr1->SetTitle("# exc J/#Psi = f(b_{exc})");
		grr1->GetYaxis()->SetRange(0,3000);
		TGraphErrors* grr2 = CreateTGraph(bExcVec, jpsiDissVec, bExcErrVec, jpsiDissErrVec);
		grr2->SetTitle("# diss J/#Psi = f(b_{exc})");
		TGraphErrors* grr3 = CreateTGraph(bExcVec, jpsiGammaPbVec, bExcErrVec, jpsiGammaPbErrVec);
		grr3->SetTitle("# #gamma-Pb J/#Psi = f(b_{exc})");
		TGraphErrors* grr4 = CreateTGraph(bExcVec, jpsiIncVec, bExcErrVec, jpsiIncErrVec);
		grr4->SetTitle("# inc J/#Psi = f(b_{exc})");
		TGraphErrors* grr5 = CreateTGraph(bExcVec, twoGammaVec, bExcErrVec, twoGammaVec);
		grr5->SetTitle("# #gamma#gamma = f(b_{exc})");
		TGraphErrors* grr6 = CreateTGraph(bExcVec, bkgVec, bExcErrVec, bkgErrVec);
		grr6->SetTitle("# bkg = f(b_{exc})");
		TGraphErrors* grr7 = CreateTGraph(bExcVec, ratioVec, bExcErrVec, ratioErrVec);
		grr7->SetTitle("N_{diss}/N_{exc} = f(b_{exc})");
		
		
		cv2->cd(1);
		grr1->Draw("ALP");
		cv2->cd(2);
		grr2->Draw("ALP");
		cv2->cd(3);
		grr3->Draw("ALP");
		cv2->cd(4);
		grr4->Draw("ALP");
		cv2->cd(5);
		grr5->Draw("ALP");
		cv2->cd(6);
		grr6->Draw("ALP");
		cv2->cd(7);
		grr7->Draw("ALP");
		
		cv2->SaveAs(Form("Plots/Yields-%s.pdf", period.c_str()));
		
		
	}
	
}
	
	
void FitResultAnalysis() {
	//MakeTree();
	PlotResults();
	
}
