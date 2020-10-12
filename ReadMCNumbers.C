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

int nEventsSignal[2] = {1000, 600};
int nEventsBkg[2] = {500, 100};


using namespace std;

void ReadMCNumbers() {
	
	vector<string> periods = {"LHC16r", "LHC16s"};
	const int nPeriods = periods.size();
	
	for (int k = 0; k<nPeriods; k++) {
		string period = periods[k];
		int nSig = nEventsSignal[k];
		int nBkg = nEventsBkg[k];
		TH1F* hSig1 = new TH1F("histSigSplot", "Number of J/Psi from sPlot", 50, nSig-100,  nSig+100);
		TH1F* hSig2 = new TH1F("histSigOriginal", "Number of J/Psi from 2D fit (original pt template)", 50, nSig-100,  nSig+100);
		TH1F* hSig3 = new TH1F("histSigTemp", "Number of J/Psi from sPlot (pt template from sPlot)", 50, nSig-100,  nSig+100);
		
		TH1F* hBkg1 = new TH1F("histBkgSplot", "Number of Bkg from sPlot", 50, nBkg-100,  nBkg+100);
		TH1F* hBkg2 = new TH1F("histBkgOriginal", "Number of Bkg from 2D fit (original pt template)", 50, nBkg-100,  nBkg+100);
		TH1F* hBkg3 = new TH1F("histBkgTemp", "Number of Bkg from sPlot (pt template from sPlot)", 50, nBkg-100,  nBkg+100);
		
		double nJpsiSplot, nJpsiOriginal, nJpsiTemp;
		double nBkgSplot, nBkgOriginal, nBkgTemp;
		
		ifstream monFluxSplot(Form("Candidate-Numbers/sPlot-Numbers-%s.txt", period.c_str()));  //Ouverture d'un fichier en lecture
		if(monFluxSplot)
		{
			string line; //Une variable pour stocker les lignes lues
			while(getline(monFluxSplot, line)) {
				stringstream stream(line);
				stream >> nJpsiSplot >> nBkgSplot;
				hSig1->Fill(nJpsiSplot);
				hBkg1->Fill(nBkgSplot);
			}
		}
		else {cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;}
		
		ifstream monFlux2dOriginal(Form("Candidate-Numbers/2Dfit-original-Numbers-%s.txt", period.c_str()));
		if(monFlux2dOriginal) {
			string line; //Une variable pour stocker les lignes lues
			
			while(getline(monFlux2dOriginal, line)) {
				stringstream stream(line);
				stream >> nJpsiOriginal >> nBkgOriginal;
				hSig2->Fill(nJpsiOriginal);
				hBkg2->Fill(nBkgOriginal);
			}
		}
		else {cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;}
		
		ifstream monFlux2dExtracted(Form("Candidate-Numbers/2Dfit-extracted-Numbers-%s.txt", period.c_str()));  //Ouverture d'un fichier en lecture
		if(monFlux2dExtracted) {
			string line; //Une variable pour stocker les lignes lues
			while(getline(monFlux2dExtracted, line)) //Tant qu'on n'est pas à la fin, on lit
			{
				stringstream stream(line);
				stream >> nJpsiTemp >> nBkgTemp;
				hSig3->Fill(nJpsiTemp);
				hBkg3->Fill(nBkgTemp);
			}
		}
		else {cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;}
		
		
		// Lines for original numbers
		double yMax1 = hSig1->GetMaximum()*1.05;
		TLine* lSig1 = new TLine(nSig, 0, nSig, yMax1);
		TLine* lBkg1 = new TLine(nBkg, 0, nBkg, yMax1);
		double yMax2 = hSig2->GetMaximum()*1.05;
		TLine* lSig2 = new TLine(nSig, 0, nSig, yMax2);
		TLine* lBkg2 = new TLine(nBkg, 0, nBkg, yMax2);
		double yMax3 = hSig3->GetMaximum()*1.05;
		TLine* lSig3 = new TLine(nSig, 0, nSig, yMax3);
		TLine* lBkg3 = new TLine(nBkg, 0, nBkg, yMax3);
		lSig1->SetLineColor(kRed);
		lSig2->SetLineColor(kRed);
		lSig3->SetLineColor(kRed);
		lBkg1->SetLineColor(kRed);
		lBkg2->SetLineColor(kRed);
		lBkg3->SetLineColor(kRed);
		
		// Original number of events
		TText* txtSig = new TText(nSig-85, yMax1*2./3, Form("%d input J/Psi", nSig));
		TText* txtBkg = new TText(nBkg-85, yMax1*2./3, Form("%d input Bkg", nBkg));
		
		TCanvas* cv = new TCanvas();
		cv->Divide(3, 2);
		// Signal on the top
		cv->cd(1);
		hSig1->Draw("hist");
		lSig1->Draw("same");
		txtSig->Draw("same");
		cv->cd(2);
		hSig2->Draw("hist");
		lSig2->Draw("same");
		cv->cd(3);
		hSig3->Draw("hist");
		lSig3->Draw("same");
		
		// Bkg on the bottom
		cv->cd(4);
		hBkg1->Draw("hist");
		lBkg1->Draw("same");
		txtBkg->Draw("same");
		cv->cd(5);
		hBkg2->Draw("hist");
		lBkg2->Draw("same");
		cv->cd(6);
		hBkg3->Draw("hist");
		lBkg3->Draw("same");
		
		cv->SaveAs(Form("Plots/CandidateFluctuationsMC-%s.pdf", period.c_str()));
	}
	
}
