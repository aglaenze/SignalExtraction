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

#include "Include/FitUtils.C"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;

void WriteResults(RooWorkspace* ws, string period, bool exp) {
	
	RooRealVar* a1 = ws->var("a1");
	RooRealVar* bExc = ws->var("bExc");
	RooRealVar* bDiss = ws->var("bDiss");
	RooRealVar* nDiss = nullptr;
	if (!exp) nDiss = ws->var("nDiss");
	else {nDiss->setVal(0); nDiss->setConstant();}
	
	RooRealVar* nInc = nullptr;
	RooRealVar* pt0 = nullptr;
	if (period == "LHC16r") {
		nInc = ws->var("nInc");
		pt0 = ws->var("pt0");
	}
	else {
		nInc->setVal(0); nInc->setConstant();
		pt0->setVal(0); pt0->setConstant();
	}
	
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
	RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
	//RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
	RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
	RooRealVar* yieldBkg = ws->var("yieldBkg");
	RooRealVar* r_diss_exc = ws->var("r_diss_exc");
	
	// Write the values in a text file
	string const file("output-"+period+".txt", ios::app);	// append to the file
	ofstream monFlux(file.c_str());
	
	if(monFlux) {
		// First parameters
		monFlux << a1->getVal() <<  " " << a1->getError() <<  " " << bExc->getVal() <<  " " << bExc->getError() <<  " " << bDiss->getVal() <<  " " << bDiss->getError() <<  " " << nDiss->getVal() <<  " " << nDiss->getError() <<  " " << nInc->getVal() <<  " " << nInc->getError() <<  " " << pt0->getVal() <<  " " << pt0->getError() << endl;
		// Then yields
		monFlux << yieldJpsiExclusive->getVal() <<  " " << yieldJpsiExclusive->getError() <<  " " << yieldJpsiDissociative->getVal() <<  " " << yieldJpsiDissociative->getError() <<  " " << yieldJpsiGammaPb->getVal() <<  " " << yieldJpsiGammaPb->getError() <<  " " << yieldJpsiInclusive->getVal() <<  " " << yieldJpsiInclusive->getError() <<  " " << yieldTwoGamma->getVal() <<  " " << yieldTwoGamma->getError() <<  " " << yieldBkg->getVal() <<  " " << yieldBkg->getError() <<  " " << r_diss_exc->getVal() <<  " " << r_diss_exc->getError() << endl;
	}
	else {
		cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
	}
}


void AddModel(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool exp, bool exclusiveOnly) {
	// Define 2D model
	// First define fits in mass
	// Second define fits in pt
	// Then make the product to get 2D PDFs
	
	LoadMassFitFunctions(ws, period);
	LoadJpsiPtFitFunctions(ws, rootfilePathMC, period, exp, exclusiveOnly);
	LoadBkgPtFitFunctions(ws, rootfilePath, rootfilePathMC, period, mMin, mMax, ptMin, ptMax, exclusiveOnly);
	// Now collect m pdf
	RooAbsPdf* jpsi = ws->pdf("jpsi");
	RooAbsPdf* psi2s = ws->pdf("psi2s");
	RooAbsPdf* bkg = ws->pdf("bkg");
	
	// Then pt pdf
	RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
	RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
	RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
	RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
	RooAbsPdf* ptPsi2s = ws->pdf("ptPsi2s");
	RooAbsPdf* ptTwoGamma = ws->pdf("ptTwoGamma");
	RooAbsPdf* ptBackground = ws->pdf("ptBackground");
	
	// Last step:
	// Product fit(m) x fit(pt)
	RooProdPdf* pdfJpsiExclusive = new RooProdPdf("pdfJpsiExclusive","jpsi*ptJpsiExclusive",RooArgList(*jpsi,*ptJpsiExclusive));
	RooProdPdf* pdfJpsiDissociative = new RooProdPdf("pdfJpsiDissociative","jpsi*ptJpsiDissociative",RooArgList(*jpsi,*ptJpsiDissociative));
	RooProdPdf* pdfJpsiGammaPb = new RooProdPdf("pdfJpsiGammaPb","jpsi*ptJpsiGammaPb",RooArgList(*jpsi,*ptJpsiGammaPb));
	RooProdPdf* pdfJpsiInclusive = new RooProdPdf("pdfJpsiInclusive","jpsi*ptJpsiInclusive",RooArgList(*jpsi,*ptJpsiInclusive));
	RooProdPdf* pdfPsi2s = new RooProdPdf("pdfPsi2s","psi2s*ptPsi2s",RooArgList(*psi2s,*ptPsi2s));
	RooProdPdf* pdfTwoGamma = new RooProdPdf("pdfTwoGamma","bkg*ptTwoGamma",RooArgList(*bkg,*ptTwoGamma));
	RooProdPdf* pdfBackground = new RooProdPdf("pdfBackground","bkg*ptBackground",RooArgList(*bkg,*ptBackground));
	
	// Load yields
	LoadYields(ws, period, exclusiveOnly);
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
	RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
	RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
	RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
	RooRealVar* yieldBkg = ws->var("yieldBkg");
	
	// Assemble all components in sets
	
	RooArgList* pdfList = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfPsi2s, *pdfTwoGamma, *pdfBackground, *pdfJpsiInclusive);
	RooArgList yieldList = RooArgList(*yieldJpsiExclusive, *yieldJpsiDissociative, *yieldJpsiGammaPb, *yieldPsi2s, *yieldTwoGamma, *yieldBkg, *yieldJpsiInclusive);
	RooArgList* pdfList2 = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfJpsiInclusive, *pdfTwoGamma, *pdfBackground);
	RooArgList yieldList2 = RooArgList(*yieldJpsiExclusive, *yieldJpsiDissociative, *yieldJpsiGammaPb, *yieldJpsiInclusive, *yieldTwoGamma, *yieldBkg);
	// Create fit model
	RooAbsPdf* fitModel;
	if (mMax > mLimitPsi2s) fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
	else fitModel = new RooAddPdf("model", "model", *pdfList2, yieldList2, kFALSE);
	
	ws->import(*fitModel);
}

void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly) {
	
	int ptBinNumber = int(10*ptMax);
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooAbsPdf* fitModel = ws->pdf("model");
	
	// Fit data
	RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
	//RooFitResult* r = fitModel->fitTo(*data, Minos(false), Strategy(0), Save());	// be quick (for testing)
	//RooFitResult* r = nullptr;
	
	RooAbsPdf* pdfJpsiExclusive = ws->pdf("pdfJpsiExclusive");
	RooAbsPdf* pdfJpsiDissociative = ws->pdf("pdfJpsiDissociative");
	RooAbsPdf* pdfJpsiGammaPb = ws->pdf("pdfJpsiGammaPb");
	RooAbsPdf* pdfTwoGamma = ws->pdf("pdfTwoGamma");
	RooAbsPdf* pdfBackground = ws->pdf("pdfBackground");
	RooAbsPdf* pdfJpsiInclusive = ws->pdf("pdfJpsiInclusive");
	
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
	RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
	RooRealVar* yieldBkg = ws->var("yieldBkg");
	RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
	
	RooRealVar* yieldPsi2s = nullptr;
	RooAbsPdf* pdfPsi2s = nullptr;
	if (mMax > mLimitPsi2s) {
		pdfPsi2s = ws->pdf("pdfPsi2s");
		yieldPsi2s = ws->var("yieldPsi2s");
	}
	
	RooRealVar* bExc = ws->var("bExc");
	RooRealVar* bDiss = ws->var("bDiss");
	RooRealVar* nDiss = nullptr;
	if (!exp) nDiss = ws->var("nDiss");
	
	// Define mass frame
	RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
	data->plotOn(mframe, Binning(50));
	fitModel->plotOn(mframe, Name("sum"), LineColor(kRed), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(2), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(1));
	if (mMax > mLimitPsi2s) fitModel->plotOn(mframe,Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kYellow+3), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(kYellow+4), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(kMagenta+2), LineWidth(1));
	
	// Define pt frame
	vector<RooPlot*> ptframes= {pt->frame(Title("Fit of p_{T} (log scale and full p_{T} range)")), pt->frame(Title("Fit of pt (zoom)"))};
	for (int k = 0; k<(int)ptframes.size(); k++) {
		data->plotOn(ptframes[k], Binning(ptBinNumber));
		fitModel->plotOn(ptframes[k], Name("sum"), LineColor(kRed), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(2), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(1));
		if (mMax > mLimitPsi2s) fitModel->plotOn(ptframes[k],Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kYellow+3), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(kYellow+4), LineWidth(1));
		fitModel->plotOn(ptframes[k],Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(kMagenta+2), LineWidth(1));
	}
	
	TCanvas* c1 = new TCanvas("2Dplot","2D fit",1800,1000) ;
	//TCanvas* c1 = new TCanvas("2Dplot","2D fit",800,300) ;
	if (drawPulls) c1->Divide(4,3) ;
	else {c1->Divide(4,2) ; c1->SetCanvasSize(1800, 600);}
	
	int nDof = fitModel->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
	std::cout << std::endl << std::endl << "Number of degrees of freedom = " << nDof << std::endl << std::endl << std::endl;
	// Mass Plot
	c1->cd(2) ;
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	if (logScale) gPad->SetLogy() ;
	double xChi = mMin + (mMax-mMin)*0.7;
	TLatex* txtChi = new TLatex(xChi, 0.5*mframe->GetMaximum(),Form("#chi^{2}/ndf = %.3f", mframe->chiSquare( "sum", "h_data", nDof)));
	mframe->addObject(txtChi) ;
	mframe->chiSquare() ;
	mframe->Draw();
	//std::cout << std::endl << "chi2 = " << chi2 << std::endl;
	
	// pt plot
	c1->cd(3) ;
	gPad->SetLogy();
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	//if (logScale) gPad->SetLogy() ;
	double yMax2 = ptframes[0]->GetMaximum();
	double y1 = 0.75*yMax2, y2 = 0.65*yMax2, y3 = 0.57*yMax2, y4 = 0.5*yMax2;
	if (logScale) {y1 = yMax2/pow(2.,1), y2 = yMax2/pow(2.,2), y3 = yMax2/pow(2.,2.8), y4 = yMax2/pow(2.,4);}
	double ptRangeMax=ptMax;
	ptframes[0]->GetXaxis()->SetRangeUser(0, ptRangeMax);
	ptframes[0]->Draw();
	ptframes[0]->Print();
	
	// Zoom pt
	c1->cd(4);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	if (!exclusiveOnly) ptRangeMax = 3;
	double xText = ptMin+(ptRangeMax-ptMin)*1/2;
	TLatex* txtChi2 = new TLatex(xText, y4,Form("#chi^{2}/ndf = %.3f", ptframes[0]->chiSquare("sum", "h_data", nDof)));
	TLatex* txtExc = new TLatex(xText, y1,Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
	ptframes[1]->addObject(txtExc) ;
	ptframes[1]->addObject(txtChi2) ;
	ptframes[1]->GetXaxis()->SetRangeUser(0, ptRangeMax);
	ptframes[1]->Draw();
	
	/*
	 // Quality plots: (data-fit)/sigma
	 c1->cd(6);
	 gPad->SetLeftMargin(0.15) ;
	 //gPad->SetTopMargin(0.15) ;
	 data->plotOn(mframe, Binning(50));
	 fitModel->plotOn(mframe, Binning(50));
	 RooHist* hpullM = mframe->pullHist();
	 hpullM->GetXaxis()->SetRangeUser(mMin, mMax);
	 hpullM->SetTitle("(data - fit)/#sigma for m distribution");
	 hpullM->Draw("");
	 
	 c1->cd(7);
	 gPad->SetLeftMargin(0.15) ;
	 data->plotOn(ptframes[k], Binning(ptBinNumber));
	 fitModel->plotOn(ptframes[0], Binning(ptBinNumber));
	 RooHist* hpullPt = ptframes[0]->pullHist();
	 //hpullPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
	 hpullPt->GetXaxis()->SetRangeUser(0,ptRangeMax);
	 hpullPt->SetTitle("(data - fit)/#sigma for p_{T} distribution");
	 hpullPt->Draw("");
	 */
	
	
	//Correlation coefficients
	c1->cd(6);
	gStyle->SetTextSize(0.05);
	RooArgList argList = r->floatParsInit();
	bool cor = false;
	int l = 0;
	vector <TText*> txtVec = {};
	for (int i = 0; i<nDof; i++) {
		RooAbsArg* arg1 = argList.at(i);
		for (int j = i+1; j < nDof; j++) {
			RooAbsArg* arg2 = argList.at(j);
			double correl = r->correlation(*arg1, *arg2);
			if (abs(correl) > 0.3) {
				cor=true;
				if (abs(correl) > 0.7) gStyle->SetTextColor(kRed);
				else if (abs(correl) > 0.5) gStyle->SetTextColor(kOrange+7);
				else gStyle->SetTextColor(kBlack);
				//cout << "\n\n\n" << Form("Correlation between %s and %s = %f",  arg1->GetName(), arg2->GetName(), correl) << "\n\n\n" << endl;
				l++;
				int l2 = l;
				if (l>=10) {l2=l-10; c1->cd(7);}
				else c1->cd(6);
				TText* txxx = new TText(0.1, 0.9-l2*0.07, Form("%s and %s = %.2f",  arg1->GetName(), arg2->GetName(), correl) );
				//TText txxx(0.2, 0.3, Form("Correlation between %d and %d", i, j) );
				txtVec.push_back(txxx);
				txxx->Draw();
			}
		}
	}
	if (cor) {
		c1->cd(6);
		TText* txxx = new TText(0.1, 0.9, "Correlation between:" );
		txxx->Draw();}
	gStyle->SetTextSize(0.07);
	gStyle->SetTextColor(kBlack);
	
	// Quality plots: data-fit
	if (drawPulls) {
		c1->cd(10);
		gPad->SetLeftMargin(0.15) ;
		RooHist* hresidM = mframe->residHist();
		hresidM->GetXaxis()->SetRangeUser(mMin, mMax);
		hresidM->SetTitle("Residuals (data - fit) for m distribution");
		hresidM->Draw("");
		
		c1->cd(11);
		gPad->SetLeftMargin(0.15) ;
		RooHist* hresidPt = ptframes[0]->residHist();
		//hresidPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
		hresidPt->GetXaxis()->SetRangeUser(0, 4);
		hresidPt->SetTitle("Residuals (data - fit) for p_{T} distribution");
		hresidPt->Draw("");
	}
	
	// Legend in a subcanvas
	c1->cd(1);
	TLegend* legend = new TLegend(0.1, 0.3, 0.9, 0.9);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	//legend->AddEntry(ptframes[k]->findObject("pdfkCohJpsiToMu"), "kCohJpsiToMu","L");
	legend->AddEntry(ptframes[0]->findObject("pdfJpsiExclusive"), "Exclusive J/Psi","L");
	legend->AddEntry(ptframes[0]->findObject("pdfJpsiDissociative"), "Dissociative J/Psi","L");
	if (!exclusiveOnly && period == "LHC16r") legend->AddEntry(ptframes[0]->findObject("pdfJpsiInclusive"), "Inclusive events","L");
	legend->AddEntry(ptframes[0]->findObject("pdfJpsiGammaPb"), "#gamma + Pb","L");
	if (mMax > mLimitPsi2s) legend->AddEntry(ptframes[0]->findObject("pdfPsi2s"), "Psi(2s)","L");
	legend->AddEntry(ptframes[0]->findObject("pdfTwoGamma"), "#gamma#gamma #rightarrow #mu^{+} #mu^{-}","L");
	legend->AddEntry(ptframes[0]->findObject("pdfBackground"), "extra background","L");
	legend->AddEntry(ptframes[0]->findObject("sum"),"sum","L");
	legend->Draw();
	
	// Number of different contributions in a subcanvas
	c1->cd(5);
	
	// Write number of candidates
	//TLatex* txt1 = new TLatex(3.3,0.9*yMax,Form("Coherent Jpsi : %.1f", yieldCohJpsi.getVal()));
	TLatex* txt2 = new TLatex(0.2,0.9,Form("Exclusive J/#Psi : %.1f #pm %.1f", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
	TLatex* txt3 = new TLatex(0.2,0.8,Form("Dissociative J/#Psi : %.1f #pm %.1f", yieldJpsiDissociative->getVal(), yieldJpsiDissociative->getError()));
	TLatex* txt4 = new TLatex(0.2,0.7,Form("#gamma-Pb J/#Psi : %.1f #pm %.1f", yieldJpsiGammaPb->getVal(), yieldJpsiGammaPb->getError()));
	TLatex* txt4bis = new TLatex(0.2,0.6,Form("Inclusive J/#Psi : %.1f #pm %.1f", yieldJpsiInclusive->getVal(), yieldJpsiInclusive->getError()));
	double yBkg = 0.5;
	if (mMax > mLimitPsi2s) {
		TLatex* txt5 = new TLatex(0.2,yBkg,Form("#Psi(2s) : %.1f #pm %.1f", yieldPsi2s->getVal(), yieldPsi2s->getError()));
		txt5->Draw();
		yBkg = 0.4;
	}
	TLatex* txt6 = new TLatex(0.2,yBkg,Form("#gamma#gamma #rightarrow #mu^{+} #mu^{-} : %.1f #pm %.1f", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
	TLatex* txt9 = new TLatex(0.2,yBkg-0.1,Form("Other background : %.1f #pm %.1f", yieldBkg->getVal(), yieldBkg->getError()));
	txt2->Draw(); txt3->Draw(); txt4->Draw(); txt6->Draw(); txt9->Draw();
	if (!exclusiveOnly && period == "LHC16r") txt4bis->Draw();
	
	// Compute ratios N_diss/N_exc and N_(gamma-Pb)/N_exc
	RooFormulaVar r_diss_exc("r_diss_exc","yieldJpsiDissociative/yieldJpsiExclusive",RooArgSet(*yieldJpsiDissociative, *yieldJpsiExclusive));
	RooFormulaVar r_gammaPb_exc("r_gammaPb_exc","yieldJpsiGammaPb/yieldJpsiExclusive",RooArgSet(*yieldJpsiGammaPb, *yieldJpsiExclusive));
	string percent = "%";
	TLatex* txt7 = new TLatex(0.2,0.3,Form("N_{diss}/N_{exc} = %.2f #pm %.2f %s", r_diss_exc.getVal()*100, r_diss_exc.getPropagatedError(*r)*100, percent.c_str()));
	TLatex* txt8 = new TLatex(0.2,0.2,Form("N_{#gamma-Pb}/N_{exc} = %.2f #pm %.2f %s", r_gammaPb_exc.getVal()*100, r_gammaPb_exc.getPropagatedError(*r)*100, percent.c_str()));
	txt7->Draw(); txt8->Draw();
	ws->import(r_diss_exc);
	
	c1->cd(8);
	if (exclusiveOnly || period == "LHC16s") {
		TLatex* tInc0 = new TLatex(0.1, 0.9, "No inclusive contribution");
		tInc0->Draw();
	}
	else {
		RooRealVar* nInc = ws->var("nInc");
		RooRealVar* pt0 = ws->var("pt0");
		TLatex* tInc0 = new TLatex(0.1, 0.9, "Inclusive: #frac{p_{T}}{(1 + (p_{T}/p_{0})^{2})^{n}} with");
		TLatex* tInc1 = new TLatex(0.1, 0.75, Form("n_{inc} = %.2f #pm %.2f", nInc->getVal(), nInc->getError()));
		TLatex* tInc2 = new TLatex(0.1, 0.65, Form("p_{0} = %.2f #pm %.2f", pt0->getVal(), pt0->getError()));
		tInc0->Draw(); tInc1->Draw(); tInc2->Draw();
	}
	
	
	TLatex* tDiss0 = nullptr;
	TLatex* tDiss1 = new TLatex(0.1, 0.35,Form("b_{diss} = %.2f #pm %.2f", bDiss->getVal(), bDiss->getError()));
	if (exp) {
		tDiss0 = new TLatex(0.1, 0.45, "Dissociative: p_{T}*exp(-b_{diss}*p_{T}^{2}) with");
	}
	else {
		tDiss0 = new TLatex(0.1, 0.45, "Dissociative: p_{T}*(1+p_{T}^{2}*b_{diss}/n_{diss})^{-n_{diss}} with");
		TLatex* tDiss2 = new TLatex(0.1, 0.25,Form("n_{diss} = %.2f #pm %.2f", nDiss->getVal(), nDiss->getError()));
		if (!(exclusiveOnly && period == "LHC16s")) tDiss2->Draw();
	}
	if (!(exclusiveOnly && period == "LHC16s")) {tDiss0->Draw(); tDiss1->Draw();}
	
	
	// save plot
	string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	string suf;
	if (exp) suf = "exp";
	else suf = "power-law";
	if (exclusiveOnly) suf += "-exclusive-only";
	c1->SaveAs(Form("Plots/%s/Fit-2D%s-%.1f-%.1f-%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str()));
	
	/*
	 RooRealVar* pt0 = ws->var("pt0");
	 RooRealVar* nInc = ws->var("nInc");
	 
	 cout << "\n\n\n\npt0 = " << pt0->getVal() << endl;
	 cout << "nInc = " << nInc->getVal() << "\n\n\n" << endl;
	 */
	
	
	
}

void TwoDPlot(string rootfilePath, string rootfilePathMC, vector<string> periods = {"LHC16r", "LHC16s"}) {
	
	Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
	bool useCuts = true, exp = false, exclusiveOnly = false;
	bool logScale = false, drawPulls = false;
	
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.07);
	gStyle->SetTextSize(.07);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
	gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
	
	const int nPeriod = periods.size();
	
	for (int k = 0; k<nPeriod; k++) {
		string period = periods[k];
		if ( ! Initiate(period, mMin, mMax, ptMin, ptMax, useCuts, logScale, drawPulls, exp, exclusiveOnly)) {cout << "Something wrong at initialisation"; return;}
		
		if (exclusiveOnly) {
			if (period == "LHC16r") ptMax = 1.8;
			else ptMax = 0.8;
		}
		else {
			if (period == "LHC16s") ptMax = 3;
		}
		
		// Define cuts
		std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
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
		//return;
		AddModel(wspace, rootfilePath, rootfilePathMC, period, mMin, mMax, ptMin, ptMax, exp, exclusiveOnly);
		wspace->Print();
		MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, drawPulls, logScale, exp, exclusiveOnly);
		
	}
}


