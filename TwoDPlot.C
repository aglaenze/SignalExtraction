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
#include "ExtendedCrystalBall.h"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;


void AddModel(RooWorkspace* ws, std::string rootfilePath, std::string rootfilePathMC, std::string period, TCut mCut, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool exp) {
	// Define 2D model
	// First define fits in mass
	// Second define fits in pt
	// Then make the product to get 2D PDFs
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
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
	
	// Then Psi(2s)
	// That is kIncohPsi2sToMu in MC data
	double alphaL2 = 1.200001, nL2 = 3.017759, alphaR2 = 2.928444, nR2 = 2.256593;
	if (period == "LHC16s") {alphaL2 = 1.192599; nL2 = 3.118819; alphaR2 = 2.927051; nR2 = 2.055075;}
	
	Double_t factorMean = 1.1902; // maybe to adjust
								  //RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07,0,0.1);
								  // scaling factor Psi(2S) / JPsi
	Double_t factorSigma = 1.05; // maybe to adjust
	RooRealVar mean_scaling("mean_scaling", "", factorMean);
	RooRealVar sigma_scaling("sigma_scaling", "", factorSigma);
	RooFormulaVar mean_psi("mean_psi","mean_jpsi*mean_scaling",RooArgSet(mean_jpsi, mean_scaling));
	RooFormulaVar sigma_psi("sigma_psi", "sigma_jpsi*sigma_scaling", RooArgSet(sigma_jpsi, sigma_scaling));
	//RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07, 0, 0.15);
	RooRealVar alpha_psi_L("alpha_psi_L","alpha_psi_L",alphaL2);
	RooRealVar n_psi_L("n_psi_L","n_psi_L",nL2);
	RooRealVar alpha_psi_R("alpha_psi_R","alpha_psi_R",alphaR2);
	RooRealVar n_psi_R("n_psi_R","n_psi_R",nR2);
	//RooCBShape *psi = new RooCBShape("psi","crystal ball PDF",m,mean_psi,sigma_psi,alpha_psi,n_psi);
	ExtendedCrystalBall *psi = new ExtendedCrystalBall("psi","crystal ball PDF", m,
													   mean_psi, sigma_psi, alpha_psi_L,
													   n_psi_L, alpha_psi_R, n_psi_R);
	
	// Finally background
	//Exponential background for mass
	RooRealVar a1("a1","a1",-2,-10,-0.5);
	RooExponential bkg("exp","exp",m,a1);
	
	// Now pt PDFs
	
	// Contribution from exclusive J/Psis
	
	// Using H1 formula (b is free or not)
	double bExcValue = 4;
	if (period == "LHC16r") {
		if (exp) bExcValue = 3.67;
		else bExcValue = 3.92;
	}
	else if (period == "LHC16s") {
		if (exp) bExcValue = 5.10;
		else bExcValue = 5.48;
	}
	//RooRealVar bExc("bExc","bExc", bExcValue, 2, 10);
	RooRealVar bExc("bExc","bExc", bExcValue);
	// H1 formula
	RooGenericPdf *ptPdfExclusive = new RooGenericPdf("jpsiExc","exclusive jPsi PDF","(2*fTrkTrkPt*exp(-bExc*(fTrkTrkPt**2)))",RooArgSet(pt,bExc)) ;
	
	/*
	 // using template pre-defined in sPlot
	 TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	 TH1F* hPtExclusive = (TH1F*)fTemplates->Get("hPtExclusive");
	 RooDataHist* ptHistExclusive = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hPtExclusive);
	 RooHistPdf* ptPdfExclusive = new RooHistPdf("jpsiExc", "ptExclusive", pt, *ptHistExclusive);
	 */
	/*
	 // using pt template from MC data
	 GetPtHistMC(ws, rootfilePathMC, period, "kIncohJpsiToMu");
	 RooAbsPdf* ptPdfExclusive = ws->pdf("ptkIncohJpsiToMu");
	 ptPdfExclusive->SetName("jpsiExc");
	 */
	
	
	// JPsi Dissociative
	/*
	 // using template pre-defined in sPlot
	 TH1F* hPtDissociative = (TH1F*)fTemplates->Get("hPtDissociative");
	 RooDataHist* ptHistDissociative = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hPtDissociative);
	 RooHistPdf* ptDissociative = new RooHistPdf("ptDissociative", "ptDissociative", pt, *ptHistDissociative);
	 */
	RooGenericPdf *ptDissociative = nullptr;
	RooRealVar* bDiss = nullptr;
	RooRealVar* nDiss = nullptr;
	if (exp) {
		// H1 formula (the first one, with the exp)
		//RooRealVar bDiss("bDiss","bDiss", 0.323027, 0, 2);
		bDiss = new RooRealVar("bDiss","bDiss", 0.323027, 0, 2);
		ptDissociative = new RooGenericPdf("jpsiDiss","Dissociative jPsi PDF","(2*fTrkTrkPt*exp(-bDiss*(fTrkTrkPt**2)))", RooArgSet(pt,*bDiss)) ;
	}
	else {
		// H1 formula (the second one, with the power law)
		bDiss = new RooRealVar("bDiss","bDiss", 2, 0, 7);
		nDiss = new RooRealVar("nDiss","nDiss", 4, 1, 10);
		ptDissociative = new RooGenericPdf("jpsiDiss","Dissociative jPsi PDF","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bDiss/nDiss))**(-nDiss))", RooArgSet(pt, *nDiss, *bDiss)) ;
	}
	
	
	// Contribution from gamma Pb events
	// using pt template from MC data
	GetPtHistMC(ws, rootfilePathMC, period, "kCohJpsiToMu");
	RooAbsPdf* ptGammaPb = ws->pdf("ptkCohJpsiToMu");
	ptGammaPb->SetName("ptGammaPb");
	
	// Psi(2s)
	GetPtHistMC(ws, rootfilePathMC, period, "kIncohPsi2sToMu");
	RooAbsPdf* ptPdfPsi2s = ws->pdf("ptkIncohPsi2sToMu");
	ptPdfPsi2s->SetName("ptPdfPsi2s");
	
	// Background
	
	// pt distribution from background is obtained with sPlot
	TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	//TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBackground");
	TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBkgSmooth");
	RooDataHist* ptHistBackground = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), hPtBackground);
	RooHistPdf* ptBackground = new RooHistPdf("ptBackground", "ptBackground", pt, *ptHistBackground);
	
	/*
	 // using sidebands
	 GetSidebandsTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptBackground = ws->pdf("ptSidebands");
	 ptBackground->SetName("ptBackground");
	 */
	
	// Inclusive events
	/*
	 // From data with additional cuts
	 GetInclusiveTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptInclusive = ws->pdf("ptInclusive");
	 ptInclusive->SetName("ptInclusive");
	 */
	
	// With a new formula
	double pt0Val, nIncVal;
	if (period == "LHC16r") {pt0Val = 2.52; nIncVal = 2.71;}
	else {pt0Val = 0.94; nIncVal = 2.71;}
	RooRealVar *pt0 = new RooRealVar("pt0","pt0", pt0Val);
	RooRealVar *nInc = new RooRealVar("nInc","nInc", nIncVal);
	/*
	 nIncVal = 2.22749; pt0Val = 1.44421;
	 RooRealVar *pt0 = new RooRealVar("pt0","pt0", pt0Val, 0.7, 6);
	 RooRealVar *nInc = new RooRealVar("nInc","nInc", nIncVal, 2, 8);
	 */
	RooGenericPdf* ptInclusive = new RooGenericPdf("ptInclusive","Inclusive jPsi PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(pt, *pt0, *nInc)) ;
	
	// Last step:
	// Product fit(m) x fit(pt)
	RooProdPdf* pdfJpsiExclusive = new RooProdPdf("pdfJpsiExclusive","jpsi*ptkIncohJpsiToMu",RooArgList(*jpsi,*ptPdfExclusive));
	RooProdPdf* pdfJpsiDissociative = new RooProdPdf("pdfJpsiDissociative","jpsi*ptDissociative",RooArgList(*jpsi,*ptDissociative));
	RooProdPdf* pdfJpsiGammaPb = new RooProdPdf("pdfJpsiGammaPb","jpsi*ptGammaPb",RooArgList(*jpsi,*ptGammaPb));
	RooProdPdf* pdfPsi2s = new RooProdPdf("pdfPsi2s","psi*ptkIncohPsi2sToMu",RooArgList(*psi,*ptPdfPsi2s));
	RooProdPdf* pdfBackground = new RooProdPdf("pdfBackground","bkg*ptBackground",RooArgList(bkg,*ptBackground));
	RooProdPdf* pdfJpsiInclusive = new RooProdPdf("pdfJpsiInclusive","jpsi*ptInclusive",RooArgList(*jpsi,*ptInclusive));
	
	// All yields
	RooRealVar yieldJpsiExclusive("yieldJpsiExclusive","yieldJpsiExclusive",1000, 1, 3.e3);
	//RooRealVar yieldJpsiExclusive("yieldJpsiExclusive","yieldJpsiExclusive",0);
	RooRealVar yieldJpsiDissociative("yieldJpsiDissociative","yieldJpsiDissociative",100, 1, 3.e3);
	//RooRealVar yieldJpsiDissociative("yieldJpsiDissociative","yieldJpsiDissociative",0);
	RooRealVar yieldJpsiGammaPb("yieldJpsiGammaPb","yieldJpsiGammaPb",10, 1, 5.e2);
	double psi2svalue = 150;
	RooRealVar yieldPsi2s("yieldPsi2s","yieldPsi2s",10, 0, 1.e3);
	if (mMax > mLimitPsi2s) yieldPsi2s.setVal(0);
	RooRealVar yieldBkg("yieldBkg","yieldBkg",500,0.,2.e3);
	//RooRealVar yieldBkg("yieldBkg","yieldBkg",122);
	/*
	RooRealVar yieldJpsiInclusive("yieldJpsiInclusive","yieldJpsiInclusive", 500, 0.,3.e3);
	if (period == "LHC16r") {yieldJpsiInclusive.setRange(30.,3.e3); }
	else if (period == "LHC16s") {yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();}
	 */
	RooRealVar yieldJpsiInclusive("yieldJpsiInclusive","yieldJpsiInclusive",0);
	
	/*
	 // Assemble all components in sets
	 RooArgList* pdfList;
	 RooArgList yieldList;
	 if (mMax > mLimitPsi2s) {
	 pdfList = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfPsi2s, *pdfBackground);
	 yieldList = RooArgList(yieldJpsiExclusive, yieldJpsiDissociative, yieldJpsiGammaPb, yieldPsi2s, yieldBkg);
	 }
	 else {
	 pdfList = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfBackground);
	 yieldList = RooArgList(yieldJpsiExclusive, yieldJpsiDissociative, yieldJpsiGammaPb, yieldBkg);
	 }
	 */
	RooArgList* pdfList = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfPsi2s, *pdfBackground, *pdfJpsiInclusive);
	RooArgList yieldList = RooArgList(yieldJpsiExclusive, yieldJpsiDissociative, yieldJpsiGammaPb, yieldPsi2s, yieldBkg, yieldJpsiInclusive);
	RooArgList* pdfList2 = new RooArgList(*pdfJpsiExclusive, *pdfJpsiDissociative, *pdfJpsiGammaPb, *pdfPsi2s, *pdfBackground, *pdfJpsiInclusive);
	RooArgList yieldList2 = RooArgList(yieldJpsiExclusive, yieldJpsiDissociative, yieldJpsiGammaPb, yieldPsi2s, yieldBkg, yieldJpsiInclusive);
	// Create fit model
	RooAbsPdf* fitModel;
	if (mMax > mLimitPsi2s) fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
	else fitModel = new RooAddPdf("model", "model", *pdfList2, yieldList2, kFALSE);
	
	ws->import(*fitModel);
}

void MakePlots(RooWorkspace* ws, std::string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool drawPulls, bool logScale, bool exp) {
	
	int ptBinNumber = int(10*ptMax);
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooAbsPdf* fitModel = ws->pdf("model");
	
	// Fit data
	RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
	//RooFitResult* r = nullptr;
	
	RooAbsPdf* pdfJpsiExclusive = ws->pdf("pdfJpsiExclusive");
	RooAbsPdf* pdfJpsiDissociative = ws->pdf("pdfJpsiDissociative");
	RooAbsPdf* pdfJpsiGammaPb = ws->pdf("pdfJpsiGammaPb");
	RooAbsPdf* pdfBackground = ws->pdf("pdfBackground");
	RooAbsPdf* pdfJpsiInclusive = ws->pdf("pdfJpsiInclusive");
	
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
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
	fitModel->plotOn(mframe,Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(3), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(1));
	if (mMax > mLimitPsi2s) fitModel->plotOn(mframe,Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(7), LineWidth(1));
	fitModel->plotOn(mframe,Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(12), LineWidth(1));
	
	// Define pt frame
	RooPlot* ptframe = pt->frame(Title("Fit of pt"));
	data->plotOn(ptframe, Binning(ptBinNumber));
	fitModel->plotOn(ptframe, Name("sum"), LineColor(kRed), LineWidth(1));
	fitModel->plotOn(ptframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(2), LineWidth(1));
	fitModel->plotOn(ptframe,Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(3), LineWidth(1));
	fitModel->plotOn(ptframe,Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(1));
	if (mMax > mLimitPsi2s) fitModel->plotOn(ptframe,Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(1));
	fitModel->plotOn(ptframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(7), LineWidth(1));
	fitModel->plotOn(ptframe,Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(12), LineWidth(1));
	
	TCanvas* c1 = new TCanvas("2Dplot","2D fit",1800,1200) ;
	//TCanvas* c1 = new TCanvas("2Dplot","2D fit",800,300) ;
	if (drawPulls) c1->Divide(3,3) ;
	else {c1->Divide(3,2) ; c1->SetCanvasSize(1800, 800);}
	
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
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	if (logScale) gPad->SetLogy() ;
	double yMax2 = ptframe->GetMaximum();
	double y1 = 0.75*yMax2, y2 = 0.65*yMax2, y3 = 0.57*yMax2, y4 = 0.4*yMax2;
	if (logScale) {y1 = yMax2/pow(2.,1), y2 = yMax2/pow(2.,2), y3 = yMax2/pow(2.,2.8), y4 = yMax2/pow(2.,4);}
	double ptRangeMax=ptMax;
	ptRangeMax = 4.;
	double xText = ptMin+(ptRangeMax-ptMin)*1/2;
	TLatex* txtExc = new TLatex(xText, y1,Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
	TLatex* txtDiss = new TLatex(xText, y2,Form("b_{diss} = %.2f #pm %.2f", bDiss->getVal(), bDiss->getError()));
	ptframe->addObject(txtExc) ;
	ptframe->addObject(txtDiss) ;
	if (!exp) {
		TLatex* txtDiss2 = new TLatex(xText, y3,Form("n_{diss} = %.2f #pm %.2f", nDiss->getVal(), nDiss->getError()));
		ptframe->addObject(txtDiss2) ;
	}
	TLatex* txtChi2 = new TLatex(xText, y4,Form("#chi^{2}/ndf = %.3f", ptframe->chiSquare("sum", "h_data", nDof)));
	ptframe->addObject(txtChi2) ;
	ptframe->GetXaxis()->SetRangeUser(0, ptRangeMax);
	ptframe->Draw();
	ptframe->Print();
	
	/*
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
	 //hpullPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
	 hpullPt->GetXaxis()->SetRangeUser(0,ptRangeMax);
	 hpullPt->SetTitle("(data - fit)/#sigma for p_{T} distribution");
	 hpullPt->Draw("");
	 */
	
	//Correlation coefficients
	c1->cd(5);
	gStyle->SetTextSize(0.05);
	TText* txxx = new TText(0.1, 0.9, "Correlation between:" );
	txxx->Draw();
	RooArgList argList = r->floatParsInit();
	int l = 0;
	vector <TText*> txtVec = {};
	for (int i = 0; i<nDof; i++) {
		RooAbsArg* arg1 = argList.at(i);
		for (int j = i+1; j < nDof; j++) {
			RooAbsArg* arg2 = argList.at(j);
			double correl = r->correlation(*arg1, *arg2);
			if (abs(correl) > 0.3) {
				if (abs(correl) > 0.7) gStyle->SetTextColor(kRed);
				else if (abs(correl) > 0.5) gStyle->SetTextColor(kOrange+7);
				else gStyle->SetTextColor(kBlack);
				//cout << "\n\n\n" << Form("Correlation between %s and %s = %f",  arg1->GetName(), arg2->GetName(), correl) << "\n\n\n" << endl;
				l++;
				int l2 = l;
				if (l>=9) {l2=l-9; c1->cd(6);}
				else c1->cd(5);
				TText* txxx = new TText(0.1, 0.9-l2*0.1, Form("%s and %s = %.2f",  arg1->GetName(), arg2->GetName(), correl) );
				//TText txxx(0.2, 0.3, Form("Correlation between %d and %d", i, j) );
				txtVec.push_back(txxx);
				txxx->Draw();
			}
		}
	}
	gStyle->SetTextSize(0.07);
	gStyle->SetTextColor(kBlack);
	
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
	//legend->AddEntry(ptframe->findObject("pdfkCohJpsiToMu"), "kCohJpsiToMu","L");
	legend->AddEntry(ptframe->findObject("pdfJpsiExclusive"), "Exclusive J/Psi","L");
	legend->AddEntry(ptframe->findObject("pdfJpsiDissociative"), "Dissociative J/Psi","L");
	legend->AddEntry(ptframe->findObject("pdfJpsiInclusive"), "Inclusive events","L");
	legend->AddEntry(ptframe->findObject("pdfJpsiGammaPb"), "#gamma + Pb","L");
	if (mMax > mLimitPsi2s) legend->AddEntry(ptframe->findObject("pdfPsi2s"), "Psi(2s)","L");
	legend->AddEntry(ptframe->findObject("pdfBackground"), "#gamma#gamma #rightarrow #mu^{+} #mu^{-} + other background","L");
	legend->AddEntry(ptframe->findObject("sum"),"sum","L");
	legend->Draw();
	
	// Number of different contributions in a subcanvas
	c1->cd(4);
	
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
	TLatex* txt6 = new TLatex(0.2,yBkg,Form("#gamma#gamma #rightarrow #mu^{+} #mu^{-} : %.1f #pm %.1f", yieldBkg->getVal(), yieldBkg->getError()));
	txt2->Draw(); txt3->Draw(); txt4->Draw(); txt6->Draw(); txt4bis->Draw();
	
	// Compute ratios N_diss/N_exc and N_(gamma-Pb)/N_exc
	RooFormulaVar r_diss_exc("r_diss_exc","yieldJpsiDissociative/yieldJpsiExclusive",RooArgSet(*yieldJpsiDissociative, *yieldJpsiExclusive));
	RooFormulaVar r_gammaPb_exc("r_gammaPb_exc","yieldJpsiGammaPb/yieldJpsiExclusive",RooArgSet(*yieldJpsiGammaPb, *yieldJpsiExclusive));
	std::string percent = "%";
	TLatex* txt7 = new TLatex(0.2,0.3,Form("N_{diss}/N_{exc} = %.2f #pm %.2f %s", r_diss_exc.getVal()*100, r_diss_exc.getPropagatedError(*r)*100, percent.c_str()));
	TLatex* txt8 = new TLatex(0.2,0.2,Form("N_{#gamma-Pb}/N_{exc} = %.2f #pm %.2f %s", r_gammaPb_exc.getVal()*100, r_gammaPb_exc.getPropagatedError(*r)*100, percent.c_str()));
	txt7->Draw(); txt8->Draw();
	
	// save plot
	std::string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	c1->SaveAs(Form("Plots/Fit-2D-%s%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax));
	
	/*
	 RooRealVar* pt0 = ws->var("pt0");
	 RooRealVar* nInc = ws->var("nInc");
	 
	 cout << "\n\n\n\npt0 = " << pt0->getVal() << endl;
	 cout << "nInc = " << nInc->getVal() << "\n\n\n" << endl;
	 */
	
	
	
}

void TwoDPlot(std::string rootfilePath, std::string rootfilePathMC, std::vector<std::string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2., Double_t mMax = 4.2, Double_t ptMin = 0., Double_t ptMax = 4, bool useCuts = true, bool logScale = false, bool exp = false, bool drawPulls = false) {
	
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
	
	const int nPeriod = periods.size();
	
	for (int k = 0; k<nPeriod; k++) {
		std::string period = periods[k];
		// Define cuts
		std::list<TCut> mCutList = DefineCuts(period);
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
		AddModel(wspace, rootfilePath, rootfilePathMC, period, mCut, mMin, mMax, ptMin, ptMax, exp);
		wspace->Print();
		MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, drawPulls, logScale, exp);
		
	}
}


