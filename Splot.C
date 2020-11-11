#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include "Include/FitUtils.C"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;


void SubtractBkg(TH1* hPtMc , TH1* hBackground, TH1* hSub) {
	
	int nBins = hBackground->GetNbinsX();

	double scale = (double)hBackground->GetMaximum()/hPtMc->GetMaximum();
	hPtMc->Scale(scale);
	
	// Create the subtraction histo between packground data and MC (two gamma)
	//TH1F* hSub = new TH1F("histSub", "hist sub", nBins, xMin, xMax);
	for (int l = 0; l<nBins; l++) {
		int binBkg = hBackground->GetBinContent(l);
		int binMc = hPtMc->GetBinContent(l);
		if (binBkg > binMc) hSub->SetBinContent(l, binBkg - binMc);
		else hSub->SetBinContent(l, 0);
	}
}

void SmoothHistLog(TH1* hPtBackground, TH1* hPtBkgSmooth, int avgSize = 3) {
	double negValue = -10;	// in case the value is negative, pb log(neg) so set log to this negvalue
	int ptBinNumber = hPtBkgSmooth->GetNbinsX();
	for (int k = 0; k<ptBinNumber; k++) {
		int binNum = k+1;
		double smoothVal = 0;
		int nSmoothBins = 1;
		if (hPtBackground->GetBinContent(binNum) > 0) {smoothVal = TMath::Log(hPtBackground->GetBinContent(binNum));}
		else smoothVal = -10;	// by default say there's 10^-2 background events (originally) in this bin
		double size = binNum*0.8;	// number of bins to smooth
		if (size > avgSize) size = avgSize;
		for (int j = 1; j < (int)size; j++ ) {
			if (binNum-j > 0) {
				nSmoothBins++;
				double neighbourVal = hPtBackground->GetBinContent(binNum-j);
				if (neighbourVal>0) {smoothVal+= TMath::Log(neighbourVal);}
				else smoothVal+= negValue;
			}
			if (binNum+j < ptBinNumber) {
				nSmoothBins++;
				double neighbourVal = hPtBackground->GetBinContent(binNum+j);
				if (neighbourVal>0) {smoothVal+= TMath::Log(hPtBackground->GetBinContent(binNum+j));}
				else smoothVal+= negValue;
			}
		}
		smoothVal /= (double)nSmoothBins;
		smoothVal = TMath::Exp(smoothVal);
		Double_t x = hPtBackground->GetXaxis()->GetBinCenter(binNum);
		hPtBkgSmooth->Fill(x, smoothVal);
		//std::cout << x << " " << smoothVal << std::endl;
	}
	hPtBkgSmooth->SetOption("hist");
}

void SmoothHist(TH1* hPtBackground, TH1* hPtBkgSmooth, int avgSize = 3) {
	int ptBinNumber = hPtBkgSmooth->GetNbinsX();
	for (int k = 0; k<ptBinNumber; k++) {
		int binNum = k+1;
		double smoothVal = hPtBackground->GetBinContent(binNum);
		int nSmoothBins = 1;
		double size = binNum*0.8;	// number of bins to smooth
		if (size > avgSize) size = avgSize;
		for (int j = 1; j < (int)size; j++ ) {
			if (binNum-j > 0) {
				nSmoothBins++;
				smoothVal+= hPtBackground->GetBinContent(binNum-j);
			}
			if (binNum+j < ptBinNumber) {
				nSmoothBins++;
				smoothVal+= hPtBackground->GetBinContent(binNum+j);
			}
		}
		smoothVal /= (double)nSmoothBins;
		Double_t x = hPtBackground->GetXaxis()->GetBinCenter(binNum);
		hPtBkgSmooth->Fill(x, smoothVal);
		//std::cout << x << " " << smoothVal << std::endl;
	}
	hPtBkgSmooth->SetOption("hist");
}

void AddModel(RooWorkspace* ws, string period, Double_t mMax) {
	// Define model
	
	//RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)",2.5,3.5);
	RooRealVar m = *ws->var("fTrkTrkM");
	
	LoadMassFitFunctions(ws, period);
	// Now collect m pdf
	RooAbsPdf* jpsi = ws->pdf("jpsi");
	RooAbsPdf* psi2s = ws->pdf("psi2s");
	RooAbsPdf* bkg = ws->pdf("bkg");
	
	//RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
	RooRealVar fsigJpsi("fsigJpsi","signalJPsi",1000,0.,1.e4);
	RooRealVar fsigPsi2s("fsigPsi2s","signalPsi2s",100,0.,1.e3);
	RooRealVar fbkg("fbkg","fbkg",500,0.,1.e7);
	
	RooAbsPdf* model;
	if (mMax > mLimitPsi2s) model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *psi2s, *bkg), RooArgList(fsigJpsi, fsigPsi2s, fbkg), kFALSE);
	else model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *bkg), RooArgList(fsigJpsi, fbkg), kFALSE);
	
	ws->import(*model);
}

void AddPtJpsiModel(RooWorkspace* ws, string rootfilePathMC, string period, bool exp, bool exclusiveOnly) {
	// Define fitting model for pt for J/Psi signal only
	// Dissociative contribution and exclusive contribution
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	LoadJpsiPtFitFunctions(ws, rootfilePathMC, period, exp, exclusiveOnly);
	
	// Then pt pdf
	RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
	RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
	RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
	RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
	
	// Load yields
	LoadYields(ws, period, exclusiveOnly);
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
	RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
	
	RooAbsPdf* model = new RooAddPdf("ptfit", "ptfit", RooArgList(*ptJpsiExclusive, *ptJpsiDissociative, *ptJpsiGammaPb, *ptJpsiInclusive), RooArgList(*yieldJpsiExclusive, *yieldJpsiDissociative, *yieldJpsiGammaPb , *yieldJpsiInclusive), kFALSE);
	
	ws->import(*model);
}

void DoSPlot(RooWorkspace* ws, Double_t mMax) {
	
	RooAbsPdf* model = ws->pdf("mfit");
	RooRealVar* jPsiYield = ws->var("fsigJpsi");
	RooRealVar* psi2sYield = ws->var("fsigPsi2s");
	RooRealVar* bkgYield = ws->var("fbkg");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	model->fitTo(*data, Extended(), Minos(true), Strategy(2));
	//The sPlot technique requires that we fix the parameters
	// of the model that are not yields after doing the fit
	// J/Psi
	RooRealVar* mean_jpsi = ws->var("mean_jpsi");
	RooRealVar* sigma_jpsi = ws->var("sigma_jpsi");
	//RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
	//RooRealVar* n_jpsi = ws->var("n_jpsi");
	RooRealVar* alpha_jpsi_L = ws->var("alpha_jpsi_L");
	RooRealVar* n_jpsi_L = ws->var("n_jpsi_L");
	RooRealVar* alpha_jpsi_R = ws->var("alpha_jpsi_R");
	RooRealVar* n_jpsi_R = ws->var("n_jpsi_R");
	
	// Psi(2s)
	RooRealVar *mean_psi2s, *sigma_psi2s, *alpha_psi2s_L, *n_psi2s_L, *alpha_psi2s_R, *n_psi2s_R;
	if (mMax > mLimitPsi2s) {
		mean_psi2s = ws->var("mean_psi2s");
		sigma_psi2s = ws->var("sigma_psi2s");
		alpha_psi2s_L = ws->var("alpha_psi2s_L");
		n_psi2s_L = ws->var("n_psi2s_L");
		alpha_psi2s_R = ws->var("alpha_psi2s_R");
		n_psi2s_R = ws->var("n_psi2s_R");
		alpha_psi2s_L->setConstant();
		n_psi2s_L->setConstant();
		alpha_psi2s_R->setConstant();
		n_psi2s_R->setConstant();
	}
	
	// Background
	RooRealVar* a1 = ws->var("a1");
	
	mean_jpsi->setConstant();
	sigma_jpsi->setConstant();
	alpha_jpsi_L->setConstant();
	n_jpsi_L->setConstant();
	alpha_jpsi_R->setConstant();
	n_jpsi_R->setConstant();
	
	a1->setConstant();
	
	RooMsgService::instance().setSilentMode(true);
	
	//Now we use the SPlot class to add SWeight to our data set
	// based on our model and our yield variables
	
	RooStats::SPlot * sData;
	if (mMax > mLimitPsi2s) sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *psi2sYield, *bkgYield) );
	else sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *bkgYield) );
	
	//Check Sweight properties
	std::cout << "Check SWeights: " << std::endl;
	
	std::cout << std::endl << "Yield of JPsi is "
	<< jPsiYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fsigJpsi") << std::endl;
	
	
	std::cout << std::endl << "Yield of bkg is "
	<< bkgYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fbkg") << std::endl;
	
	if (mMax > mLimitPsi2s) {
		std::cout << std::endl << "Yield of Psi(2s) is "
		<< psi2sYield->getVal() << ". From sWeights it is "
		<< sData->GetYieldFromSWeight("fsigPsi2s") << std::endl;
		for (Int_t i=0; i < 10; i ++ )
		{
			std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
			<< " psi(2s) Yield "<< sData->GetSWeight(i, "fsigPsi2s")
			<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
			<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
			<< std::endl;
		}
	}
	else {
		for (Int_t i=0; i < 10; i ++ )
		{
			std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
			<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
			<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
			<< std::endl;
		}
	}
	
	//import the new data set with Sweight
	
	std::cout << "import new dataset with sWeight" << std::endl;
	ws->import(*data, Rename("dataWithSWeights"));
	
	std::cout << "Sigma = " << sigma_jpsi->getVal() << "\n\n\n\n\n" << std::endl;
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exp, bool exclusiveOnly, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax){
	
	
	int ptBinNumber = int(10*ptMax);
	
	//make some plots
	TCanvas* cv = new TCanvas("splot","splot", 800, 800) ;
	cv->Divide(3,3);
	
	//get what we need of the workspace
	RooAbsPdf* model = ws->pdf("mfit");
	RooAbsPdf* jPsiModel = ws->pdf("jpsi");
	RooAbsPdf* bkgModel = ws->pdf("bkg");
	
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	
	RooDataSet* sdata = (RooDataSet*) ws->data("dataWithSWeights");
	
	RooRealVar* jPsiYield = ws->var("fsigJpsi");
	RooRealVar* psi2sYield = ws->var("fsigPsi2s");
	RooRealVar* bkgYield = ws->var("fbkg");
	
	// Draw mass fit
	cv->cd(1);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frame = m->frame();
	sdata->plotOn(frame, Binning(50));
	model->plotOn(frame);
	model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kBlue));
	model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen));
	frame->SetTitle("Fit to model to discriminating variable");
	double yMax0 = frame->GetMaximum();
	double xPos = mMin+(mMax-mMin)*0.1;
	//if (mMax < mLimitPsi2s) xPos = mMin+(mMax-mMin)/3;
	TText* txt0 = new TText(xPos,0.7*yMax0,Form("%.1f J/Psi", jPsiYield->getVal()));
	TText* txt1 = new TText(xPos,0.55*yMax0,Form("%.1f bkg", bkgYield->getVal()));
	frame->addObject(txt0) ;
	frame->addObject(txt1) ;
	if (mMax > mLimitPsi2s) {
		RooAbsPdf* psi2sModel = ws->pdf("psi2s");
		model->plotOn(frame, Components(*psi2sModel), LineStyle(kDashed), LineColor(kOrange));
		TText* txtPsi2s = new TText((mMin+mMax+0.1)/2,0.45*yMax0,Form("%.1f Psi(2s)", psi2sYield->getVal()));
		frame->addObject(txtPsi2s) ;
	}
	frame->Draw();
	
	// Draw quality plots
	cv->cd(4);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	sdata->plotOn(frame, Binning(ptBinNumber));
	model->plotOn(frame);
	RooHist* hpull = frame->pullHist();
	hpull->GetXaxis()->SetRangeUser(mMin, mMax);
	hpull->SetTitle("(data - fit)/#sigma of the mass fit");
	hpull->Draw("");
	cv->cd(5);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooHist* hresid = frame->residHist();
	hresid->GetXaxis()->SetRangeUser(mMin, mMax);
	hresid->SetTitle("Residuals (data - fit) of the mass fit");
	hresid->Draw("");
	
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	// create weighted data set for JPsi
	RooDataSet * dataw_jpsi = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fsigJpsi_sw") ;
	//RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),"fTrkTrkM < 3.2 && fTrkTrkM > 2.8","fsigJpsi_sw") ;
	// create weighted data set for Background
	RooDataSet * dataw_bkg = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fbkg_sw") ;
	
	
	// Draw J/Psi pt with weights
	cv->cd(2);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frameJpsi = pt->frame() ;
	int binningNum = 100;
	if (exclusiveOnly) {binningNum = 40;}
	dataw_jpsi->plotOn(frameJpsi, Name("jPsiData"), DataError(RooAbsData::SumW2), Binning(binningNum) ) ;
	frameJpsi->SetTitle("p_{T} distribution for J/#Psi with weights");
	
	// Fit pt distrib of J/Psi with exclusive / dissociative components
	RooAbsPdf* ptmodel = ws->pdf("ptfit");
	RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
	RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
	RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
	RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
	ptmodel->fitTo(*dataw_jpsi, Extended(), Minos(true), Strategy(1), SumW2Error(true));
	ptmodel->plotOn(frameJpsi, Name("ptfit"));
	ptmodel->plotOn(frameJpsi, Components(*ptJpsiExclusive), LineStyle(kDashed), LineColor(kRed));
	ptmodel->plotOn(frameJpsi, Components(*ptJpsiDissociative), LineStyle(kDashed), LineColor(kGreen));
	ptmodel->plotOn(frameJpsi, Components(*ptJpsiGammaPb), LineStyle(kDashed), LineColor(kBlue));
	ptmodel->plotOn(frameJpsi, Components(*ptJpsiInclusive), LineStyle(kDashed), LineColor(kGray));
	
	// Write chi2
	double xPosJpsi = ptMax/3;
	int nDof = ptmodel->getParameters(dataw_jpsi)->selectByAttrib("Constant",kFALSE)->getSize();
	double yMax = frameJpsi->GetMaximum();
	TLatex* txtChi = new TLatex(xPosJpsi/2, 0.93*yMax,Form("#chi^{2}/ndf = %.3f", frameJpsi->chiSquare( "ptfit", "jPsiData", nDof)));
	frameJpsi->addObject(txtChi);
	cout << endl << endl << endl << endl;
	//frameJpsi->Print();
	
	// Write parameters on plot
	RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
	RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
	RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
	RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
	TLatex* txt2 = new TLatex(xPosJpsi,0.85*yMax,Form("%.1f #pm %.1f exclusive J/Psi", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
	frameJpsi->addObject(txt2) ;
	TLatex* txt = new TLatex(xPosJpsi,0.77*yMax,Form("%.1f #pm %.1f dissociative J/Psi", yieldJpsiDissociative->getVal(), yieldJpsiDissociative->getError()));
	frameJpsi->addObject(txt) ;
	TLatex* txt4 = new TLatex(xPosJpsi,0.69*yMax,Form("%.1f #pm %.1f #gamma-Pb J/Psi", yieldJpsiGammaPb->getVal(), yieldJpsiGammaPb->getError()));
	frameJpsi->addObject(txt4) ;
	TLatex* txt3 = new TLatex(xPosJpsi,0.61*yMax,Form("%.1f #pm %.1f inclusive J/Psi", yieldJpsiInclusive->getVal(), yieldJpsiInclusive->getError()));
	if (!exclusiveOnly && period == "LHC16r") frameJpsi->addObject(txt3) ;
	
	
	// Write fit function parameters on the plot
	double xPosJpsi2 = xPosJpsi*1.2;
	RooRealVar* bExc = ws->var("bExc");
	TLatex* txtExc = new TLatex(xPosJpsi2,0.45*yMax,Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
	frameJpsi->addObject(txtExc) ;
	RooRealVar* bDiss = ws->var("bDiss");
	
	TLatex* txtDiss = new TLatex(xPosJpsi2,0.35*yMax,Form("b_{diss} = %.2f #pm %.2f", bDiss->getVal(), bDiss->getError()));
	frameJpsi->addObject(txtDiss) ;
	
	if (!exp) {
		RooRealVar* nDiss = ws->var("nDiss");
		TLatex* txtDiss2 = new TLatex(xPosJpsi2,0.25*yMax,Form("n_{diss} = %.2f #pm %.2f", nDiss->getVal(), nDiss->getError()));
		frameJpsi->addObject(txtDiss2) ;
	}
	RooRealVar* pt0 = ws->var("pt0");
	RooRealVar* nInc = ws->var("nInc");
	TLatex* txtInc1 = new TLatex(xPosJpsi2,0.15*yMax,Form("p_{0} = %.2f #pm %.2f", pt0->getVal(), pt0->getError()));
	TLatex* txtInc2 = new TLatex(xPosJpsi2,0.05*yMax,Form("n_{inc} = %.2f #pm %.2f", nInc->getVal(), nInc->getError()));
	if (!exclusiveOnly && period == "LHC16r") {frameJpsi->addObject(txtInc1) ; frameJpsi->addObject(txtInc2) ;}
	
	
	frameJpsi->Draw() ;
	
	// Draw Psi(2s) pt with sweights
	if (mMax > mLimitPsi2s) {
		cv->cd(6);
		gPad->SetLeftMargin(0.15) ;
		gPad->SetBottomMargin(0.15) ;
		// create weighted data set for Psi(2s)
		RooDataSet * dataw_psi2s = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fsigPsi2s_sw") ;
		RooPlot* framePsi2s = pt->frame() ;
		dataw_psi2s->plotOn(framePsi2s, DataError(RooAbsData::SumW2), Binning(50) ) ;
		framePsi2s->SetTitle("p_{T} distribution for #Psi(2s) with weights");
		framePsi2s->Draw() ;
	}
	
	// Draw background pt with sweights
	cv->cd(3);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frameBkg = pt->frame() ;
	frameBkg->SetTitle("p_{T} distribution for background with weights");
	frameBkg->GetXaxis()->SetRangeUser(0,3);
	
	TLegend* legend = new TLegend(0.4, 0.6, 0.9, 0.9);
	legend->SetTextSize(0.05);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	
	dataw_bkg->plotOn(frameBkg,DataError(RooAbsData::SumW2), Binning(50)) ;
	
	// Draw fit of pt
	// Add MC templates for (gammma gamma to mu mu) for comparison
	TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
	GetPtHistMC(ws, rootfilePathMc, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
	GetV0Template(ws, rootfilePath, period);
	RooAbsPdf* pdfV0 = ws->pdf("ptV0C7");
	RooAbsPdf* pdfMC = ws->pdf("ptkTwoGammaToMuLow");
	pdfMC->SetName("ptMC");
	RooRealVar* yieldMC = new RooRealVar("yieldMC","yieldMC", 10, 0, 1000);
	RooRealVar* yieldV0 = new RooRealVar("yieldV0","yieldV0", 40, 0, 1000);
	
	RooArgList* pdfList = new RooArgList(*pdfMC, *pdfV0);
	RooArgList yieldList = RooArgList(*yieldMC, *yieldV0);
	
	//RooArgList* pdfList = new RooArgList(*pdfMC);
	//RooArgList yieldList = RooArgList(*yieldMC);
	
	RooAbsPdf* fitModel = new RooAddPdf("modelBkg", "modelBkg", *pdfList, yieldList, kFALSE);
	
	frameBkg->Draw() ;
	frameBkg->Print("v");
	
	legend->Draw("same");
	
	
	// Plot isolation for QCD component.
	// Eg. plot all events weighted by the sWeight for the QCD component.
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	
	
	// Save plots
	string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	string suffix = "";
	if (exp) suffix = "exp";
	else suffix = "power-law";
	if (exclusiveOnly) suffix += "-exclusive-only";
	cv->SaveAs(Form("Plots/%s/Splot%s-%.1f-%.1f-%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suffix.c_str()));
	
	/*
	 // Look at V0CCounts
	 cv->cd(2);
	 gPad->SetLeftMargin(0.15) ;
	 RooRealVar* v0Ccounts = ws->var("fV0CCounts");
	 RooPlot* frame6 = v0Ccounts->frame() ;
	 //dataw_jpsi->plotOn(frame6, DataError(RooAbsData::SumW2) ) ;
	 data->plotOn(frame6, DataError(RooAbsData::SumW2) ) ;
	 
	 frame6->SetTitle("V0C Counts");
	 frame6->Draw() ;
	 
	 // Look at 2D hist of V0CCounts et pt
	 TCanvas* c4 = new TCanvas("ptV0", "ptV0", 600, 600);
	 TH1* hh_data = dataw_jpsi->createHistogram("fTrkTrkPt,fV0CCounts",40,10) ;
	 gPad->SetLeftMargin(0.15) ; gPad->SetRightMargin(0.15) ; hh_data->Draw("colz") ;
	 
	 c4->SaveAs(Form("Plots/Pt-V0-2D-%s.pdf", period.c_str()));
	 */
	
	// And save pt shapes
	string suffix2 = "";
	if (exclusiveOnly) suffix2 = "-exclusive-only";
	TFile* fNewTemplates = new TFile(Form("%s/sPlotTemplates-%s%s.root", rootfilePath.c_str(), period.c_str(), suffix2.c_str()),"RECREATE");
	TH1* hPtDissociative = ptJpsiDissociative->createHistogram("hPtDissociative",*pt, IntrinsicBinning(false), Binning(ptBinNumber, 0, ptMax)) ;
	hPtDissociative->Write();
	
	TH1* hPtExclusive = ptJpsiExclusive->createHistogram("hPtExclusive",*pt, IntrinsicBinning(false), Binning(ptBinNumber, 0, ptMax)) ;
	TH1* hPtBackground = dataw_bkg->createHistogram("hPtBackground",*pt, IntrinsicBinning(false), Binning(ptBinNumber, 0, ptMax)) ;
	
	TH1* hPtMc = pdfMC->createHistogram("hPtMc",*pt, IntrinsicBinning(false), Binning(ptBinNumber, 0, ptMax)) ;
	
	// Create template of extra background = not two gamma to mu
	TH1F* hSubtraction = new TH1F("hSubtraction", "Bkg - MC", ptBinNumber, 0, ptMax);
	SubtractBkg(hPtMc, hPtBackground, hSubtraction);
	TH1* hPtExtraBackground = pdfV0->createHistogram("hPtExtraBackground",*pt, IntrinsicBinning(false), Binning(ptBinNumber, 0, ptMax)) ;
	TH1F* hSubSmooth = new TH1F("hSubSmooth", "Smoothed Bkg - MC", ptBinNumber, 0, ptMax);
	SmoothHist(hSubtraction, hSubSmooth, 4);
	
	// Now write all these histograms
	fNewTemplates->cd();
	hPtBackground->SetOption("hist");
	hPtExtraBackground->SetOption("hist");
	hPtExclusive->Write();
	hPtBackground->Write();
	hPtExtraBackground->Write();
	hSubtraction->Write();
	hSubSmooth->Write();
	hPtMc->SetOption("hist");
	hPtMc->Write();
	
	// Create a smooth histogram (the size of the smooting depends on pt, should be logarithmic so that the smoothing is a lot more agressive for high pt)
	TH1* hPtBkgSmooth = new TH1F("hPtBkgSmooth", "smoothed background pt template", ptBinNumber, 0, ptMax);
	SmoothHistLog(hPtBackground, hPtBkgSmooth, 7);
	
	fNewTemplates->cd();
	hPtBkgSmooth->Write();
	fNewTemplates->Close();
}


void DrawBkgPt(string rootfilePath, string period, bool exclusiveOnly) {
	
	gStyle->SetTextSize(.04);
	
	string suffix = "";
	if (exclusiveOnly) suffix = "-exclusive-only";
	TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s%s.root", rootfilePath.c_str(), period.c_str(), suffix.c_str()),"READ");
	
	TH1* hPtBackground = (TH1*)fTemplates->Get("hPtBackground__fTrkTrkPt");
	hPtBackground->SetTitle("Background p_{T} distribution with sPlot");
	TH1* hSubtraction = (TH1*)fTemplates->Get("hSubtraction");
	TH1* hSubSmooth = (TH1*)fTemplates->Get("hSubSmooth");
	TH1* hPtMc = (TH1*)fTemplates->Get("hPtMc__fTrkTrkPt");
	TCanvas* cv3 = new TCanvas();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	hPtBackground->SetLineColor(kBlack);
	hPtBackground->SetLineWidth(3);
	hPtBackground->Draw("hist");
	hPtMc->SetLineColor(kRed);
	hPtMc->SetLineWidth(2);
	hSubtraction->SetLineColor(kRed+2);
	hSubtraction->SetLineWidth(2);
	hSubSmooth->SetLineColor(kBlue);
	hSubSmooth->SetLineWidth(3);
	hPtMc->Draw("hist same");
	hSubtraction->Draw("hist same");
	hSubSmooth->Draw("hist same");
	/*
	 hPtBkgSmooth->SetLineColor(kRed);
	 hPtBkgSmooth->Draw("hist same");
	 */
	
	// Write number of candidates
	Double_t yMax = hPtBackground->GetMaximum();
	TLatex* t1 = new TLatex(0.5, 0.8*yMax, Form("Number of bkg = %.1f", hPtBackground->Integral()));
	TLatex* t2 = new TLatex(0.5, 0.7*yMax, Form("Number of #gamma#gamma = %.1f", hPtMc->Integral()));
	t1->Draw(); t2->Draw();
	
	TLegend* lgd = new TLegend(0.5, 0.5, 0.9, 0.9);
	lgd->AddEntry(hPtBackground, "Data with sPlot", "l");
	lgd->AddEntry(hPtMc, "MC #gamma#gamma", "l");
	lgd->AddEntry(hSubtraction, "Subtraction: Data - MC", "l");
	lgd->AddEntry(hSubSmooth, "Smoothed subtraction", "l");
	lgd->Draw();
	cv3->SaveAs(Form("Plots/%s/pt-splotsmooth%s.pdf", period.c_str(), suffix.c_str()));
}

//____________________________________
void DrawWeights(RooWorkspace *ws){
	
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	
	RooDataSet* sData = (RooDataSet*) ws->data("dataWithSWeights");
	
	//std::cout << std::endl << std::endl << std::endl << std::endl;
	const int nEntries = (int) sData->sumEntries();
	Double_t massForWeights[nEntries], ptForWeights[nEntries], jPsiWeight[nEntries], bkgWeight[nEntries];
	
	for (Int_t i=0; i < (Int_t)nEntries; i ++ ) {
		massForWeights[i] = sData->get(i)->getRealValue("fTrkTrkM");
		ptForWeights[i] = sData->get(i)->getRealValue("fTrkTrkPt");
		
		jPsiWeight[i] = sData->get(i)->getRealValue("fsigJpsi_sw");
		bkgWeight[i] = sData->get(i)->getRealValue("fbkg_sw");
		
	}
	//std::cout << std::endl << std::endl << std::endl << std::endl;
	
	// TGraphs of sWeights = f(m) or sWeights = f(pt)
	TGraph* gWeightsMJpsi = new TGraph(nEntries,massForWeights, jPsiWeight);
	TGraph* gWeightsPtJpsi = new TGraph(nEntries,ptForWeights, jPsiWeight);
	TGraph* gWeightsMBkg = new TGraph(nEntries,massForWeights, bkgWeight);
	TGraph* gWeightsPtBkg = new TGraph(nEntries,ptForWeights, bkgWeight);
	
	gWeightsMJpsi->SetTitle("JPsi sWeights = f(m)");
	gWeightsPtJpsi->SetTitle("JPsi sWeights = f(pt)");
	gWeightsMBkg->SetTitle("Background sWeights = f(m)");
	gWeightsPtBkg->SetTitle("Background sWeights = f(pt)");
	
	gWeightsMJpsi->GetXaxis()->SetTitle("m");
	gWeightsPtJpsi->GetXaxis()->SetTitle("pt");
	gWeightsMBkg->GetXaxis()->SetTitle("m");
	gWeightsPtBkg->GetXaxis()->SetTitle("pt");
	
	gWeightsMJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
	gWeightsPtJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
	gWeightsMBkg->GetYaxis()->SetTitle("Bkg sWeights");
	gWeightsPtBkg->GetYaxis()->SetTitle("Bkg sWeights");
	
	/*
	 TCanvas* cv = new TCanvas("sweights","sweights", 800, 1200) ;
	 cv->Divide(2,3);
	 
	 cv->cd(1);
	 gPad->SetLeftMargin(0.2);
	 gWeightsMJpsi->Draw("A*");
	 
	 cv->cd(2);
	 gPad->SetLeftMargin(0.2);
	 gWeightsPtJpsi->Draw("A*");
	 
	 cv->cd(3);
	 gPad->SetLeftMargin(0.2);
	 gWeightsMBkg->Draw("A*");
	 
	 cv->cd(4);
	 gPad->SetLeftMargin(0.2);
	 gWeightsPtBkg->Draw("A*");
	 
	 // 2d histograms
	 TH2F* hWeightsBkg2d = new TH2F("hWeightsBkg2d", "Background sWeights)", 50, mMin, mMax, 100, ptMin, 4.5);
	 TH2F* hWeightsJpsi2d = new TH2F("hWeightsJpsi2d", "J/Psi sWeights", 50, mMin, mMax, 100, ptMin, 4.5);
	 for (int i = 0; i<nEntries; i++) {
	 hWeightsBkg2d->Fill(massForWeights[i], ptForWeights[i], bkgWeight[i]);
	 hWeightsJpsi2d->Fill(massForWeights[i], ptForWeights[i], jPsiWeight[i]);
	 }
	 cv->cd(5);
	 gPad->SetLeftMargin(0.2);
	 hWeightsBkg2d->GetXaxis()->SetTitle("Mass");
	 hWeightsBkg2d->GetYaxis()->SetTitle("Pt");
	 //hWeightsBkg2d->Draw("CONTZ");
	 hWeightsBkg2d->Draw("COLZ");
	 
	 cv->cd(6);
	 gPad->SetLeftMargin(0.2);
	 hWeightsJpsi2d->GetXaxis()->SetTitle("Mass");
	 hWeightsJpsi2d->GetYaxis()->SetTitle("Pt");
	 hWeightsJpsi2d->Draw("COLZ");
	 
	 // Save plots
	 string cutType = "";
	 if (!useCuts) cutType = "-nocuts";
	 cv->SaveAs(Form("Plots/Splot-weights-%s%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax));
	 */
}

void Splot(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}) {
	
	Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
	bool useCuts = true, exp = false, exclusiveOnly = false;
	bool logScale = false, drawPulls = false;

	gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.05);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
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
		
		std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
		// Define cuts
		TCut mCut = "";
		
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		if (!useCuts) mCut = "";
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()), "READ");
		
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax);
		AddModel(wspace, period, mMax);
		AddPtJpsiModel(wspace, rootfilePathMc, period, exp, exclusiveOnly);
		wspace->Print();
		
		DoSPlot(wspace, mMax);
		MakePlots(wspace, rootfilePath, rootfilePathMc, period, exp, exclusiveOnly, useCuts, mMin, mMax, ptMin, ptMax);
		//DrawWeights(wspace);
		
		//cleanup
		delete wspace;
		DrawBkgPt(rootfilePath, period, exclusiveOnly);
	}
	
	
}


