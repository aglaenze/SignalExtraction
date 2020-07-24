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
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include "_FitUtils.C"
#include "ExtendedCrystalBall.h"

using namespace RooFit;


void AddModel(RooWorkspace* ws, std::string period) {
	// Define model
	
	//RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)",2.5,3.5);
	RooRealVar m = *ws->var("fTrkTrkM");
	
	//Crystal ball for the J/psi
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3.1,2.9,3.2);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.07,0,0.1);
	/*
	 // case of a simple CB
	 RooRealVar alpha_jpsi("alpha_jpsi","alpha_jpsi",1);
	 RooRealVar n_jpsi("n_jpsi","n_jpsi",5);
	 RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
	 */
	// 2-sided CB
	double alphaL = 0.961839, nL = 7.521515, alphaR = 2.641260, nR = 3.325886;
	if (period == "LHC16s") {alphaL = 0.993482; nL = 6.845735; alphaR = 2.669157; nR = 3.078395;}
	RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L", alphaL);
	RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L", nL);
	RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R", alphaR);
	RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R", nR);
	ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", m,
														mean_jpsi, sigma_jpsi, alpha_jpsi_L,
														n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
	
	// Then Psi(2s)
	double alphaL2 = 1.200001, nL2 = 3.017759, alphaR2 = 2.928444, nR2 = 2.256593;
	if (period == "LHC16s") {alphaL2 = 1.192599; nL2 = 3.118819; alphaR2 = 2.927051; nR2 = 2.055075;}
	
	Double_t factorMean = 1.1902; // maybe to adjust
								  //RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07,0,0.1);
								  // scaling factor Psi(2S) / JPsi
	Double_t factorSigma = 1.05; // maybe to adjust
	RooRealVar mean_scaling("mean_scaling", "", factorMean);
	RooRealVar sigma_scaling("sigma_scaling", "", factorSigma);
	RooFormulaVar mean_psi2s("mean_psi2s","mean_jpsi*mean_scaling",RooArgSet(mean_jpsi, mean_scaling));
	RooFormulaVar sigma_psi2s("sigma_psi2s", "sigma_jpsi*sigma_scaling", RooArgSet(sigma_jpsi, sigma_scaling));
	//RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07, 0, 0.15);
	RooRealVar alpha_psi2s_L("alpha_psi2s_L","alpha_psi2s_L",alphaL2);
	RooRealVar n_psi2s_L("n_psi2s_L","n_psi2s_L",nL2);
	RooRealVar alpha_psi2s_R("alpha_psi2s_R","alpha_psi2s_R",alphaR2);
	RooRealVar n_psi2s_R("n_psi2s_R","n_psi2s_R",nR2);
	//RooCBShape *psi = new RooCBShape("psi","crystal ball PDF",m,mean_psi,sigma_psi,alpha_psi,n_psi);
	ExtendedCrystalBall *psi2s = new ExtendedCrystalBall("psi2s","crystal ball PDF", m,
														 mean_psi2s, sigma_psi2s, alpha_psi2s_L,
														 n_psi2s_L, alpha_psi2s_R, n_psi2s_R);
	
	//Exponential background
	RooRealVar a1("a1","a1",-2,-10,-0.5);
	RooExponential *bkg = new RooExponential("exp","exp",m,a1);
	
	//RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
	RooRealVar fsigJpsi("fsigJpsi","signalJPsi",1000,0.,1.e4);
	RooRealVar fsigPsi2s("fsigPsi2s","signalPsi2s",100,0.,1.e3);
	RooRealVar fbkg("fbkg","fbkg",500,0.,1.e7);
	
	RooAbsPdf* model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *psi2s, *bkg), RooArgList(fsigJpsi, fsigPsi2s, fbkg), kFALSE);
	
	ws->import(*model);
}

void AddPtModel(RooWorkspace* ws, std::string period, bool useMCtemplate) {
	// Define fitting model for pt for J/Psi signal only
	// Dissociative contribution and exclusive contribution
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	
	//JPsi Exclusive
	RooAbsPdf* jpsiExclusive = nullptr;
	RooRealVar bExc("bExc","bExc",5, 2, 8);
	if (useMCtemplate){
		GetPtHistMC(ws, period, "kIncohJpsiToMu");
		jpsiExclusive = ws->pdf("ptkIncohJpsiToMu");
		jpsiExclusive->SetName("jpsiExc");
	}
	else {
		// H1 formula
		jpsiExclusive = new RooGenericPdf("jpsiExc","exclusive jPsi PDF","(fTrkTrkPt*exp(-bExc*(fTrkTrkPt**2)))",RooArgSet(pt,bExc)) ;
	}
	
	// JPsi Dissociative
	// H1 formula
	RooRealVar bDiss("bDiss","bDiss",0.32, 0, 2);
	RooGenericPdf *jpsiDissociative = new RooGenericPdf("jpsiDiss","dissociative jPsi PDF","(fTrkTrkPt*exp(-bDiss*(fTrkTrkPt**2)))",RooArgSet(pt,bDiss)) ;
	
	
	//RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
	RooRealVar fsigJpsiDiss("fsigJpsiDiss","signalJPsiDiss",100,0.,1.e5);
	RooRealVar fsigJpsiExc("fsigJpsiExc","fsigJpsiExc",100,0.,1.e5);
	
	RooAbsPdf* model = new RooAddPdf("ptfit", "ptfit", RooArgList(*jpsiDissociative, *jpsiExclusive), RooArgList(fsigJpsiDiss, fsigJpsiExc), kFALSE);
	
	ws->import(*model);
}

void DoSPlot(RooWorkspace* ws) {
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
	RooRealVar* mean_psi2s = ws->var("mean_psi2s");
	RooRealVar* sigma_psi2s = ws->var("sigma_psi2s");
	//RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
	//RooRealVar* n_jpsi = ws->var("n_jpsi");
	RooRealVar* alpha_psi2s_L = ws->var("alpha_psi2s_L");
	RooRealVar* n_psi2s_L = ws->var("n_psi2s_L");
	RooRealVar* alpha_psi2s_R = ws->var("alpha_psi2s_R");
	RooRealVar* n_psi2s_R = ws->var("n_psi2s_R");
	
	// Background
	RooRealVar* a1 = ws->var("a1");
	
	mean_jpsi->setConstant();
	sigma_jpsi->setConstant();
	alpha_jpsi_L->setConstant();
	n_jpsi_L->setConstant();
	alpha_jpsi_R->setConstant();
	n_jpsi_R->setConstant();
	
	//mean_psi2s->setConstant();
	//sigma_psi2s->setConstant();
	alpha_psi2s_L->setConstant();
	n_psi2s_L->setConstant();
	alpha_psi2s_R->setConstant();
	n_psi2s_R->setConstant();
	
	a1->setConstant();
	
	
	RooMsgService::instance().setSilentMode(true);
	
	//Now we use the SPlot class to add SWeight to our data set
	// based on our model and our yield variables
	
	RooStats::SPlot * sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *psi2sYield, *bkgYield) );
	
	//Check Sweight properties
	std::cout << "Check SWeights: " << std::endl;
	
	std::cout << std::endl << "Yield of JPsi is "
	<< jPsiYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fsigJpsi") << std::endl;
	
	std::cout << std::endl << "Yield of Psi(2s) is "
	<< psi2sYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fsigPsi2s") << std::endl;
	
	std::cout << std::endl << "Yield of bkg is "
	<< bkgYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fbkg") << std::endl;
	
	
	for (Int_t i=0; i < 10; i ++ )
	{
		std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
		<< " psi(2s) Yield "<< sData->GetSWeight(i, "fsigPsi2s")
		<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
		<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
		<< std::endl;
	}
	
	//import the new data set with Sweight
	
	std::cout << "import new dataset with sWeight" << std::endl;
	ws->import(*data, Rename("dataWithSWeights"));
	
	std::cout << "Sigma = " << sigma_jpsi->getVal() << "\n\n\n\n\n" << std::endl;
}


//____________________________________
void MakePlots(RooWorkspace *ws, std::string rootfilePath, std::string period, bool useCuts, Double_t mMin, Double_t mMax){
	
	bool drawMC = false;
	
	//make some plots
	TCanvas* cv = new TCanvas("splot","splot",800,800) ;
	cv->Divide(2,3);
	
	//get what we need of the workspace
	RooAbsPdf* model = ws->pdf("mfit");
	RooAbsPdf* jPsiModel = ws->pdf("jpsi");
	RooAbsPdf* psi2sModel = ws->pdf("psi2s");
	RooAbsPdf* bkgModel = ws->pdf("exp");
	
	RooAbsPdf* ptmodel = ws->pdf("ptfit");
	RooAbsPdf* jpsiDiss = ws->pdf("jpsiDiss");
	RooAbsPdf* jpsiExc = ws->pdf("jpsiExc");
	
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* pt = ws->var("fTrkTrkPt");
	
	RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
	
	cv->cd(1);
	gPad->SetLeftMargin(0.15) ;
	RooRealVar* jPsiYield = ws->var("fsigJpsi");
	RooRealVar* psi2sYield = ws->var("fsigPsi2s");
	RooPlot* frame = m->frame();
	data->plotOn(frame);
	model->plotOn(frame);
	model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kBlue));
	model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame, Components(*psi2sModel), LineStyle(kDashed), LineColor(kOrange));
	model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen));
	frame->SetTitle("Fit to model to discriminating variable");
	double yMax0 = frame->GetMaximum();
	TText* txt0 = new TText((mMin+mMax+0.1)/2,0.55*yMax0,Form("%.1f J/Psi", jPsiYield->getVal()));
	TText* txtPsi2s = new TText((mMin+mMax+0.1)/2,0.45*yMax0,Form("%.1f Psi(2s)", psi2sYield->getVal()));
	//txt3->SetTextSize(0.05) ;
	frame->addObject(txt0) ;
	frame->addObject(txtPsi2s) ;
	frame->Draw();
	
	
	cv->cd(3);
	gPad->SetLeftMargin(0.15) ;
	data->plotOn(frame, Binning(40));
	model->plotOn(frame, Binning(40));
	RooHist* hpull = frame->pullHist();
	hpull->GetXaxis()->SetRangeUser(mMin, mMax);
	hpull->SetTitle("(data - fit)/#sigma of the mass fit");
	hpull->Draw("");
	cv->cd(5);
	gPad->SetLeftMargin(0.15) ;
	RooHist* hresid = frame->residHist();
	hresid->GetXaxis()->SetRangeUser(mMin, mMax);
	hresid->SetTitle("Residuals (data - fit) of the mass fit");
	hresid->Draw("");
	
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	// create weighted data set for JPsi
	RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fsigJpsi_sw") ;
	// create weighted data set for Psi(2s)
	RooDataSet * dataw_psi2s = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fsigPsi2s_sw") ;
	// create weighted data set for Background
	RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fbkg_sw") ;
	
	
	// Draw pt with weights
	cv->cd(2);
	RooPlot* frameJpsi = pt->frame() ;
	// Fit pt distrib of J/Psi with exclusive / dissociative components
	ptmodel->fitTo(*dataw_jpsi, Extended(), Minos(true), Strategy(2));
	dataw_jpsi->plotOn(frameJpsi, DataError(RooAbsData::SumW2) ) ;
	ptmodel->plotOn(frameJpsi);
	ptmodel->plotOn(frameJpsi, Components(*jpsiDiss), LineStyle(kDashed), LineColor(kRed));
	ptmodel->plotOn(frameJpsi, Components(*jpsiExc), LineStyle(kDashed), LineColor(kGreen));
	
	frameJpsi->SetTitle("pt distribution for J/#Psi with weights");
	
	RooRealVar* jPsiDissYield = ws->var("fsigJpsiDiss");
	RooRealVar* jPsiExclYield = ws->var("fsigJpsiExc");
	double yMax = frameJpsi->GetMaximum();
	TLatex* txt2 = new TLatex(1.5,0.75*yMax,Form("%.1f exclusive J/Psi", jPsiExclYield->getVal()));
	TText* txt = new TText(1.5,0.65*yMax,Form("%.1f dissociative J/Psi", jPsiDissYield->getVal()));
	//txt3->SetTextSize(0.05) ;
	frameJpsi->addObject(txt) ;
	frameJpsi->addObject(txt2) ;
	// Write b parameters on the plot
	RooRealVar* bDiss = ws->var("bDiss");
	RooRealVar* bExc = ws->var("bExc");
	TLatex* txtExc = new TLatex(1.5,0.45*yMax,Form("b_{exc} = %.2f", bExc->getVal()));
	TLatex* txtDiss = new TLatex(1.5,0.35*yMax,Form("b_{diss} = %.2f", bDiss->getVal()));
	frameJpsi->addObject(txtExc) ;
	frameJpsi->addObject(txtDiss) ;
	
	frameJpsi->Draw() ;
	
	cv->cd(4);
	RooPlot* framePsi2s = pt->frame() ;
	dataw_psi2s->plotOn(framePsi2s, DataError(RooAbsData::SumW2) ) ;
	framePsi2s->SetTitle("pt distribution for #Psi(2s) with weights");
	framePsi2s->Draw() ;
	
	
	// Plot isolation for QCD component.
	// Eg. plot all events weighted by the sWeight for the QCD component.
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	cv->cd(6);
	RooPlot* frameBkg = pt->frame() ;
	dataw_bkg->plotOn(frameBkg,DataError(RooAbsData::SumW2) ) ;
	
	frameBkg->SetTitle("pt distribution for background with weights");
	
	// Add MC templates for (gammma gamma to mu mu) for comparison
	if (drawMC) {
		GetPtHistMC(ws, period, "kTwoGammaToMuLow");
		GetPtHistMC(ws, period, "kTwoGammaToMuMedium");
		RooAbsPdf* ptMuLow = ws->pdf("ptkTwoGammaToMuLow");
		RooAbsPdf* ptMuMedium = ws->pdf("ptkTwoGammaToMuMedium");
		RooRealVar fsigMuLow("fsigMuLow","signalMuLow",100,0.,2.e3);
		RooRealVar fsigMuMedium("fsigMuMedium","signalMuMedium",100,0.,2.e3);
		
		RooAbsPdf* modelMuLow = new RooAddPdf("ptfitMuLow", "ptfitMuLow", RooArgList(*ptMuLow), RooArgList(fsigMuLow), kFALSE);
		RooAbsPdf* modelMuMedium = new RooAddPdf("ptfitMuMedium", "ptfitMuMedium", RooArgList(*ptMuMedium), RooArgList(fsigMuMedium), kFALSE);
		modelMuLow->fitTo(*dataw_bkg, Extended());
		modelMuMedium->fitTo(*dataw_bkg, Extended());
		modelMuLow->plotOn(frameBkg, LineStyle(kDashed), LineColor(8));
		modelMuMedium->plotOn(frameBkg, LineStyle(kDashed), LineColor(9));
	}
	frameBkg->Draw() ;
	frameBkg->Print("v");
	//return;
	
	TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.9);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	if (drawMC) {
		legend->AddEntry(frameBkg->findObject("ptfitMuLow_Norm[fTrkTrkPt]"), "MuLow","L");
		legend->AddEntry(frameBkg->findObject("ptfitMuMedium_Norm[fTrkTrkPt]"), "MuMedium","L");
	}
	legend->Draw("same");
	
	// Save plots
	std::string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	cv->SaveAs(Form("Plots/Splot-%s%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax));
	
	/*
	 // Look at V0CCounts
	 cv->cd(2);
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
	TH1* hPt = dataw_jpsi->createHistogram("fTrkTrkPt",50) ;
	TH1* hPtDissociative = jpsiDiss->createHistogram("fTrkTrkPt",50) ;
	TH1* hPtExclusive = jpsiExc->createHistogram("fTrkTrkPt",50) ;
	TH1* hPtBackground = dataw_bkg->createHistogram("fTrkTrkPt",50) ;
	
	
	TFile* fNewTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"RECREATE");
	hPtExclusive->SetName("hPtExclusive");
	hPtExclusive->Write();
	hPtDissociative->SetName("hPtDissociative");
	hPtDissociative->Write();
	hPtBackground->SetName("hPtBackground");
	hPtBackground->Write();
	
	fNewTemplates->Close();
	
}


void Splot(std::string rootfilePath = "", std::vector<std::string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2., Double_t mMax = 3.5, Double_t ptMin = 0., Double_t ptMax = 3.5, bool useCuts = true) {
	
	gROOT->ProcessLine(".L ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./ExtendedCrystalBall_cxx.so") ;
	
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
		std::string period = periods[k];
		
		std::list<TCut> mCutList = DefineCuts(period);
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
		AddModel(wspace, period);
		AddPtModel(wspace, period, false);	// don't use MC template but H1 function
		wspace->Print();
		
		DoSPlot(wspace);
		MakePlots(wspace, rootfilePath, period, useCuts, mMin, mMax);

		//cleanup
		delete wspace;
	}
	
	
}


