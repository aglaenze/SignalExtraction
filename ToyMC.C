#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>
#include <TROOT.h>

#include "TSystem.h"
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
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBasket.h"

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

using namespace RooFit;

Double_t mMin = 2.6, mMax = 3.8;
Double_t ptMin = 0., ptMax = 3.;

Int_t nEventsSignal = 5000, nEventsBkg = 3000;

void AddModel(RooWorkspace* ws) {
	// Define model
	
	//Crystal ball for the J/psi

	RooRealVar m = *ws->var("m");
	
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3.1,2.9,3.2);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.07,0,0.15);
	RooRealVar alpha_jpsi("alpha_jpsi","alpha_jpsi",1, 0, 5);
	RooRealVar n_jpsi("n_jpsi","n_jpsi",5, 0, 10);
	RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
	
	//Exponential background
	RooRealVar a1("a1","a1",-1.,-5,-0.);
	RooExponential *bkg = new RooExponential("exp","exp",m,a1);
	
	//RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
	RooRealVar fsigJPsi("fsigJPsi","signalJPsi",5.e3,1.e3,1.e5);
	RooRealVar fbkg("fbkg","fbkg",5.e3,1.e3,1.e5);
	
	RooAbsPdf* model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *bkg), RooArgList(fsigJPsi, fbkg), kFALSE);
	
	ws->import(*model);
}

void GetPtHistMC(RooWorkspace* ws, std::string rootfilePathMC, std::string period){
	TFile *fSimu = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"READ");
	TTree* tSignal = (TTree*)fSimu->Get("tSignal");
	TTree* tBkg = (TTree*)fSimu->Get("tBkg");
	
	RooRealVar ptSignal("pt","Dimuon p_{T} (GeV/c)",0,4);
	RooRealVar ptBkg("pt","Dimuon p_{T} (GeV/c)",0,4);
	RooDataSet* dataSignal = new RooDataSet("data","data",RooArgSet(ptSignal),Import(*tSignal));
	RooDataSet* dataBkg = new RooDataSet("data","data",RooArgSet(ptBkg),Import(*tBkg));
	
	TH1F* histSignalPt = new TH1F("hPtSignal", "hPtSignal", 100, 0, 4);
	TH1F* histBkgPt = new TH1F("hPtBkg", "hPtBkg", 100, 0, 4);
	tSignal->Draw("pt>>hPtSignal");
	tBkg->Draw("pt>>hPtBkg");
	RooDataHist* ptSignalHist = new RooDataHist("ptSignalHist","ptSignalHist", RooArgList(ptSignal),histSignalPt);
	RooDataHist* ptBkgHist = new RooDataHist("ptBkgHist","ptBkgHist", RooArgList(ptBkg),histBkgPt);
	RooHistPdf* ptSignalPdf = new RooHistPdf("ptSignalPdf", "ptSignalPdf", ptSignal, *ptSignalHist);
	RooHistPdf* ptBkgPdf = new RooHistPdf("ptBkgPdf", "ptBkgPdf", ptBkg, *ptBkgHist);
	
	ws->import(*ptSignalPdf);
	ws->import(*ptBkgPdf);
}

void DoSPlot(RooWorkspace* ws) {
	RooAbsPdf* model = ws->pdf("mfit");
	RooRealVar* jPsiYield = ws->var("fsigJPsi");
	RooRealVar* bkgYield = ws->var("fbkg");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	model->fitTo(*data, Extended(), Minos(true), Strategy(2));
	
	//The sPlot technique requires that we fix the parameters
	// of the model that are not yields after doing the fit
	RooRealVar* mean_jpsi = ws->var("mean_jpsi");
	RooRealVar* sigma_jpsi = ws->var("sigma_jpsi");
	RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
	RooRealVar* n_jpsi = ws->var("n_jpsi");
	
	RooRealVar* a1 = ws->var("a1");
	
	mean_jpsi->setConstant();
	sigma_jpsi->setConstant();
	alpha_jpsi->setConstant();
	n_jpsi->setConstant();
	a1->setConstant();
	
	RooMsgService::instance().setSilentMode(true);
	
	//Now we use the SPlot class to add SWeight to our data set
	// based on our model and our yield variables
	
	RooStats::SPlot * sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield,*bkgYield) );
	
	//Check Sweight properties
	std::cout << "Check SWeights: " << std::endl;
	
	std::cout << std::endl << "Yield of JPsi is "
	<< jPsiYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fsigJPsi") << std::endl;
	
	std::cout << std::endl << "Yield of bkg is "
	<< bkgYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fbkg") << std::endl;
	
	for (Int_t i=0; i < 10; i ++ )
	{
		std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJPsi")
		<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
		<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
		<< std::endl;
	}
	
	//import the new data set with Sweight
	
	std::cout << "import new dataset with sWeight" << std::endl;
	ws->import(*data, Rename("dataWithSWeights"));
	
}


void MakePlots(RooWorkspace *ws, std::string period, bool drawPulls) {
	
	//make some plots
	TCanvas* cv = new TCanvas("splot","splot",900,900) ;
	if (drawPulls) cv->Divide(3, 3);
	else {cv->Divide(3, 2); cv->SetCanvasSize(900, 600);}
	
	//get what we need of the workspace
	RooAbsPdf* model = ws->pdf("mfit");
	RooAbsPdf* jPsiModel = ws->pdf("jpsi");
	RooAbsPdf* bkgModel = ws->pdf("exp");
	
	RooAbsPdf* ptSignalPdf = ws->pdf("ptSignalPdf");
	RooAbsPdf* ptBkgPdf = ws->pdf("ptBkgPdf");
	
	RooRealVar* m = ws->var("m");
	RooRealVar* pt = ws->var("pt");
	
	RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

	// First draw mass fit
	cv->cd(1);
	gPad->SetLeftMargin(0.15);
	RooRealVar* jPsiYield = ws->var("fsigJPsi");
	RooRealVar* bkgYield = ws->var("fbkg");
	RooPlot* frame = m->frame();
	data->plotOn(frame);
	model->plotOn(frame);
	model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen));
	frame->SetTitle("Fit to model to discriminating variable");
	double yMax0 = frame->GetMaximum();
	TText* txt0 = new TText((mMax+mMin+0.1)/2,0.6*yMax0,Form("%.1f J/Psi", jPsiYield->getVal()));
	TText* txt1 = new TText((mMax+mMin+0.1)/2,0.5*yMax0,Form("%.1f Bkg", bkgYield->getVal()));
	//txt3->SetTextSize(0.05) ;
	frame->addObject(txt0) ;
	frame->addObject(txt1) ;
	frame->Draw();
	
	
	// Get weighted data
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	RooDataSet *dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fsigJPsi_sw") ;
	RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fbkg_sw") ;

	// Extract pt templates to compare them with reconstructed data
	RooRealVar yieldSignal("yieldSignal","yieldSignal",nEventsSignal);
	RooRealVar yieldBkg("yieldBkg","yieldBkg",nEventsBkg);
	RooAbsPdf* ptSignalModel = new RooAddPdf("ptSignalMC", "ptSignalMC", RooArgList(*ptSignalPdf), RooArgList(yieldSignal), kFALSE);
	RooAbsPdf* ptBkgModel = new RooAddPdf("ptBkgMC", "ptBkgMC", RooArgList(*ptBkgPdf), RooArgList(yieldBkg), kFALSE);
	//ptSignalModel->fitTo(*dataw_jpsi, Extended(), Minos(true), Strategy(2));
	//ptBkgModel->fitTo(*dataw_bkg, Extended(), Minos(true), Strategy(2));

	// Draw pt with weights
	// For signal
	cv->cd(2);
	gPad->SetLeftMargin(0.15);
	RooPlot* frame2 = pt->frame() ;
	dataw_jpsi->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
	ptSignalModel->plotOn(frame2, LineStyle(kDashed), LineColor(kRed));
	frame2->SetTitle("p_{T} distribution for J/#Psi with weights");
	frame2->Draw() ;
	TLegend* lgdJpsi = new TLegend(0.3, 0.7, 0.9, 0.9);
	lgdJpsi->AddEntry(frame2->findObject("h_dataWithSWeights"), "reconstructed signal", "LP");
	lgdJpsi->AddEntry(frame2->findObject("ptSignalMC_Norm[pt]"), "original distribution", "L");
	lgdJpsi->Draw();
	
	// And for background
	cv->cd(3);
	gPad->SetLeftMargin(0.15);
	RooPlot* frame3 = pt->frame() ;
	dataw_bkg->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;
	ptBkgModel->plotOn(frame3, LineStyle(kDashed), LineColor(kRed));
	frame3->SetTitle("p_{T} distribution for #gamma#gamma #rightarrow #mu#mu with weights");
	frame3->Draw() ;
	frame3->Print("v");
	TLegend* lgdBkg = new TLegend(0.3, 0.7, 0.9, 0.9);
	lgdBkg->AddEntry(frame3->findObject("h_dataWithSWeights"), "reconstructed signal", "LP");
	lgdBkg->AddEntry(frame3->findObject("ptBkgMC_Norm[pt]"), "original distribution", "L");
	lgdBkg->Draw();
	//return;
	
	// Then finally, draw quality plots
	// For the mass fit
	cv->cd(4);
	gPad->SetLeftMargin(0.15);
	data->plotOn(frame, Binning(40));
	model->plotOn(frame, Binning(40));
	RooHist* hpull = frame->pullHist();
	hpull->GetXaxis()->SetRangeUser(mMin, mMax);
	hpull->SetTitle("(data - fit)/#sigma of the mass fit");
	hpull->Draw("");
	cv->cd(7);
	gPad->SetLeftMargin(0.15);
	RooHist* hresid = frame->residHist();
	hresid->GetXaxis()->SetRangeUser(mMin, mMax);
	hresid->SetTitle("Residuals (data - fit) of the mass fit");
	hresid->Draw("");
	
	// Draw quality histograms
	// First (data-fit)/sigma
	// For signal pt
	cv->cd(5);
	gPad->SetLeftMargin(0.15);
	dataw_jpsi->plotOn(frame2, Binning(40));
	ptSignalModel->plotOn(frame2, Binning(40));
	RooHist* hpull2 = frame2->pullHist();
	hpull2->GetXaxis()->SetRangeUser(ptMin, ptMax);
	hpull2->SetTitle("(reconstructed - original)/#sigma of J/#Psi p_{T}");
	hpull2->Draw("");
	
	// For background pt
	cv->cd(6);
	gPad->SetLeftMargin(0.15);
	dataw_bkg->plotOn(frame3, Binning(40));
	ptBkgModel->plotOn(frame3, Binning(40));
	RooHist* hpull3 = frame3->pullHist();
	hpull3->GetXaxis()->SetRangeUser(ptMin, ptMax);
	hpull3->SetTitle("(reconstructed - original)/#sigma of bkg p_{T}");
	hpull3->Draw("");
	
	// Draw pull histograms = data-fit
	if (drawPulls) {
	cv->cd(8);
	gPad->SetLeftMargin(0.15);
	RooHist* hresid2 = frame2->residHist();
	hresid2->GetXaxis()->SetRangeUser(ptMin, ptMax);
	hresid2->SetTitle("Residuals (reconstructed - original) of J/#Psi p_{T}");
	hresid2->Draw("");
	
	cv->cd(9);
	gPad->SetLeftMargin(0.15);
	RooHist* hresid3 = frame3->residHist();
	hresid3->GetXaxis()->SetRangeUser(ptMin, ptMax);
	hresid3->SetTitle("Residuals (reconstructed - original) of bkg p_{T}");
	hresid3->Draw("");
	}

	cv->SaveAs(Form("Plots/toyMC-splot-%s.pdf", period.c_str()));

}



bool CreateTrees(std::string rootfilePathMC, std::string period, bool draw = false) {
	
	gStyle->SetOptStat(0);
	
	// 1ere étape : générer des "fausses" données en utilisant les MC
	
	// MC files
	TFile *fSignal = new TFile(Form("%s/AnalysisResults_%s_MC_kIncohJpsiToMu.root", rootfilePathMC.c_str(), period.c_str()),"READ");
	TTree* fSignalTree = (TTree*)fSignal->Get("MyTask/fAnaTree");
	
	TFile *fBkg = new TFile(Form("%s/AnalysisResults_%s_MC_kTwoGammaToMuLow.root", rootfilePathMC.c_str(), period.c_str()),"READ");
	TTree* fBkgTree = (TTree*)fBkg->Get("MyTask/fAnaTree");
	
	Double_t mSignal = 0, mBkg = 0, ptSignal = 0, ptBkg = 0;
	fSignalTree->SetBranchAddress("fTrkTrkM", &mSignal);
	fSignalTree->SetBranchAddress("fTrkTrkPt", &ptSignal);
	fBkgTree->SetBranchAddress("fTrkTrkM", &mBkg);
	fBkgTree->SetBranchAddress("fTrkTrkPt", &ptBkg);
	fBkgTree->SetBranchStatus("*",0);	// Select the branches to look at
	fBkgTree->SetBranchStatus("fTrkTrkM",1);	// Select the branches to look at
	fBkgTree->SetBranchStatus("fTrkTrkPt",1);	// Select the branches to look at
	
	const int nSignal = fSignalTree->GetEntries();
	const int nBkg = fBkgTree->GetEntries();
	
	
	// new file
	TFile *fAna = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"RECREATE");
	
	// new TTrees
	TTree* tSignal = new TTree("tSignal", "tSignal");
	tSignal->Branch("m", &mSignal, "m/D");
	tSignal->Branch("pt", &ptSignal, "pt/D");
	TTree* tBkg = new TTree("tBkg", "tBkg");
	tBkg->Branch("m", &mBkg, "m/D");
	tBkg->Branch("pt", &ptBkg, "pt/D");
	Double_t mMixed = 0, ptMixed = 0;
	TTree* tMixed = new TTree("tMixed", "tMixed");
	tMixed->Branch("m", &mMixed, "m/D");
	tMixed->Branch("pt", &ptMixed, "pt/D");
	
	int j = 0;
	
	//for (int i = 0; i<nSignal; i++) {
	for (int i = 0; i<nEventsSignal; i++) {
		j++;
		fSignalTree->GetEntry(j);
		if (j == nSignal) {std::cout << "not enough stats! bye" << std::endl; return false;}
		if (mSignal < mMin || mSignal > mMax || ptSignal < ptMin || ptSignal > ptMax) {
			std::cout << mSignal << std::endl;
			std::cout << "event (signal) ignored because not in the right range of pt or m" << std::endl;
			i--;
			continue;}
		mMixed = mSignal;
		ptMixed = ptSignal;
		
		tSignal->Fill();
		tMixed->Fill();
		/*
		 histM->Fill(mSignal);
		 histPt->Fill(ptSignal);
		 histSignalM->Fill(mSignal);
		 histSignalPt->Fill(ptSignal);
		 */
	}
	
	j=0;
	//for (int i = 0; i<nBkg; i++) {
	for (int i = 0; i<nEventsBkg; i++) {
		j++;
		fBkgTree->GetEntry(j);
		if (j == nBkg) {std::cout << "not enough stats! bye" << std::endl; return false;}
		if (mBkg < mMin || mBkg > mMax || ptBkg < ptMin || ptBkg > ptMax) {
			std::cout << mBkg << std::endl;
			std::cout << "event (background) ignored because not in the right range of pt or m" << std::endl;
			i--;
			continue;}
		mMixed = mBkg;
		ptMixed = ptBkg;
		/*
		 histM->Fill(mBkg);
		 histPt->Fill(ptBkg);
		 histBkgM->Fill(mBkg);
		 histBkgPt->Fill(ptBkg);
		 */
		tBkg->Fill();
		tMixed->Fill();
	}
	
	if (draw) {
		TCanvas* cv = new TCanvas("ToyMC","Toy MC",600,300) ;
		//TCanvas* cv = new TCanvas("2Dplot","2D fit",800,300) ;
		TH1F* histM = new TH1F("histM", "M distribution", 100, mMin, mMax);
		TH1F* histPt = new TH1F("histPt", "Pt distribution", 100, ptMin, ptMax);
		
		TH1F* histSignalM = new TH1F("histSignalM", "M Signal distribution", 100, mMin, mMax);
		TH1F* histSignalPt = new TH1F("histSignalPt", "Pt Signal distribution", 100, ptMin, ptMax);
		TH1F* histBkgM = new TH1F("histBkgM", "M Bkg distribution", 100, mMin, mMax);
		TH1F* histBkgPt = new TH1F("histBkgPt", "Pt Bkg distribution", 100, ptMin, ptMax);
		// m distributions
		cv->Divide(3,2) ;
		cv->cd(1);
		tSignal->Draw("m>>histSignalM");
		histSignalM->Draw();
		cv->cd(2);
		tBkg->Draw("m>>histBkgM");
		histBkgM->Draw();
		cv->cd(3);
		tMixed->Draw("m>>histM");
		histM->Draw();
		
		// pt distributions
		cv->cd(4);
		tSignal->Draw("pt>>histSignalPt");
		histSignalPt->Draw();
		cv->cd(5);
		tBkg->Draw("pt>>histBkgPt");
		histBkgPt->Draw();
		cv->cd(6);
		tMixed->Draw("pt>>histPt");
		histPt->Draw();
		
		cv->SaveAs(Form("Plots/ToyMC-%s.pdf", period.c_str()));
	}
	
	// Finally, write TTrees
	fAna->cd();
	tSignal->Write();
	tBkg->Write();
	tMixed->Write();
	fAna->Close();
	return true;
}

// 2e étape : faire un sPlot pour extraire ma distribution en pt
void ToyMC(std::string rootfilePathMC, bool drawPulls) {
	
	std::vector <std::string> periods = {"LHC16r", "LHC16s"};
	const int nPeriods = periods.size();
	
	for (int k = 0; k<nPeriods; k++) {
		std::string period = periods[k];
		bool newTrees = CreateTrees(rootfilePathMC, period, true);	// false = don't draw
		if (!newTrees) {std::cout << "Trees not created! bye" << std::endl; return;}
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()), "READ");
		// Connect to the tree
		TTree* t = (TTree*)fAna->Get("tMixed");
		
		// New working space
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		// Import data
		RooRealVar m("m","M_{#mu#mu} (GeV/c2)", mMin, mMax);
		RooRealVar pt("pt","Dimuon p_{T} (GeV/c)", ptMin, ptMax);
		RooArgSet variables(m, pt);
		RooDataSet* data = new RooDataSet("data","data",variables,Import(*t));
		wspace->import(*data, Rename("data"));
		data->Print();

		// Add model for m distribution
		AddModel(wspace);
		
		wspace->Print();
		
		DoSPlot(wspace);
		
		GetPtHistMC(wspace, rootfilePathMC, period);
		MakePlots(wspace, period, drawPulls);
		
		//cleanup
		delete wspace;
	}
}
