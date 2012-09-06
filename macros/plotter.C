#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TRandom.h>
#include "../src/MultiPlot.cxx"
#include "../src/MultiFitPlot.cxx"

#define DEBUG false

// TODO - POLYMORPHYSM

const Double_t particleMasses[3] = { // in GeV / c^2
	0.493677,// kaon
	0.1349766, // pion
	0.938272029 // proton
};

const int numberOfCentralities = 11;

void createLcmsPlot(std::string *graphNames, std::string *prefixes, const char *fileName, MultiFitPlot &Rinv, MultiFitPlot &Rlcms);
void createPrfPlotAllInONe(MultiFitPlot Rinv[], MultiFitPlot Rlcms[], int nCentralities);
void createPrfPlots(MultiFitPlot Rinv[], MultiFitPlot Rlcms[], int nCentralities);
void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant = kFALSE);

int plotter()
{	
	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.06, "xyz");
	gStyle->SetPadTopMargin(0.023);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetPadRightMargin(0.001);

	MultiFitPlot::defaultFunction  = new TF1("fit","[0]*TMath::Power(x,-[1])", 0.1, 1.5);
	MultiFitPlot::defaultFunction->SetParameter(0,1);
	MultiFitPlot::defaultFunction->SetParameter(1,1);

	std::string prefixes[graphCount];
	std::string graphNames[graphCount];

	MultiFitPlot Rinv[numberOfCentralities];
	MultiFitPlot Rlcms[numberOfCentralities];

	prefixes[0] = "b2/kk";
	prefixes[1] = "b2/pipi";
	prefixes[2] = "b2/pp";
	graphNames[0] = "K-K b = 2 fm";
	graphNames[1] = "\\pi-\\pi b = 2 fm";
	graphNames[2] = "p-p b = 2 fm";
	createLcmsPlot(graphNames, prefixes, "b2", Rinv[0], Rlcms[0]);

	prefixes[0] = "b3/kk";
	prefixes[1] = "b3/pipi";
	prefixes[2] = "b3/pp";
	graphNames[0] = "K-K b = 3 fm";
	graphNames[1] = "\\pi-\\pi b = 3 fm";
	graphNames[2] = "p-p b = 3 fm";
	createLcmsPlot(graphNames, prefixes, "b3", Rinv[1], Rlcms[1]);

	prefixes[0] = "b5/kk";
	prefixes[1] = "b5/pipi";
	prefixes[2] = "b5/pp";
	graphNames[0] = "K-K b = 5 fm";
	graphNames[1] = "\\pi-\\pi b = 5 fm";
	graphNames[2] = "p-p b = 5 fm";
	createLcmsPlot(graphNames, prefixes, "b5", Rinv[2], Rlcms[2]);

	// old EPOS results
	// prefixes[0] = "bb3m6/kk";
	// prefixes[1] = "bb3m6/pipi";
	// prefixes[2] = "bb3m6/pp";
	// graphNames[0] = "K-K EPOS";
	// graphNames[1] = "\\pi-\\pi EPOS";
	// graphNames[2] = "p-p EPOS";
	// createLcmsPlot(graphNames, prefixes, "bb3m6_EPOS", Rinv[3], Rlcms[3]);

	prefixes[0] = "lhc0005/kk";
	prefixes[1] = "lhc0005/pipi";
	prefixes[2] = "lhc0005/pp";
	graphNames[0] = "K-K 0-5%";
	graphNames[1] = "\\pi-\\pi 0-5%";
	graphNames[2] = "p-p 0-5%";
	createLcmsPlot(graphNames, prefixes, "lhc0005_hHKM", Rinv[3], Rlcms[3]);

	prefixes[0] = "lhc1020/kk";
	prefixes[1] = "lhc1020/pipi";
	prefixes[2] = "lhc1020/pp";
	graphNames[0] = "K-K 10-20%";
	graphNames[1] = "\\pi-\\pi 10-20%";
	graphNames[2] = "p-p 10-20%";
	createLcmsPlot(graphNames, prefixes, "lhc1020_hHKM", Rinv[4], Rlcms[4]);

	prefixes[0] = "lhc2030/kk";
	prefixes[1] = "lhc2030/pipi";
	prefixes[2] = "lhc2030/pp";
	graphNames[0] = "K-K 20-30%";
	graphNames[1] = "\\pi-\\pi 20-30%";
	graphNames[2] = "p-p 20-30%";
	createLcmsPlot(graphNames, prefixes, "lhc2030_hHKM", Rinv[5], Rlcms[5]);

	prefixes[0] = "lhc3040/kk";
	prefixes[1] = "lhc3040/pipi";
	prefixes[2] = "lhc3040/pp";
	graphNames[0] = "K-K 30-40%";
	graphNames[1] = "\\pi-\\pi 30-40%";
	graphNames[2] = "p-p 30-40%";
	createLcmsPlot(graphNames, prefixes, "lhc3040_hHKM", Rinv[6], Rlcms[6]);

	prefixes[0] = "epos_0005/kk";
	prefixes[1] = "epos_0005/pipi";
	prefixes[2] = "epos_0005/pp";
	graphNames[0] = "K-K 0-5%";
	graphNames[1] = "\\pi-\\pi 0-5%";
	graphNames[2] = "p-p 0-5%";
	createLcmsPlot(graphNames, prefixes, "epos_0005", Rinv[7], Rlcms[7]);

	prefixes[0] = "epos_1020/kk";
	prefixes[1] = "epos_1020/pipi";
	prefixes[2] = "epos_1020/pp";
	graphNames[0] = "K-K 10-20%";
	graphNames[1] = "\\pi-\\pi 10-20%";
	graphNames[2] = "p-p 10-20%";
	createLcmsPlot(graphNames, prefixes, "epos_1020", Rinv[8], Rlcms[8]);

	prefixes[0] = "epos_2030/kk";
	prefixes[1] = "epos_2030/pipi";
	prefixes[2] = "epos_2030/pp";
	graphNames[0] = "K-K 20-30%";
	graphNames[1] = "\\pi-\\pi 20-30%";
	graphNames[2] = "p-p 20-30%";
	createLcmsPlot(graphNames, prefixes, "epos_2030", Rinv[9], Rlcms[9]);

	prefixes[0] = "epos_3040/kk";
	prefixes[1] = "epos_3040/pipi";
	prefixes[2] = "epos_3040/pp";
	graphNames[0] = "K-K 30-40%";
	graphNames[1] = "\\pi-\\pi 30-40%";
	graphNames[2] = "p-p 30-40%";
	createLcmsPlot(graphNames, prefixes, "epos_3040", Rinv[10], Rlcms[10]);

	createPrfPlots(Rinv, Rlcms, numberOfCentralities);
	return 0;
}

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant)
{
	char buffer[256];
	std::ifstream infile(fileName.c_str(), std::ifstream::in);
	double kT, kTmin, kTmax, R, dR, mT, gamma;
	Int_t i;
	if(DEBUG) std::cout << std::endl << fileName << std::endl << "mT\tR\tdR" << std::endl;
	while(infile.good())
	{
		for(i=0; i < 256; ++i)
			buffer[i] = '\0';
    	infile >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	i = graph->GetN();

    	std::stringstream(buffer) >> kTmin;
    	infile >> kTmax;
    	infile >> R;
    	infile >> dR;
    	kT = (kTmax + kTmin)/2;
		mT = TMath::Sqrt(
			TMath::Power(kT,2)
			+ TMath::Power(particleMasses[iParticle],2));
		if(isInvariant)
		{
			gamma = 1/TMath::Sqrt(1-(kT/mT));
			R /= TMath::Sqrt( (TMath::Sqrt(gamma) + 2) / 3 );
		}
		if(DEBUG) std::cout << mT << "\t" << R << "\t+/- " << dR << std::endl;
		graph->SetPoint(i, mT, R);
		graph->SetPointError(i, 0, dR);
	}
}

void createLcmsPlot(std::string *graphNames, std::string *prefixes, const char *fileName, MultiFitPlot &Rinv, MultiFitPlot &Rlcms)
{
	TCanvas *canvas = new TCanvas("canvas", "R_LCMS", 900, 800);
	canvas->Divide(2,2);

	MultiFitPlot Rout(";m_{T} [GeV/c^{2}];R_{out} [fm]");

	Rout.graphNames[0] = graphNames[0];
	Rout.graphNames[1] = graphNames[1];
	Rout.graphNames[2] = graphNames[2];

	Rout.theme[0].markerColor = 1;
	Rout.theme[0].lineColor = 1;
	Rout.theme[0].markerSize = 1.5;
	Rout.theme[0].markerStyle = 20;

	Rout.theme[1].markerColor = 4;
	Rout.theme[1].lineColor = 4;
	Rout.theme[1].markerSize = 1.4;
	Rout.theme[1].markerStyle = 21;

	Rout.theme[2].markerColor = 2;
	Rout.theme[2].lineColor = 2;
	Rout.theme[2].markerSize = 1.5;
	Rout.theme[2].markerStyle = 22;

	MultiFitPlot Rside(Rout);
	MultiFitPlot Rlong(Rout);
	Rlcms = MultiFitPlot(Rout);
	Rinv = MultiFitPlot(Rout);

	Rside.labels = ";m_{T} [GeV/c^{2}];R_{side} [fm]";
	Rlong.labels = ";m_{T} [GeV/c^{2}];R_{long} [fm]";
	Rlcms.labels = ";m_{T} [GeV/c^{2}];R_{LCMS} [fm]";
	Rinv.labels = ";m_{T} [GeV/c^{2}];R_{inv}/[(\\sqrt{\\gamma}+2)/3]^{1/2} [fm]";

	for( int j = 0; j < graphCount; ++j)
	{
		if(DEBUG) std::cout << std::endl << "  ###  " << prefixes[j] << "  ###   " << std::endl;
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rout.out"), Rout.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rside.out"), Rside.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rlong.out"), Rlong.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rinv.out"), Rinv.graphs[j], j, kTRUE);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rlcms.out"), Rlcms.graphs[j], j);
	}

	Rout.Fit();
	Rside.Fit();
	Rlong.Fit();
	Rlcms.Fit();
	Rinv.Fit();

	canvas->cd(1);
	Rout.Draw();

	canvas->cd(2);
	Rside.Draw();

	canvas->cd(3);
	Rlong.Draw();

	canvas->cd(4);
	Rlcms.Draw();


	canvas->SaveAs((std::string("output/")+fileName+std::string(".png")).c_str());
	delete canvas;

	// normalized
	Rout.labels = ";m_{T} [GeV/c^{2}];R_{out} / R_{out}^{FIT}";
	Rside.labels = ";m_{T} [GeV/c^{2}];R_{side} / R_{side}^{FIT}";
	Rlong.labels = ";m_{T} [GeV/c^{2}];R_{long} / R_{long}^{FIT}";
	Rlcms.labels = ";m_{T} [GeV/c^{2}];R_{LCMS} / R_{LCMS}^{FIT}";

	canvas = new TCanvas("canvas", "R_LCMS", 900, 800);
	canvas->SetGrid();
	gStyle->SetPadGridX(true);
	gStyle->SetPadGridY(true);
	canvas->Divide(2,2);
	canvas->cd(1);
	Rout.GetNormalizedPlot().Draw();
	canvas->cd(2);
	Rside.GetNormalizedPlot().Draw();
	canvas->cd(3);
	Rlong.GetNormalizedPlot().Draw();
	canvas->cd(4);
	Rlcms.GetNormalizedPlot().Draw();
	canvas->SaveAs((std::string("output/")+fileName+std::string("_div.png")).c_str());
	delete canvas;
	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
}

void createPrfPlotAllInONe(MultiFitPlot *Rinv, MultiFitPlot *Rlcms, int nCentralities)
{
	TCanvas *canvas = new TCanvas("canvas", "R_LCMS", 2500, 600);
	canvas->Divide(nCentralities, 2);

	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1].Draw();
		Rinv[i-1].Fit();
		canvas->cd(i+nCentralities);
		Rlcms[i-1].Draw();
	}
	canvas->SaveAs("output/all.png");
	delete canvas;
}

void createPrfPlots(MultiFitPlot *Rinv, MultiFitPlot *Rlcms, int nCentralities)
{
	TCanvas *canvas = 0;

	// therminator
	nCentralities = 3;
	canvas = new TCanvas("canvas", "R_LCMS", 900, 600);
	canvas->Divide(nCentralities, 2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1].Draw();
		Rinv[i-1].Fit();
		canvas->cd(i+nCentralities);
		Rlcms[i-1].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} [fm]";
		Rlcms[i-1].Draw();
	}
	canvas->SaveAs("output/all_therminator.png");
	delete canvas;

	// hHKM
	nCentralities = 4;
	canvas = new TCanvas("canvas", "R_LCMS", 1200, 600);
	canvas->Divide(nCentralities, 2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1+3].Draw();
		Rinv[i-1+3].Fit();
		canvas->cd(i+nCentralities);
		Rlcms[i-1+3].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} [fm]";
		Rlcms[i-1+3].Draw();
	}
	canvas->SaveAs("output/all_hHKM.png");
	delete canvas;

	// EPOS
	nCentralities = 4;
	canvas = new TCanvas("canvas", "R_LCMS", 1200, 600);
	canvas->Divide(nCentralities, 2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1+7].Draw();
		Rinv[i-1+7].Fit();
		canvas->cd(i+nCentralities);
		Rlcms[i-1+7].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} [fm]";
		Rlcms[i-1+7].Draw();
	}
	canvas->SaveAs("output/all_EPOS.png");
	delete canvas;


	// normalized
	gStyle->SetPadGridX(true);
	gStyle->SetPadGridY(true);
	// therminator
	nCentralities = 3;
	canvas = new TCanvas("canvas", "R_LCMS", 900, 600);
	canvas->Divide(nCentralities, 2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1].labels = ";m_{T} [GeV/c^{2}];R_{inv}/ R_{inv}^{FIT}";
		Rinv[i-1].GetNormalizedPlot().Draw();
		canvas->cd(i+nCentralities);
		Rlcms[i-1].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} / R_{LCMS}^{FIT}";
		Rlcms[i-1].GetNormalizedPlot().Draw();
	}
	canvas->SaveAs("output/all_therminator_div.png");
	delete canvas;


	// hHKM
	nCentralities = 4;
	canvas = new TCanvas("canvas", "R_LCMS", 1200, 600);
	canvas->Divide(nCentralities, 2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1+3].labels = ";m_{T} [GeV/c^{2}];R_{inv}/ R_{inv}^{FIT}";
		Rinv[i-1+3].GetNormalizedPlot().Draw();
		canvas->cd(i+nCentralities);
		Rlcms[i-1+3].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} / R_{LCMS}^{FIT}";
		Rlcms[i-1+3].GetNormalizedPlot().Draw();
	}
	canvas->SaveAs("output/all_hHKM_div.png");
	delete canvas;

	// EPOS
	nCentralities = 4;
	canvas = new TCanvas("canvas", "R_LCMS", 1200, 600);
	canvas->Divide(nCentralities,2);
	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1+7].labels = ";m_{T} [GeV/c^{2}];R_{inv}/ R_{inv}^{FIT}";
		Rinv[i-1+7].GetNormalizedPlot().Draw();
		canvas->cd(i+nCentralities);
		Rlcms[i-1+7].labels = ";m_{T} [GeV/c^{2}];R_{LCMS} / R_{LCMS}^{FIT}";
		Rlcms[i-1+7].GetNormalizedPlot().Draw();
	}
	canvas->SaveAs("output/all_EPOS_div.png");
	delete canvas;

	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
}
