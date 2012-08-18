#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TRandom.h>
#include "../src/Multiplot.cxx"

#define DEBUG false

const Double_t particleMasses[3] = { // in GeV / c^2
	0.493677,// kaon
	0.1349766, // pion
	0.938272029 // proton
};

const int numberOfCentralities = 6;

void createLcmsPlot(std::string *graphNames, std::string *prefixes, const char *fileName, Multiplot &Rinv, Multiplot &Rlcms);
void createPrfPlot(Multiplot Rinv[], Multiplot Rlcms[], int nCentralities);
void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant = kFALSE);

int plotter()
{	
	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.06, "xyz");
	gStyle->SetPadTopMargin(0.023);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.11);
	gStyle->SetPadRightMargin(0.001);

	std::string prefixes[graphCount];
	std::string graphNames[graphCount];

	Multiplot Rinv[numberOfCentralities];
	Multiplot Rlcms[numberOfCentralities];

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

	prefixes[0] = "bb3m6/kk";
	prefixes[1] = "bb3m6/pipi";
	prefixes[2] = "bb3m6/pp";
	graphNames[0] = "K-K bb3m6";
	graphNames[1] = "\\pi-\\pi bb3m6";
	graphNames[2] = "p-p bb3m6";
	createLcmsPlot(graphNames, prefixes, "bb3m6", Rinv[3], Rlcms[3]);

	prefixes[0] = "lhc0005/kk";
	prefixes[1] = "lhc0005/pipi";
	prefixes[2] = "lhc0005/pp";
	graphNames[0] = "K-K lhc0005";
	graphNames[1] = "\\pi-\\pi lhc0005";
	graphNames[2] = "p-p lhc0005";
	createLcmsPlot(graphNames, prefixes, "lhc0005", Rinv[4], Rlcms[4]);

	prefixes[0] = "lhc1020/kk";
	prefixes[1] = "lhc1020/pipi";
	prefixes[2] = "lhc1020/pp";
	graphNames[0] = "K-K lhc1020";
	graphNames[1] = "\\pi-\\pi lhc1020";
	graphNames[2] = "p-p lhc1020";
	createLcmsPlot(graphNames, prefixes, "lhc1020", Rinv[5], Rlcms[5]);

	createPrfPlot(Rinv, Rlcms, numberOfCentralities);
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

void createLcmsPlot(std::string *graphNames, std::string *prefixes, const char *fileName, Multiplot &Rinv, Multiplot &Rlcms)
{
	TCanvas *canvas = new TCanvas("canvas", "R_LCMS", 900, 800);
	canvas->Divide(2,2);

	Multiplot Rout(";m_{T} [GeV/c^{2}];R_{out} [fm]");

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

	Multiplot Rside(Rout);
	Multiplot Rlong(Rout);
	Rlcms = Multiplot(Rout);
	Rinv = Multiplot(Rout);

	Rside.labels = ";m_{T} [GeV/c^{2}];R_{side} [fm]";
	Rlong.labels = ";m_{T} [GeV/c^{2}];R_{long} [fm]";
	Rlcms.labels = ";m_{T} [GeV/c^{2}];R_{LCMS} [fm]";
	Rinv.labels = ";m_{T} [GeV/c^{2}];R_{inv}/[(\\sqrt{\\gamma}+2)/3]^{1/2} [fm]";

	canvas->cd(1);
	Rout.Draw();

	canvas->cd(2);
	Rside.Draw();

	canvas->cd(3);
	Rlong.Draw();

	canvas->cd(4);
	Rlcms.Draw();

	for( int j = 0; j < graphCount; ++j)
	{
		if(DEBUG) std::cout << std::endl << "  ###  " << prefixes[j] << "  ###   " << std::endl;
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rout.out"), Rout.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rside.out"), Rside.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rlong.out"), Rlong.graphs[j], j);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rinv.out"), Rinv.graphs[j], j, kTRUE);
		fillGraph(std::string("data/") + prefixes[j] + std::string("_Rlcms.out"), Rlcms.graphs[j], j);
	}
	canvas->SaveAs((std::string("output/")+fileName+std::string(".png")).c_str());
	delete canvas;
}

void createPrfPlot(Multiplot *Rinv, Multiplot *Rlcms, int nCentralities)
{
	TCanvas *canvas = new TCanvas("canvas", "R_LCMS", 2500, 600);
	canvas->Divide(nCentralities, 2);

	for(int i = 1; i <= nCentralities; ++i)
	{
		canvas->cd(i);
		Rinv[i-1].Draw();
		canvas->cd(i+nCentralities);
		Rlcms[i-1].Draw();
	}
	canvas->SaveAs("output/all.png");
	delete canvas;
}
