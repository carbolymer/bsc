#include <iostream>

using namespace std;

const unsigned int 
	nParticles = 3,
	nCentralities = 3,
	nPlots = nCentralities*nParticles;

const Double_t
		titleSize = 0.06,
		labelSize = 0.04,
		markerSize = 1.6,
		mTRangeMin = 0,
		mTRangeMax = 1.5,
		RRangeMin = 1.6,
		RRangeMax = 11;

const Double_t particleMasses[3] = { // in GeV / c^2
	0.493677,// kaon
	0.1349766, // pion
	0.938272029 // proton
};

const int colors[9] = {1,2,4,6,8,9,13,28,46};
const int markers[9] = {20,21,22,23,29,33,34};

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant = kFALSE);

int plotter()
{
	Double_t Rout, Rside, Rlong, Rlcms, Rinv, lambda,
			dRout, dRside, dRlong, dRlcms, dRinv, dlambda, mT;

	int markerColor, markerStyle;

	const string labels[nPlots] = 
	{
		"K-K b2", "\\pi-\\pi b2", "p-p b2",
		"K-K b3", "\\pi-\\pi b3", "p-p b3",
		"K-K b5", "\\pi-\\pi b5", "p-p b5"
		// "K-K bb3m6", "\\pi-\\pi bb3m6", "p-p bb3m6",
	};
	const string prefixes[nPlots] = 
	{
		"b2/kk", "b2/pipi", "b2/pp",
		"b3/kk", "b3/pipi", "b3/pp",
		"b5/kk", "b5/pipi", "b5/pp"
		// "bb3m6/kk", "bb3m6/pipi", "bb3m6/pp"
	};

	TCanvas *canvasRcomponents = new TCanvas("canvasRcomponents", "R_components", 1500, 900);
	canvasRcomponents->Divide(nCentralities,2);
	
	TLegend *legendLCMS[nCentralities];
	legendLCMS[0] = new TLegend(0.75,0.65,0.9,0.9);
	for(int i = 1; i < nCentralities; ++i)
		legendLCMS[i] = new TLegend(*legendLCMS[0]);
	
	TLegend *legendInv[nCentralities];
	for(int i = 0; i < nCentralities; ++i)
		legendInv[i] = new TLegend(*legendLCMS[0]);
	TLegend *legendOut = new TLegend(*legendLCMS[0]);
	TLegend *legendSide = new TLegend(*legendLCMS[0]);
	TLegend *legendLong = new TLegend(*legendLCMS[0]);

	TGraphErrors *RoutGraph[nPlots];
	TGraphErrors *RsideGraph[nPlots];
	TGraphErrors *RlongGraph[nPlots];
	TGraphErrors *RlcmsGraph[nPlots];
	TGraphErrors *RinvGraph[nPlots];

	for(int i = 0; i < nPlots; ++i)
	{
		RlcmsGraph[i] = new TGraphErrors();
		RinvGraph[i] = new TGraphErrors();
		RoutGraph[i] = new TGraphErrors();
		RsideGraph[i] = new TGraphErrors();
		RlongGraph[i] = new TGraphErrors();
	}

	for(int j = 0; j < nPlots; ++j)
	{	
		markerStyle = markers[(int)(j/nCentralities)%3];
		markerColor = colors[j%3];

		// std::cout << std::endl << "  ###  " << prefixes[j] << "  ###   " << std::endl;
		fillGraph(string("data/") + prefixes[j] + std::string("_Rout.out"), RoutGraph[j], j%3);
		fillGraph(string("data/") + prefixes[j] + std::string("_Rside.out"), RsideGraph[j], j%3);
		fillGraph(string("data/") + prefixes[j] + std::string("_Rlong.out"), RlongGraph[j], j%3);
		fillGraph(string("data/") + prefixes[j] + std::string("_Rinv.out"), RinvGraph[j], j%3, kTRUE);
		fillGraph(string("data/") + prefixes[j] + std::string("_Rlcms.out"), RlcmsGraph[j], j%3);

		// /* R_out */
		// canvasRcomponents->cd(1);
		// RoutGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		// RoutGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		// if(j == 0)
		// 	RoutGraph[j]->Draw("A*");
		// else
		// 	RoutGraph[j]->Draw("SAMEP*");
		// RoutGraph[j]->SetMarkerColor(markerColor);
		// RoutGraph[j]->SetLineColor(markerColor);
		// RoutGraph[j]->SetMarkerSize(markerSize);
		// RoutGraph[j]->SetMarkerStyle(markerStyle);
		// RoutGraph[j]->SetLineWidth(1);
		// RoutGraph[j]->GetXaxis()->SetTitleSize(titleSize);
		// RoutGraph[j]->GetXaxis()->SetLabelSize(labelSize);
		// RoutGraph[j]->GetYaxis()->SetTitleSize(titleSize);
		// RoutGraph[j]->GetYaxis()->SetLabelSize(labelSize);
		// legendOut->AddEntry(RoutGraph[j], labels[j].c_str(),"P");

		// /* R_side */
		// canvasRcomponents->cd(2);
		// RsideGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		// RsideGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		// if(j == 0)
		// 	RsideGraph[j]->Draw("A*");
		// else
		// 	RsideGraph[j]->Draw("SAMEP*");
		// RsideGraph[j]->SetMarkerColor(markerColor);
		// RsideGraph[j]->SetLineColor(markerColor);
		// RsideGraph[j]->SetMarkerSize(markerSize);
		// RsideGraph[j]->SetMarkerStyle(markerStyle);
		// RsideGraph[j]->SetLineWidth(1);
		// RsideGraph[j]->GetXaxis()->SetTitleSize(titleSize);
		// RsideGraph[j]->GetXaxis()->SetLabelSize(labelSize);
		// RsideGraph[j]->GetYaxis()->SetTitleSize(titleSize);
		// RsideGraph[j]->GetYaxis()->SetLabelSize(labelSize);
		// legendSide->AddEntry(RsideGraph[j], labels[j].c_str(),"P");

		// /* R_long */
		// canvasRcomponents->cd(3);
		// RlongGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		// RlongGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		// if(j == 0)
		// 	RlongGraph[j]->Draw("A*");
		// else
		// 	RlongGraph[j]->Draw("SAMEP*");
		// RlongGraph[j]->SetMarkerColor(markerColor);
		// RlongGraph[j]->SetLineColor(markerColor);
		// RlongGraph[j]->SetMarkerSize(markerSize);
		// RlongGraph[j]->SetMarkerStyle(markerStyle);
		// RlongGraph[j]->SetLineWidth(1);
		// RlongGraph[j]->GetXaxis()->SetTitleSize(titleSize);
		// RlongGraph[j]->GetXaxis()->SetLabelSize(labelSize);
		// RlongGraph[j]->GetYaxis()->SetTitleSize(titleSize);
		// RlongGraph[j]->GetYaxis()->SetLabelSize(labelSize);
		// legendLong->AddEntry(RlongGraph[j], labels[j].c_str(),"P");

		/* R_LCMS */
		canvasRcomponents->cd(nCentralities+1+(int)(j/nCentralities));
		RlcmsGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RlcmsGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j % nCentralities == 0)
			RlcmsGraph[j]->Draw("A*");
		else
			RlcmsGraph[j]->Draw("SAMEP*");
		RlcmsGraph[j]->SetMarkerColor(markerColor);
		RlcmsGraph[j]->SetLineColor(markerColor);
		RlcmsGraph[j]->SetMarkerSize(markerSize);
		RlcmsGraph[j]->SetMarkerStyle(markerStyle);
		RlcmsGraph[j]->SetLineWidth(1);
		RlcmsGraph[j]->GetXaxis()->SetTitleSize(titleSize);
		RlcmsGraph[j]->GetXaxis()->SetLabelSize(labelSize);
		RlcmsGraph[j]->GetYaxis()->SetTitleSize(titleSize);
		RlcmsGraph[j]->GetYaxis()->SetLabelSize(labelSize);
		legendLCMS[(int)j/nCentralities]->AddEntry(RlcmsGraph[j], labels[j].c_str(),"P");

		/* R_invariant */
		canvasRcomponents->cd(1+(int)(j/nCentralities));
		RinvGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RinvGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j % nCentralities == 0)
			RinvGraph[j]->Draw("A*");
		else
			RinvGraph[j]->Draw("SAMEP*");
		RinvGraph[j]->SetMarkerColor(markerColor);
		RinvGraph[j]->SetLineColor(markerColor);
		RinvGraph[j]->SetMarkerSize(markerSize);
		RinvGraph[j]->SetMarkerStyle(markerStyle);
		RinvGraph[j]->SetLineWidth(1);
		RinvGraph[j]->GetXaxis()->SetTitleSize(titleSize);
		RinvGraph[j]->GetXaxis()->SetLabelSize(labelSize);
		RinvGraph[j]->GetYaxis()->SetTitleSize(titleSize);
		RinvGraph[j]->GetYaxis()->SetLabelSize(labelSize);
		legendInv[(int)j/nCentralities]->AddEntry(RinvGraph[j], labels[j].c_str(),"P");
	}
	// /* R_out */
	// RoutGraph[0]->SetTitle(";m_{T} (GeV/c^{2});R_{out} (fm)");
	// legendOut->SetFillColor(0);
	// canvasRcomponents->cd(1);
	// legendOut->Draw();

	// /* R_side */
	// RsideGraph[0]->SetTitle(";m_{T} (GeV/c^{2});R_{side} (fm)");
	// legendSide->SetFillColor(0);
	// canvasRcomponents->cd(2);
	// legendSide->Draw();

	// /* R_long */
	// RlongGraph[0]->SetTitle(";m_{T} (GeV/c^{2});R_{long} (fm)");
	// legendLong->SetFillColor(0);
	// canvasRcomponents->cd(3);
	// legendLong->Draw();

	/* R_LCMS */
	for(int i = 0; i < nPlots; ++i)
		if(i%3 == 0)
			RlcmsGraph[i]->SetTitle(";m_{T} (GeV/c^{2});R_{LCMS} (fm)");
		
	for(int i = 0; i < nCentralities; ++i)
	{
		legendLCMS[i]->SetFillColor(0);
		canvasRcomponents->cd(nCentralities+1+i);
		legendLCMS[i]->Draw();
	}


	/* R_invariant */
	for(int i = 0; i < nPlots; ++i)
		if(i%3 == 0)
			RinvGraph[i]->SetTitle(";m_{T} (GeV/c^{2});R_{inv}/[(\\sqrt{\\gamma}+2)/3]^{1/2} (fm)");

	for(int i = 0; i < nCentralities; ++i)
	{
		legendInv[i]->SetFillColor(0);
		canvasRcomponents->cd(i+1);
		legendLCMS[i]->Draw();
	}


	canvasRcomponents->SaveAs("output/all.png");
	return 0;
}

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant)
{
	char buffer[256];
	std::ifstream infile(fileName.c_str(), std::ifstream::in);
	double kT, kTmin, kTmax, R, dR, mT, gamma;
	Int_t i;
	// std::cout << std::endl << fileName << std::endl << "mT\tR\tdR" << std::endl;
	while(infile.good())
	{
		for(int i=0; i < 256; ++i)
			buffer[i] = '\0';
    	infile >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	i = graph->GetN();

    	stringstream(buffer) >> kTmin;
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
		// std::cout << mT << "\t" << R << "\t+/- " << dR << std::endl;
		graph->SetPoint(i, mT, R);
		graph->SetPointError(i, 0, dR);
	}
}