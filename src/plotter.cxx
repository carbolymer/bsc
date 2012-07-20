#include "plotter.hxx"

using namespace std;

const unsigned int nParticles = 3;
const unsigned int nCentralities = 3;
const unsigned int nPlots = nParticles*nCentralities;
const Double_t markerSize = 1.6;
const Double_t mTRangeMin = 0;
const Double_t mTRangeMax = 1.5;
const Double_t RRangeMin = 1.6;
const Double_t RRangeMax = 11;

const Double_t particleMasses[3] = { // in GeV / c^2
	0.493677,// kaon
	0.1349766, // pion
	0.938272029 // proton
};

const int colors[9] = {1,2,4,6,8,9,13,28,46};
const int markers[9] = {20,21,22,23,29,33,34};

int main()
{
	Double_t Rout, Rside, Rlong, Rlcms, Rinv, lambda,
			dRout, dRside, dRlong, dRlcms, dRinv, dlambda, mT;

	int markerColor, markerStyle;

	// string labels[nPlots] = 
	// {
	// 	"K-K b5", "\\pi-\\pi b5", "p-p b5"
	// };

	// string prefixes[nPlots] = 
	// {
	// 	"b5/kk", "b5/pipi", "b5/pp"
	// };

	string labels[nPlots] = 
	{
		"K-K b2", "\\pi-\\pi b2", "p-p b2",
		"K-K b3", "\\pi-\\pi b3", "p-p b3",
		"K-K b5", "\\pi-\\pi b5", "p-p b5"
	};
	string prefixes[nPlots] = 
	{
		"b2/kk", "b2/pipi", "b2/pp",
		"b3/kk", "b3/pipi", "b3/pp",
		"b5/kk", "b5/pipi", "b5/pp"
	};

	TCanvas *canvasRcomponents = new TCanvas("canvasRcomponents", "R_components", 1300, 900);
	canvasRcomponents->Divide(3,2);
	
	TLegend *legendLCMS = new TLegend(0.65,0.55,0.9,0.9);
	TLegend *legendInv = new TLegend(0.65,0.55,0.9,0.9);
	TLegend *legendOut = new TLegend(0.65,0.55,0.9,0.9);
	TLegend *legendSide = new TLegend(0.65,0.55,0.9,0.9);
	TLegend *legendLong = new TLegend(0.65,0.55,0.9,0.9);

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
		markerStyle = markers[(int)(j/nCentralities)];
		markerColor = colors[j];

		std::cout << std::endl << "  ###  " << prefixes[j] << "  ###   " << std::endl;
		fillGraph(prefixes[j] + std::string("_Rout.out"), RoutGraph[j], j%3);
		fillGraph(prefixes[j] + std::string("_Rside.out"), RsideGraph[j], j%3);
		fillGraph(prefixes[j] + std::string("_Rlong.out"), RlongGraph[j], j%3);
		fillGraph(prefixes[j] + std::string("_Rinv.out"), RinvGraph[j], j%3, kTRUE);
		fillGraph(prefixes[j] + std::string("_Rlcms.out"), RlcmsGraph[j], j%3);

		/* R_out */
		canvasRcomponents->cd(1);
		RoutGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RoutGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RoutGraph[j]->Draw("A*");
		else
			RoutGraph[j]->Draw("SAMEP*");
		RoutGraph[j]->SetMarkerColor(markerColor);
		RoutGraph[j]->SetLineColor(markerColor);
		RoutGraph[j]->SetMarkerSize(markerSize);
		RoutGraph[j]->SetMarkerStyle(markerStyle);
		RoutGraph[j]->SetLineWidth(1);
		legendOut->AddEntry(RoutGraph[j], labels[j].c_str(),"P");

		/* R_side */
		canvasRcomponents->cd(2);
		RsideGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RsideGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RsideGraph[j]->Draw("A*");
		else
			RsideGraph[j]->Draw("SAMEP*");
		RsideGraph[j]->SetMarkerColor(markerColor);
		RsideGraph[j]->SetLineColor(markerColor);
		RsideGraph[j]->SetMarkerSize(markerSize);
		RsideGraph[j]->SetMarkerStyle(markerStyle);
		RsideGraph[j]->SetLineWidth(1);
		legendSide->AddEntry(RsideGraph[j], labels[j].c_str(),"P");

		/* R_long */
		canvasRcomponents->cd(3);
		RlongGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RlongGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RlongGraph[j]->Draw("A*");
		else
			RlongGraph[j]->Draw("SAMEP*");
		RlongGraph[j]->SetMarkerColor(markerColor);
		RlongGraph[j]->SetLineColor(markerColor);
		RlongGraph[j]->SetMarkerSize(markerSize);
		RlongGraph[j]->SetMarkerStyle(markerStyle);
		RlongGraph[j]->SetLineWidth(1);
		legendLong->AddEntry(RlongGraph[j], labels[j].c_str(),"P");

		/* R_LCMS */
		canvasRcomponents->cd(4);
		RlcmsGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RlcmsGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RlcmsGraph[j]->Draw("A*");
		else
			RlcmsGraph[j]->Draw("SAMEP*");
		RlcmsGraph[j]->SetMarkerColor(markerColor);
		RlcmsGraph[j]->SetLineColor(markerColor);
		RlcmsGraph[j]->SetMarkerSize(markerSize);
		RlcmsGraph[j]->SetMarkerStyle(markerStyle);
		RlcmsGraph[j]->SetLineWidth(1);
		legendLCMS->AddEntry(RlcmsGraph[j], labels[j].c_str(),"P");

		/* R_invariant */
		canvasRcomponents->cd(5);
		RinvGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RinvGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RinvGraph[j]->Draw("A*");
		else
			RinvGraph[j]->Draw("SAMEP*");
		RinvGraph[j]->SetMarkerColor(markerColor);
		RinvGraph[j]->SetLineColor(markerColor);
		RinvGraph[j]->SetMarkerSize(markerSize);
		RinvGraph[j]->SetMarkerStyle(markerStyle);
		RinvGraph[j]->SetLineWidth(1);
		legendInv->AddEntry(RinvGraph[j], labels[j].c_str(),"P");
	}
	/* R_out */
	RoutGraph[0]->SetTitle("R_{out};m_{T};R_{out}(m_{T})");
	legendOut->SetFillColor(0);
	canvasRcomponents->cd(1);
	legendOut->Draw();

	/* R_side */
	RsideGraph[0]->SetTitle("R_{side};m_{T};R_{side}(m_{T})");
	legendSide->SetFillColor(0);
	canvasRcomponents->cd(2);
	legendSide->Draw();

	/* R_long */
	RlongGraph[0]->SetTitle("R_{long};m_{T};R_{long}(m_{T})");
	legendLong->SetFillColor(0);
	canvasRcomponents->cd(3);
	legendLong->Draw();

	/* R_LCMS */
	RlcmsGraph[0]->SetTitle("R_{LCMS};m_{T};R_{LCMS}(m_{T})");
	legendLCMS->SetFillColor(0);
	canvasRcomponents->cd(4);
	legendLCMS->Draw();

	/* R_invariant */
	RinvGraph[0]->SetTitle("R_{inv}/[(\\sqrt{\\gamma}+2)/3]^{1/2};m_{T} (GeV/c^2);R_{inv}(m_{T}) (fm)");
	legendInv->SetFillColor(0);
	canvasRcomponents->cd(5);
	legendInv->Draw();

	canvasRcomponents->SaveAs("output.png");
	canvasRcomponents->SaveAs("output.root");
	return 0;
}

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle, bool isInvariant)
{
	char buffer[256];
	std::ifstream infile(fileName.c_str(), std::ifstream::in);
	double kT, kTmin, kTmax, R, dR, mT, gamma;
	Int_t i;
	std::cout << std::endl << fileName << std::endl << "mT\tR\tdR" << std::endl;
	while(infile.good())
	{
		std::fill(buffer, buffer+255, 0);
    	infile >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	i = graph->GetN();
    	kTmin = boost::lexical_cast<double>(buffer);
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
		std::cout << mT << "\t" << R << "\t+/- " << dR << std::endl;
		graph->SetPoint(i, mT, R);
		graph->SetPointError(i, 0, dR);
	}
}