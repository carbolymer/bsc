#include "plotter.hxx"

using namespace std;

const unsigned int nParticles = 3;
const unsigned int markerSize = 1;
const Double_t mTRangeMin = 0;
const Double_t mTRangeMax = 1.5;
const Double_t RRangeMin = 2;
const Double_t RRangeMax = 8.5;

const Double_t particleMasses[3] = { // in GeV / c^2
	0.493677,// kaon
	0.1349766, // pion
	0.938272029 // proton
};

int main()
{
	Double_t Rout, Rside, Rlong, Rlcms, Rinv, lambda,
			dRout, dRside, dRlong, dRlcms, dRinv, dlambda, mT;

	// tuple filename, ktbinmin, ktbin max
	string labels[nParticles] = {"K-K b2", "\\pi-\\pi b2", "p-p b2"};
	string prefixes[nParticles] = {"kk", "pipi", "pp"};

	TCanvas *canvasRcomponents = new TCanvas("canvasRcomponents", "R_components", 1300, 900);
	canvasRcomponents->Divide(3,2);
	// TCanvas *canvasR = new TCanvas("canvasR", "R", 1400, 800);
	// canvasRcomponents->Divide(2,3);
	
	TLegend *legendLCMS = new TLegend(0.8,0.75,0.9,0.9);
	TLegend *legendInv = new TLegend(0.8,0.75,0.9,0.9);
	TLegend *legendOut = new TLegend(0.8,0.75,0.9,0.9);
	TLegend *legendSide = new TLegend(0.8,0.75,0.9,0.9);
	TLegend *legendLong = new TLegend(0.8,0.75,0.9,0.9);

	TGraphErrors *RoutGraph[nParticles];
	TGraphErrors *RsideGraph[nParticles];
	TGraphErrors *RlongGraph[nParticles];
	TGraphErrors *RlcmsGraph[nParticles];
	TGraphErrors *RinvGraph[nParticles];

	for(int i = 0; i < nParticles; ++i)
	{
		RlcmsGraph[i] = new TGraphErrors();
		RinvGraph[i] = new TGraphErrors();
		RoutGraph[i] = new TGraphErrors();
		RsideGraph[i] = new TGraphErrors();
		RlongGraph[i] = new TGraphErrors();
	}

	for(int j = 0; j < nParticles; ++j)
	{	
		std::cout << prefixes[j] << std::endl;
		fillGraph(prefixes[j] + std::string("_Rout.out"), RoutGraph[j], j);
		fillGraph(prefixes[j] + std::string("_Rside.out"), RsideGraph[j], j);
		fillGraph(prefixes[j] + std::string("_Rlong.out"), RlongGraph[j], j);
		fillGraph(prefixes[j] + std::string("_Rinv.out"), RinvGraph[j], j);
		fillGraph(prefixes[j] + std::string("_Rlcms.out"), RlcmsGraph[j], j);

		/* R_out */
		canvasRcomponents->cd(1);
		RoutGraph[j]->GetYaxis()->SetRangeUser(RRangeMin,RRangeMax);
		RoutGraph[j]->GetXaxis()->SetLimits(mTRangeMin,mTRangeMax);
		if(j == 0)
			RoutGraph[j]->Draw("A*");
		else
			RoutGraph[j]->Draw("SAMEP*");
		RoutGraph[j]->SetMarkerColor(1+j);
		RoutGraph[j]->SetLineColor(1+j);
		RoutGraph[j]->SetMarkerSize(markerSize);
		RoutGraph[j]->SetMarkerStyle(25);
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
		RsideGraph[j]->SetMarkerColor(1+j);
		RsideGraph[j]->SetLineColor(1+j);
		RsideGraph[j]->SetMarkerSize(markerSize);
		RsideGraph[j]->SetMarkerStyle(25);
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
		RlongGraph[j]->SetMarkerColor(1+j);
		RlongGraph[j]->SetLineColor(1+j);
		RlongGraph[j]->SetMarkerSize(markerSize);
		RlongGraph[j]->SetMarkerStyle(25);
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
		RlcmsGraph[j]->SetMarkerColor(1+j);
		RlcmsGraph[j]->SetLineColor(1+j);
		RlcmsGraph[j]->SetMarkerSize(markerSize);
		RlcmsGraph[j]->SetMarkerStyle(25);
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
		RinvGraph[j]->SetMarkerColor(1+j);
		RinvGraph[j]->SetLineColor(1+j);
		RinvGraph[j]->SetMarkerSize(markerSize);
		RinvGraph[j]->SetMarkerStyle(25);
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
	RinvGraph[0]->SetTitle("R_{inv};m_{T};R_{inv}(m_{T})");
	legendInv->SetFillColor(0);
	canvasRcomponents->cd(5);
	legendInv->Draw();

	canvasRcomponents->SaveAs("output.png");
	canvasRcomponents->SaveAs("output.root");
	return 0;
}

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle)
{
	char buffer[256];
	std::ifstream infile(fileName.c_str(), std::ifstream::in);
	double kTmin, kTmax, R, dR, mT;
	Int_t i;
	std::cout << "mT\tR\tdR" << std::endl;
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
		mT = TMath::Sqrt(
			TMath::Power((kTmin+kTmax)/2,2)
			+ TMath::Power(particleMasses[iParticle],2));
		std::cout << mT << "\t" << R << "\t+/- " << dR << std::endl;
		graph->SetPoint(i, mT, R);
		graph->SetPointError(i, 0, dR);
	}
}
