#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <TLegend.h>
#include <TRandom.h>
#include "fitshanalyticaaabackshdircovcoulpars.hxx"
#include "merger.hxx"
#include "fit1dcould.hxx"

using namespace std;

int main()
{
	const unsigned int nParticles = 3;
	const unsigned int markerSize = 1;

	const Double_t particleMasses[3] = { // in GeV / c^2
		0.493677,// kaon
		0.1349766, // pion
		0.938272029 // proton
	};

	Double_t Rout, Rside, Rlong, Rlcms, Rinv, lambda,
			dRout, dRside, dRlong, dRlcms, dRinv, dlambda, mT;

	// tuple filename, ktbinmin, ktbin max
	vector<ktBinFile> particleFiles[nParticles];

	string labels[nParticles] = {"K - K", "\\pi - \\pi", "p - p"};

	ifstream kaonFilesList("filelist.kk.in", ifstream::in);
	ifstream pionFilesList("filelist.pipi.in", ifstream::in);
	ifstream protonFilesList("filelist.pp.in", ifstream::in);

	TCanvas *canvasRcomponents = new TCanvas("canvasRcomponents", "R_components", 1400, 800);
	canvasRcomponents->Divide(3,2);
	// TCanvas *canvasR = new TCanvas("canvasR", "R", 1400, 800);
	// canvasRcomponents->Divide(2,3);
	
	TLegend *legendLCMS = new TLegend(0.80,0.15,0.9,0.4);
	TLegend *legendInv = new TLegend(0.80,0.15,0.9,0.4);
	TLegend *legendOut = new TLegend(0.80,0.15,0.9,0.4);
	TLegend *legendSide = new TLegend(0.80,0.15,0.9,0.4);
	TLegend *legendLong = new TLegend(0.80,0.15,0.9,0.4);

	TGraphErrors *RlcmsGraph[nParticles];
	TGraphErrors *RinvGraph[nParticles];
	TGraphErrors *RoutGraph[nParticles];
	TGraphErrors *RsideGraph[nParticles];
	TGraphErrors *RlongGraph[nParticles];

	for(int i = 0; i < nParticles; ++i)
	{
		RlcmsGraph[i] = new TGraphErrors();
		RinvGraph[i] = new TGraphErrors();
		RoutGraph[i] = new TGraphErrors();
		RsideGraph[i] = new TGraphErrors();
		RlongGraph[i] = new TGraphErrors();
	}

	cout << endl << "Loading kaons..." << endl;
	loadFileList(kaonFilesList, particleFiles[0]);
	cout << endl << "Loading pions..." << endl;
	loadFileList(pionFilesList, particleFiles[1]);
	cout << endl << "Loading protons..." << endl;
	loadFileList(protonFilesList, particleFiles[2]);

	kaonFilesList.close();
	pionFilesList.close();
	protonFilesList.close();

	// TRandom *gRandom = new TRandom();

	for(int j = 0; j < nParticles; ++j)
	{	
		for(int i = 0, size = particleFiles[j].size(); i < size; ++i)
		{
			mT = TMath::Sqrt(
				TMath::Power((boost::get<1>(particleFiles[j][i])+boost::get<2>(particleFiles[j][i]))/2,2)
				+ TMath::Power(particleMasses[i],2));

			if(!fit1dcould(boost::get<0>(particleFiles[j][i]).c_str(), Rinv, dRinv))
			{
				cout << "1D Fit failed!" << endl;
				return 1;
			}

			// if(!fitshanalyticaaabackshdircovcoulpars(boost::get<0>(particleFiles[j][i]).c_str(),
			// 	Rout, Rside, Rlong, Rlcms, lambda,
	  //  			dRout, dRside, dRlong, dRlcms, dlambda))
			// {
			// 	cout << "3D Fit failed!" << endl;
			// 	return 1;
			// }

			/*R_long*/
			// Rlong = gRandom->Uniform(2.5/(i+1.0)+j,3.0/(i+1.0)+j);
			// dRlong = gRandom->Uniform(0.1,1.0/(i+1));

			RlongGraph[j]->SetPoint(i, mT, Rlong);
			RlongGraph[j]->SetPointError(i, 0, dRlong);

			/*R_side*/
			// Rside = gRandom->Uniform(2.5/(i+1.0)+j,3.0/(i+1.0)+j);
			// dRside = gRandom->Uniform(0.1,1.0/(i+1));

			RsideGraph[j]->SetPoint(i, mT, Rside);
			RsideGraph[j]->SetPointError(i, 0, dRside);	

			/*R_out*/
			// Rout = gRandom->Uniform(2.5/(i+1.0)+j,3.0/(i+1.0)+j);
			// dRout = gRandom->Uniform(0.1,1.0/(i+1));

			RoutGraph[j]->SetPoint(i, mT, Rout);
			RoutGraph[j]->SetPointError(i, 0, dRout);

			/* R_LCMS */
			// Rlcms = gRandom->Uniform(2.5/(i+1.0)+j,3.0/(i+1.0)+j);
			// dRlcms = gRandom->Uniform(0.1,1.0/(i+1));

			RlcmsGraph[j]->SetPoint(i, mT, Rlcms);
			RlcmsGraph[j]->SetPointError(i, 0, dRlcms);

			/* R_invariant */
			cout << "RINV" << Rinv << "+/-" << dRinv << endl;
			RinvGraph[j]->SetPoint(i, mT, Rinv);
			RinvGraph[j]->SetPointError(i, 0, dRinv);

		}
		/* R_out */
		canvasRcomponents->cd(1);
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

	canvasRcomponents->SaveAs("Routsidelong.png");

	return 0;
}