#include <iostream>
#include <TFile.h>

TH1D *bDistribution;

double intgrl(double x)
{
	int nbin = bDistribution->GetBin(x);
	return  bDistribution->Integral(0,nbin);
}

int EPOScreateCumulativebDistribution()
{
	TFile *file = new TFile("epos_b_dist.root");
	TCanvas *canv = (TCanvas*) file->Get("bdist");
	bDistribution =(TH1D*) canv->FindObject("bDIST");

	bDistribution->Draw();

	return 0;
}