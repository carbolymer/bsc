#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>

using namespace std;

TH1D *bDistribution, *bCumulativeDistribution;

double intgrl(double x)
{
	int nbin = bDistribution->GetBin(x);
	return  bDistribution->Integral(0,nbin);
}

int EPOScreateCumulativebDistribution()
{
	double percent = 0, cumulativeValue = 0, integral;
	TFile *file = new TFile("epos_b_dist.root");
	TCanvas *canv = (TCanvas*) file->Get("bdist");
	bDistribution =(TH1D*) canv->FindObject("bDIST");
	canv = new TCanvas("bdist","bdist",800,800);

	bCumulativeDistribution = new TH1D("bCumulativeDist", "Impact parameter cumulative distribution", bDistribution->GetNbinsX(), 0,20);
	integral = bDistribution->Integral(0, bDistribution->GetNbinsX());

	for(int i = 0; i < bDistribution->GetNbinsX(); ++i)
	{
		cumulativeValue += bDistribution->GetBinContent(i)/integral;
		bCumulativeDistribution->SetBinContent(i,cumulativeValue);
	}
	bCumulativeDistribution->Draw();
	canv->SaveAs("epos_b_cumulative_dist.root");
	cout.precision(5);
	cout << "\tbmin\tbmax" << endl;
	for(percent = 0; percent < 1; percent+= 0.05)
	{
		cout << percent*100 << "-" << percent*100+5 << "%\t" << bCumulativeDistribution->GetBinCenter(bCumulativeDistribution->FindFirstBinAbove(percent))
			<< "\t" <<  bCumulativeDistribution->GetBinCenter(bCumulativeDistribution->FindFirstBinAbove(percent+0.05)) << endl;
	}

	return 0;
}