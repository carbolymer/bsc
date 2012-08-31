#include <iostream>
#include <string>
#include <TFile.h>

void extract(const char* cFfileName, const char* outFileName, const char *title);


int extractSH()
{
	extract("/data/disk1/Models/tpi_input/bb3m6/piku/outfilepikulcf31a.root",
		"piku31.png",
		""); // "\\pi-K (+,-)   \\beta_{T} #in #[]{0.35, 0.5}"

	extract("/data/disk1/Models/tpi_input/bb3m6/piku/outfilepikulcf32a.root",
		"piku32.png",
		""); // "\\pi-K (+,-)   \\beta_{T} #in #[]{0.5, 0.65}"

	extract("/data/disk1/Models/tpi_input/bb3m6/piku/outfilepikulcf33a.root",
		"piku33.png",
		""); // "\\pi-K (+,-)   \\beta_{T} #in #[]{0.65, 0.8}"

	extract("/data/disk1/Models/tpi_input/bb3m6/piku/outfilepikulcf34a.root",
		"piku34.png",
		""); // "\\pi-K (+,-)   \\beta_{T} #in #[]{0.8, 0.95}"
	return 1;
}

void extract(const char* cFfileName, const char* outFileName, const char *title)
{
	int i;
	const int nSH = 3;
	TH1D *c[nSH];
	double
		titleSize = 0.06,
		titleOffset = 1.3,
		labelOffset = 0.005,
		labelSize = 0.04;

	TFile *inFile = new TFile(cFfileName);
	TCanvas *canv = new TCanvas("crfncn","crfncn", 1800,600);


	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.06, "xyz");
	gStyle->SetPadTopMargin(0.07);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadRightMargin(0.001);

	c[0] = (TH1D*) inFile->Get("CfnReYlm00NonIdCYlmTrue");
	c[1] = (TH1D*) inFile->Get("CfnReYlm10NonIdCYlmTrue");
	c[2] = (TH1D*) inFile->Get("CfnReYlm11NonIdCYlmTrue");

	canv->Divide(3,1);
	for(i = 0;i < nSH; ++i)
	{
		canv->cd(i+1);
		c[i]->Rebin(3);
		c[i]->Draw();
		c[i]->SetAxisRange(0,1.6);
		c[i]->SetMarkerStyle(20);
		c[i]->GetYaxis()->SetTitleSize(titleSize);
		c[i]->GetYaxis()->SetTitleOffset(titleOffset);
		c[i]->GetYaxis()->SetLabelOffset(labelOffset);
		c[i]->GetYaxis()->SetLabelSize(labelSize);


		c[i]->GetXaxis()->SetTitleSize(titleSize);
		c[i]->GetXaxis()->SetTitleOffset(titleOffset-0.4);
		c[i]->GetXaxis()->SetLabelOffset(labelOffset);
		c[i]->GetXaxis()->SetLabelSize(labelSize);
	}


	c[0]->SetTitle((std::string(title)+std::string(";k*;C^{0}_{0}")).c_str());
	c[1]->SetTitle(";k*;#RgothicC^{0}_{1}");
	c[2]->SetTitle(";k*;#RgothicC^{1}_{1}");

	canv->SaveAs(outFileName);
	delete canv;
	// delete inFile;
}