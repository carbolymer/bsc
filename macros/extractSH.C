#include <iostream>
#include <string>
#include <TFile.h>

void extract(const char* cFfileName, const char* outFileName, const char *title);

int extractSH()
{
	// "\\pi-K (+,-)   \\beta_{T} #in #[]{0.35, 0.5}"

	// hkm 0-5
	extract("/data/disk1/Models/tpi_input/lhc0005/piku/outfilepikulcf35a.root",
		"hhkm0_piku35.png",
		"\\pi-K unlike sign hHKM 0-5%");
	extract("/data/disk1/Models/tpi_input/lhc0005/pik/outfilepikcf35a.root",
		"hhkm0_pik35.png",
		"\\pi-K like sign hHKM 0-5%");
	extract("/data/disk1/Models/tpi_input/lhc0005/pipu/outfilepipulcf35a.root",
		"hhkm0_pipu35.png",
		"\\pi-p unlike sign hHKM 0-5%");
	extract("/data/disk1/Models/tpi_input/lhc0005/pip/outfilepipcf35a.root",
		"hhkm0_pip35.png",
		"\\pi-p like sign hHKM 0-5%");
	extract("/data/disk1/Models/tpi_input/lhc0005/kp/outfilekpcf35a.root",
		"hhkm0_kp35.png",
		"p-K like sign hHKM 0-5%");


	// hkm 10-20
	extract("/data/disk1/Models/tpi_input/lhc1020/piku/outfilepikulcf35a.root",
		"hhkm1_piku35.png",
		"\\pi-K unlike sign hHKM 10-20%");
	extract("/data/disk1/Models/tpi_input/lhc1020/pik/outfilepikcf35a.root",
		"hhkm1_pik35.png",
		"\\pi-K like sign hHKM 10-20%");
	extract("/data/disk1/Models/tpi_input/lhc1020/pipu/outfilepipulcf35a.root",
		"hhkm1_pipu35.png",
		"\\pi-p unlike sign hHKM 10-20%");
	extract("/data/disk1/Models/tpi_input/lhc1020/pip/outfilepipcf35a.root",
		"hhkm1_pip35.png",
		"\\pi-p like sign hHKM 10-20%");
	extract("/data/disk1/Models/tpi_input/lhc1020/kp/outfilekpcf35a.root",
		"hhkm1_kp35.png",
		"p-K like sign hHKM 10-20%");


	// hhkm 20-30
	extract("/data/disk1/Models/tpi_input/lhc2030/piku/outfilepikulcf35a.root",
		"hhkm2_piku35.png",
		"\\pi-K unlike sign hHKM 20-30%");
	extract("/data/disk1/Models/tpi_input/lhc2030/pik/outfilepikcf35a.root",
		"hhkm2_pik35.png",
		"\\pi-K like sign hHKM 20-30%");
	extract("/data/disk1/Models/tpi_input/lhc2030/pipu/outfilepipulcf35a.root",
		"hhkm2_pipu35.png",
		"\\pi-p unlike sign hHKM 20-30%");
	extract("/data/disk1/Models/tpi_input/lhc2030/pip/outfilepipcf35a.root",
		"hhkm2_pip35.png",
		"\\pi-p like sign hHKM 20-30%");
	extract("/data/disk1/Models/tpi_input/lhc2030/kp/outfilekpcf35a.root",
		"hhkm2_kp35.png",
		"p-K like sign hHKM 20-30%");

	// hhkm 30-40
	extract("/data/disk1/Models/tpi_input/lhc3040/piku/outfilepikulcf35a.root",
		"hhkm3_piku35.png",
		"\\pi-K unlike sign hHKM 30-40%");
	extract("/data/disk1/Models/tpi_input/lhc3040/pik/outfilepikcf35a.root",
		"hhkm3_pik35.png",
		"\\pi-K like sign hHKM 30-40%");
	extract("/data/disk1/Models/tpi_input/lhc3040/pipu/outfilepipulcf35a.root",
		"hhkm3_pipu35.png",
		"\\pi-p unlike sign hHKM 30-40%");
	extract("/data/disk1/Models/tpi_input/lhc3040/pip/outfilepipcf35a.root",
		"hhkm3_pip35.png",
		"\\pi-p like sign hHKM 30-40%");
	extract("/data/disk1/Models/tpi_input/lhc3040/kp/outfilekpcf35a.root",
		"hhkm3_kp35.png",
		"p-K like sign hHKM 30-40%");

	return 1;
}

void extract(const char* cFfileName, const char* outFileName, const char *title)
{
	int i;
	const int nSH = 5;
	TH1D *c[nSH];
	double
		titleSize = 0.06,
		titleOffset = 1.3,
		labelOffset = 0.005,
		labelSize = 0.04;

	TFile *inFile = new TFile(cFfileName);
	TCanvas *canv = new TCanvas("crfncn","crfncn", 1800,1200);


	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.06, "xyz");
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadRightMargin(0.001);

	c[0] = (TH1D*) inFile->Get("CfnReYlm00NonIdCYlmTrue");
	c[1] = (TH1D*) inFile->Get("CfnReYlm10NonIdCYlmTrue");
	c[2] = (TH1D*) inFile->Get("CfnReYlm11NonIdCYlmTrue");
	c[3] = (TH1D*) inFile->Get("CfnReYlm20NonIdCYlmTrue");
	c[4] = (TH1D*) inFile->Get("CfnReYlm22NonIdCYlmTrue");


	canv->Divide(3,2);
	for(i = 0;i < nSH; ++i)
	{
		canv->cd(i+1);
		// c[i]->Rebin(3);
		// c[i]->Scale(1./3.);
		c[i]->Draw();
		c[i]->SetAxisRange(0,0.4);
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

	c[0]->SetAxisRange(0.5,1.5,"Y");
	c[1]->SetAxisRange(-0.01,0.015,"Y");
	c[2]->SetAxisRange(-0.03,0.03,"Y");
	c[3]->SetAxisRange(-0.02,0.015,"Y");
	c[4]->SetAxisRange(-0.02,0.02,"Y");

	c[0]->SetTitle(";k*;C^{0}_{0}");
	c[1]->SetTitle((std::string(title)+std::string(";k*;#RgothicC^{0}_{1}")).c_str());
	c[2]->SetTitle(";k*;#RgothicC^{1}_{1}");
	c[3]->SetTitle(";k*;#RgothicC^{0}_{2}");
	c[4]->SetTitle(";k*;#RgothicC^{2}_{2}");

	canv->SaveAs(outFileName);
	delete canv;
	// delete inFile;
}