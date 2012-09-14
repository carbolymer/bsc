#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TStyle.h>

#define LABELTEST false

const int pairTypesNumber = 8;

typedef struct
{
	TH1D
		*cf1D,
		*cfSH00,
		*cfSH10,
		*cfSH11,
		*cfSH20,
		*cfSH22;
	std::string
		centrality,
		model,
		particle;
} CorrelationFunction;

const Int_t
		markerStyle = 0,
		markerColor = 0,
		lineColor = 4;
const Float_t
		titleSize = 0.055,
		titleOffset = 0.85,
		labelOffset = 0.005,
		labelSize = 0.04,
		legendX1 = 0.7,
		legendY1 = 0.8,
		legendX2 = 0.999,
		legendY2 = 0.977,
		markerSize = 1.3;

void loadCentralityBin(
		const char *directory,
		const char *model,
		const char *centrality,
		std::vector<CorrelationFunction> &pipi,
		std::vector<CorrelationFunction> &pp,
		std::vector<CorrelationFunction> &kk,
		std::vector<CorrelationFunction> &kp,
		std::vector<CorrelationFunction> &pik,
		std::vector<CorrelationFunction> &piku,
		std::vector<CorrelationFunction> &pip,
		std::vector<CorrelationFunction> &pipu
	);

void Draw(const char* canvasName, std::vector<CorrelationFunction> &plots);
void DrawSinglePlot(TH1D* correlationFunction, const char *title, unsigned int&);
// void Draw(const char* model, DrawConfiguration configuration, vector<CorrelationFunction> &plots);

int createCorrelationPlots()
{
	std::vector<CorrelationFunction> 
		pipi,
		pp,
		kk,
		kp,
		pik,
		piku,
		pip,
		pipu;

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/lhc0005",
		"hHKM",
		"0-5%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/lhc1020",
		"hHKM",
		"10-20%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/lhc2030",
		"hHKM",
		"10-20%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/lhc3040",
		"hHKM",
		"30-30%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/epos_0005",
		"EPOS",
		"0-5%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/epos_1020",
		"EPOS",
		"10-20%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/epos_2030",
		"EPOS",
		"20-30%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);


	loadCentralityBin(
		"/home/mgalazyn/workspace/tpi_output/epos_3040",
		"EPOS",
		"30-40%",
		pipi, pp, kk, kp, pik, piku, pip, pipu);

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(labelSize, "xyz");
	gStyle->SetPadTopMargin(0.09);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.11);
	gStyle->SetPadRightMargin(0.05);

	Draw("pipi",pipi);

	return 0;
}

void loadCentralityBin(
		const char *directory,
		const char *model,
		const char *centrality,
		std::vector<CorrelationFunction> &pipi,
		std::vector<CorrelationFunction> &pp,
		std::vector<CorrelationFunction> &kk,
		std::vector<CorrelationFunction> &kp,
		std::vector<CorrelationFunction> &pik,
		std::vector<CorrelationFunction> &piku,
		std::vector<CorrelationFunction> &pip,
		std::vector<CorrelationFunction> &pipu
	)
{
	CorrelationFunction cf;

	std::vector<CorrelationFunction> *correlationFunctions[pairTypesNumber];

	TChain 
		*chPipi,
		*chPp,
		*chKk,
		*chKp,
		*chPik,
		*chPiku,
		*chPip,
		*chPipu,
		*chains[pairTypesNumber];
	std::string sDirectory = directory;
	cf.model = model;
	cf.centrality = centrality;

	for(int i = 0; i < pairTypesNumber; ++i)
		chains[i] = new TChain;

	chPipi = chains[0];
	chPp = chains[1];
	chKk = chains[2];
	chKp = chains[3];
	chPik = chains[4];
	chPiku = chains[5];
	chPip = chains[6];
	chPipu = chains[7];

	correlationFunctions[0] = &pipi;
	correlationFunctions[1] = &pp;
	correlationFunctions[2] = &kk;
	correlationFunctions[3] = &kp;
	correlationFunctions[4] = &pik;
	correlationFunctions[5] = &piku;
	correlationFunctions[6] = &pip;
	correlationFunctions[7] = &pipu;

	if(LABELTEST)
	{
		for(int i = 0; i < pairTypesNumber; ++i)
			for(int j = 0; j < 1; ++j)
			{
				cf.cf1D = new TH1D;
				cf.cfSH00 = new TH1D;
				cf.cfSH10 = new TH1D;
				cf.cfSH11 = new TH1D;
				cf.cfSH20 = new TH1D;
				cf.cfSH22 = new TH1D;
				correlationFunctions[i]->push_back(cf);
			}
		return;
	}

	chPipi->Add( (sDirectory+std::string("/outfilecf*.root")).c_str() );
	chPp->Add( (sDirectory+std::string("/outfileppcf*.root")).c_str() );
	chKk->Add( (sDirectory+std::string("/outfilekkcf*.root")).c_str() );
	chKp->Add( (sDirectory+std::string("/outfilekpcf*.root")).c_str() );
	chPik->Add( (sDirectory+std::string("/outfilepikcf*.root")).c_str() );
	chPiku->Add( (sDirectory+std::string("/outfilepikulcf*.root")).c_str() );
	chPip->Add( (sDirectory+std::string("/outfilepipcf*.root")).c_str() );
	chPipu->Add( (sDirectory+std::string("/outfilepipulcf*.root")).c_str() );

	for(int i =0; i < pairTypesNumber ; ++i)
	{
		TObjArray *fileElements = chains[i]->GetListOfFiles();
        TIter next(fileElements);
        TChainElement *chEl=0;
		TFileMerger merger(kFALSE,kFALSE);
		merger.SetMsgPrefix("[merger]");
		merger.SetNotrees(kFALSE);
        while( (chEl = (TChainElement*)next()) )
		{
			if(!merger.AddFile(chEl->GetTitle()))
				std::cerr << "Could not add file: " << chEl->GetTitle() << std::endl;
			merger.Merge();
			// break; //////// DEBUG!!!
        }
		TFile *file = new TFile(merger.GetOutputFileName());

		if(file == 0)
		{
			std::cout << "NULL POINTER: *" << merger.GetOutputFileName() << "* !" << std::endl;
			continue;
		}
		if(file->IsZombie())
		{
			std::cout << "ZOMBIE: *" << merger.GetOutputFileName() << "* !" << std::endl;
			continue;
		}

		cf.cf1D = (TH1D*) file->Get("cnumn1da");
		cf.cf1D->Divide((TH1D*) file->Get("cdenn1da"));

		file->GetObject("CfnReYlm00NonIdCYlmTrue", cf.cfSH00);
		if(cf.cfSH00)
		{
			file->GetObject("CfnReYlm00NonIdCYlmTrue", cf.cfSH00);
			file->GetObject("CfnReYlm10NonIdCYlmTrue", cf.cfSH10);
			file->GetObject("CfnReYlm11NonIdCYlmTrue", cf.cfSH11);
			file->GetObject("CfnReYlm20NonIdCYlmTrue", cf.cfSH20);
			file->GetObject("CfnReYlm22NonIdCYlmTrue", cf.cfSH22);
		}

		file->GetObject("CfnReYlm00NonIdCYlm", cf.cfSH00);
		if(cf.cfSH00)
		{
			file->GetObject("CfnReYlm00NonIdCYlm", cf.cfSH00);
			file->GetObject("CfnReYlm10NonIdCYlm", cf.cfSH10);
			file->GetObject("CfnReYlm11NonIdCYlm", cf.cfSH11);
			file->GetObject("CfnReYlm20NonIdCYlm", cf.cfSH20);
			file->GetObject("CfnReYlm22NonIdCYlm", cf.cfSH22);
		}

		file->GetObject("CfnReYlm00IdLCYlm", cf.cfSH00);
		if(cf.cfSH00)
		{
			file->GetObject("CfnReYlm00IdLCYlm", cf.cfSH00);
			file->GetObject("CfnReYlm10IdLCYlm", cf.cfSH10);
			file->GetObject("CfnReYlm11IdLCYlm", cf.cfSH11);
			file->GetObject("CfnReYlm20IdLCYlm", cf.cfSH20);
			file->GetObject("CfnReYlm22IdLCYlm", cf.cfSH22);
		}

		correlationFunctions[i]->push_back(cf);
	}
}

void Draw(const char* canvasName, std::vector<CorrelationFunction> &plots)
{
	unsigned int i, k;
	TCanvas *canv = new TCanvas(canvasName, canvasName, 1800, 600);
	canv->Divide(6,2);

	// hHKM
	k = 0;
	canv->cd(1);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cf1D, ";q_{inv} (GeV/c);C(q_{inv})", k);
		}
			
	k = 0;
	canv->cd(2);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cfSH00, ";k* (GeV/c);#RgothicC^{0}_{0}", k);
		}

	k = 0;
	canv->cd(3);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cfSH10, ";k* (GeV/c);#RgothicC^{0}_{1}", k);
		}

	k = 0;
	canv->cd(4);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cfSH11, ";k* (GeV/c);#RgothicC^{1}_{1}", k);
		}

	k = 0;
	canv->cd(5);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cfSH20, ";k* (GeV/c);#RgothicC^{0}_{2}", k);
		}

	k = 0;
	canv->cd(6);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("hHKM") == 0)
		{
			DrawSinglePlot(plots[i].cfSH22, ";k* (GeV/c);#RgothicC^{2}_{2}", k);
		}

	// EPOS
	k = 0;
	canv->cd(7);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cf1D, ";q_{inv} (GeV/c);C(q_{inv})", k);
		}
			
	k = 0;
	canv->cd(8);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cfSH00, ";k* (GeV/c);#RgothicC^{0}_{0}", k);
		}

	k = 0;
	canv->cd(9);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cfSH10, ";k* (GeV/c);#RgothicC^{0}_{1}", k);
		}

	k = 0;
	canv->cd(10);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cfSH11, ";k* (GeV/c);#RgothicC^{1}_{1}", k);
		}

	k = 0;
	canv->cd(11);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cfSH20, ";k* (GeV/c);#RgothicC^{0}_{2}", k);
		}

	k = 0;
	canv->cd(12);
	for(i = 0 ; i < plots.size(); ++i)
		if(plots[i].model.compare("EPOS") == 0)
		{
			DrawSinglePlot(plots[i].cfSH22, ";k* (GeV/c);#RgothicC^{2}_{2}", k);
		}

	canv->SaveAs((std::string(canvasName)+std::string(".png")).c_str());
	canv->SaveAs((std::string(canvasName)+std::string(".root")).c_str());
}

void DrawSinglePlot(TH1D* correlationFunction, const char* title, unsigned int &i)
{
	if(i == 0)
	{
		correlationFunction->Draw();
		correlationFunction->GetYaxis()->SetTitleSize(titleSize);
		correlationFunction->GetYaxis()->SetTitleOffset(titleOffset);
		correlationFunction->GetYaxis()->SetLabelOffset(labelOffset);
		correlationFunction->GetXaxis()->SetLimits(0,0.2);
		correlationFunction->GetXaxis()->SetTitleSize(titleSize);
		correlationFunction->GetXaxis()->SetTitleOffset(titleOffset);
		correlationFunction->GetXaxis()->SetLabelOffset(labelOffset);
		correlationFunction->GetXaxis()->SetLabelSize(labelSize);
		correlationFunction->SetTitle(title);
		i = 1;
	}
	else
		correlationFunction->Draw("SAME");
	correlationFunction->SetMarkerColor(i++);
	correlationFunction->SetLineColor(1);
	correlationFunction->SetMarkerStyle(34);
}

// void Draw(const char* model, DrawConfiguration configuration, vector<CorrelationFunction> &plots)
// {
// 	configuration.axes->Draw();
// 	configuration.axes->GetYaxis()->SetRangeUser(configuration.yMin,configuration.yMax);
// 	configuration.axes->GetYaxis()->SetTitleSize(configuration.titleSize);
// 	configuration.axes->GetYaxis()->SetTitleOffset(configuration.titleOffset);
// 	configuration.axes->GetYaxis()->SetLabelOffset(configuration.labelOffset);
// 	configuration.axes->GetYaxis()->SetLabelSize(configuration.labelSize);
// 	configuration.axes->GetXaxis()->SetLimits(configuration.xMin,xMax);
// 	configuration.axes->GetXaxis()->SetTitleSize(configuration.titleSize);
// 	configuration.axes->GetXaxis()->SetTitleOffset(configuration.titleOffset+0.05);
// 	configuration.axes->GetXaxis()->SetLabelOffset(configuration.labelOffset);
// 	configuration.axes->GetXaxis()->SetLabelSize(configuration.labelSize);
// 	configuration.axes->SetTitle(labels.c_str());
	
// 	TLegend *_legend = new TLegend(configuration.legendX1, configuration.legendY1, configuration.legendX2, configuration.legendY2);
// 		_legend->Clear();
// 		for(int i = 0; i < graphCount; ++i)
// 		{
// 			graphs[i]->Draw("SAMEP*");
// 			graphs[i]->SetMarkerColor(theme[i].markerColor);
// 			graphs[i]->SetLineColor(theme[i].lineColor);
// 			graphs[i]->SetMarkerSize(theme[i].markerSize);
// 			graphs[i]->SetMarkerStyle(theme[i].markerStyle);
// 			_legend->AddEntry(graphs[i], graphNames[i].c_str(),"P");
// 		}

// 		_legend->Draw();
// 		_legend->SetFillColor(0);

// }
