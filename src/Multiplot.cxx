#include <TH1F.h>
#include <TLegend.h>
#include <TGraphErrors.h>
const int graphCount = 3;

typedef struct
{
	Int_t markerColor;
	Int_t lineColor;
	Float_t markerSize;
	Int_t markerStyle;
} TGraphErrorsPreset;


class Multiplot
{
private:
	TH1F *_axes;
	TLegend *_legend;
public:
	Float_t
		titleSize,
		titleOffset,
		labelOffset,
		labelSize,
		xMin,
		xMax,
		yMin,
		yMax,
		legendX1,
		legendY1,
		legendX2,
		legendY2;
	TGraphErrors *kaons;
	TGraphErrors *pions;
	TGraphErrors *protons;
	TGraphErrors *graphs[graphCount];
	std::string labels;		// axes labels
	std::string graphNames[graphCount];
	TGraphErrorsPreset theme[graphCount];

	Multiplot(const char _labels[] = "") :
		titleSize(0.08),
		titleOffset(0.6),
		labelOffset(0.005),
		labelSize(0.07),
		xMin(0.1),
		xMax(1.5),
		yMin(2),
		yMax(10),
		legendX1(0.7),
		legendY1(0.75),
		legendX2(0.999),
		legendY2(0.977)
	{
		labels = _labels;
		_axes = new TH1F;
		_legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
		for(int i = 0; i < graphCount; ++i)
			graphs[i] = new TGraphErrors;
		kaons = graphs[0];
		pions = graphs[1];
		protons = graphs[2];
	}

	Multiplot(const Multiplot &rhs) :
		titleSize(rhs.titleSize),
		titleOffset(rhs.titleOffset),
		labelOffset(rhs.labelOffset),
		labelSize(rhs.labelSize),
		xMin(rhs.xMin),
		xMax(rhs.xMax),
		yMin(rhs.yMin),
		yMax(rhs.yMax),
		legendX1(rhs.legendX1),
		legendY1(rhs.legendY1),
		legendX2(rhs.legendX2),
		legendY2(rhs.legendY2),
		labels(rhs.labels)
	{
		_axes = new TH1F(*rhs._axes);
		_legend = new TLegend(*rhs._legend);
		for(int i = 0; i < graphCount; ++i)
		{
			graphNames[i] = rhs.graphNames[i];
			graphs[i] = new TGraphErrors(*rhs.graphs[i]);
			theme[i] = rhs.theme[i];
		}
		kaons = graphs[0];
		pions = graphs[1];
		protons = graphs[2];
	}

	Multiplot& operator=(const Multiplot &rhs)
	{
		titleSize = rhs.titleSize;
		titleOffset = rhs.titleOffset;
		labelOffset = rhs.labelOffset;
		labelSize = rhs.labelSize;
		xMin = rhs.xMin;
		xMax = rhs.xMax;
		yMin = rhs.yMin;
		yMax = rhs.yMax;
		legendX1 = rhs.legendX1;
		legendY1 = rhs.legendY1;
		legendX2 = rhs.legendX2;
		legendY2 = rhs.legendY2;
		labels = rhs.labels;

		_axes = new TH1F(*rhs._axes);
		_legend = new TLegend(*rhs._legend);
		for(int i = 0; i < graphCount; ++i)
		{
			graphNames[i] = rhs.graphNames[i];
			graphs[i] = new TGraphErrors(*rhs.graphs[i]);
			theme[i] = rhs.theme[i];
		}
		kaons = graphs[0];
		pions = graphs[1];
		protons = graphs[2];
		return *this;
	}

	void Draw()
	{
		_axes->Draw();
		_axes->GetYaxis()->SetRangeUser(yMin,yMax);
		_axes->GetYaxis()->SetTitleSize(titleSize);
		_axes->GetYaxis()->SetTitleOffset(titleOffset);
		_axes->GetYaxis()->SetLabelOffset(labelOffset);
		_axes->GetYaxis()->SetLabelSize(labelSize);
		_axes->GetXaxis()->SetLimits(xMin,xMax);
		_axes->GetXaxis()->SetTitleSize(titleSize);
		_axes->GetXaxis()->SetTitleOffset(titleOffset+0.3);
		_axes->GetXaxis()->SetLabelOffset(labelOffset);
		_axes->GetXaxis()->SetLabelSize(labelSize);
		_axes->SetTitle(labels.c_str());

		_legend->Clear();
		for(int i = 0; i < graphCount; ++i)
		{
			graphs[i]->Draw("SAMEP*");
			graphs[i]->SetMarkerColor(theme[i].markerColor);
			graphs[i]->SetLineColor(theme[i].lineColor);
			graphs[i]->SetMarkerSize(theme[i].markerSize);
			graphs[i]->SetMarkerStyle(theme[i].markerStyle);
			_legend->AddEntry(graphs[i], graphNames[i].c_str(),"P");
		}

		_legend->SetFillColor(0);
		_legend->Draw();
	}
};