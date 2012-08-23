#ifndef _MULTIFITPLOT_CXX_
#define _MULTIFITPLOT_CXX_

#include <TF1.h>
#include "MultiPlot.cxx"


class MultiFitPlot: public MultiPlot
{
private:
	TF1 *_fittingFunctions[graphCount];

protected:
	MultiFitPlot& _copyIntoMe(const MultiFitPlot &rhs)
	{
		for(int i = 0; i < graphCount; ++i)
			_fittingFunctions[i] = new TF1(*rhs._fittingFunctions[i]);
		MultiPlot::_copyIntoMe(rhs);
		return *this;
	}

public:
	static TF1* defaultFunction;

	MultiFitPlot() : MultiPlot("")
	{
		for(int i = 0; i < graphCount; ++i)
			_fittingFunctions[i] = new TF1(*defaultFunction);
	}

	MultiFitPlot(const char _labels[]) : MultiPlot(_labels)
	{
		for(int i = 0; i < graphCount; ++i)
			_fittingFunctions[i] = new TF1(*defaultFunction);
	}

	MultiFitPlot(const TF1* fittingFunctions[], const char _labels[] = "") : MultiPlot(_labels)
	{
		for(int i = 0; i < graphCount; ++i)
			_fittingFunctions[i] = (TF1*) fittingFunctions[i];
	}

	MultiFitPlot(const MultiFitPlot &rhs) : MultiPlot("")
	{
		_copyIntoMe(rhs);
	}

	MultiPlot& operator=(const MultiFitPlot &rhs)
	{
		return _copyIntoMe(rhs);
	}

	void Fit()
	{
		Double_t *ey, *x, *y, it =0, nPoints = 0;
		for(int i = 0; i < graphCount; ++i)
			nPoints += graphs[i]->GetN();

		TGraphErrors *allPoints = new TGraphErrors();
		for(int i = 0; i < graphCount; ++i)
		{
			nPoints = graphs[i]->GetN();
			x = graphs[i]->GetX();
			ey = graphs[i]->GetEY();
			y = graphs[i]->GetY();
			for(int j = 0; j < nPoints; ++j)
			{
				allPoints->SetPoint(it,x[j],y[j]);
				allPoints->SetPointError(it,0,y[j]);
				++it;
			}
		}
		_fittingFunctions[0]->SetParameter(0,1);
		_fittingFunctions[0]->SetParameter(1,1);
		allPoints->Fit(_fittingFunctions[0],"S");
	}

	void Draw()
	{
		MultiPlot::Draw();
		// for(int i = 0; i < graphCount; ++i)
		// {
			_fittingFunctions[0]->Draw("SAME");
			_fittingFunctions[0]->SetLineColor(1);
			// if(graphs[i]->GetN() < 2)
			// 	_fittingFunctions[i]->SetLineWidth(0);
			_fittingFunctions[0]->SetLineWidth(2);
		// }
	}

	MultiPlot& GetNormalizedPlot()
	{
		MultiPlot *mp = new MultiPlot();
		int i, j, nPoints;
		Double_t *yErrors;
		Double_t x,y;
		mp->labels = labels;
		mp->yMin = 0.7;
		mp->yMax = 1.4;
		mp->legendY1 = 0.95;
		for(j=0 ; j < graphCount; ++j )
		{
			mp->theme[j] = theme[j];
			mp->graphNames[j] = graphNames[j];
		}
		for(i = 0; i < graphCount; ++i)
		{
			nPoints = graphs[i]->GetN();
			yErrors = graphs[i]->GetEY();
			for(j = 0; j < nPoints; ++j)
			{
				// if(nPoints >= 2)
				{
					graphs[i]->GetPoint(j, x, y);
					mp->graphs[i]->SetPoint(j, x, y/_fittingFunctions[0]->Eval(x));
					mp->graphs[i]->SetPointError(j, 0, yErrors[i]/_fittingFunctions[0]->Eval(x));
				}
			}
		}
		return *mp;
	}
};

TF1*  MultiFitPlot::defaultFunction = 0;

#endif