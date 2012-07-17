#ifndef _PLOTTER_HXX_
#define _PLOTTER_HXX_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TMath.h>

void fillGraph(std::string fileName, TGraphErrors *graph, unsigned int iParticle);

#endif