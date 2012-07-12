#ifndef _FIT1DCOULD_HXX_
#define _FIT1DCOULD_HXX_

#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMath.h>

Double_t fungek(Double_t *x, Double_t *par);
bool fit1dcould(const char *fileName, Double_t &Rinv, Double_t &RinvE);

#endif
