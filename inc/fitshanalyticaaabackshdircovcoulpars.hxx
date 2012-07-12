#ifndef _FITSHANALYTICAAABACKSHDIRCOVCOULPARS_HXX_
#define _FITSHANALYTICAAABACKSHDIRCOVCOULPARS_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFitter.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>

#define _NOSAVE_
#define BINSC 30
#define BINSP 30

double getavk(TGraph *grk, Double_t minb, Double_t maxb);
double erfi(double ax);
Double_t rqval(Double_t Rpar, Double_t qval, int fsw);
Double_t funfull(Double_t *x, Double_t *par, Double_t *cv00, Double_t *cv20, Double_t *cv22);
Double_t myfunctionegg(Double_t *x, Double_t *par);
Double_t myfun20egg(Double_t *x, Double_t *par);
Double_t myfun22egg(Double_t *x, Double_t *par);
Double_t myfunctionggg(Double_t *x, Double_t *par);
Double_t myfun20ggg(Double_t *x, Double_t *par);
Double_t myfun22ggg(Double_t *x, Double_t *par);
void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag);
void GetPar(ifstream *inf, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax);
void fitshanalyticreal( char *pref,
   Double_t &Ro, Double_t &Rs, Double_t &Rl, Double_t &Ri, Double_t &lambda,
   Double_t &RoE, Double_t &RsE, Double_t &RlE, Double_t &RiE, Double_t &dlambda);

bool fitshanalyticaaabackshdircovcoulpars(const char *filname, 
  Double_t &Rout, Double_t &Rside, Double_t &Rlong, Double_t &Rinv, Double_t &lambda,
   Double_t &dRout, Double_t &dRside, Double_t &dRlong, Double_t &dRinv, Double_t &dlambda);


#endif