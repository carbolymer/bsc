#include "fit1dcould.hxx"

TGraph *calckcoulggg;

Double_t fungek(Double_t *x, Double_t *par)
{
  Double_t qinv2 = x[0]*x[0];
  
  Double_t norm = par[0];
  Double_t lam  = par[1];
  Double_t rad  = par[2]/0.197327;
  Double_t expf = par[3];
  Double_t exps = par[4];
  Double_t q2sl = par[5];
  Double_t qbgs = par[6];

  Int_t nkbin = (int) (x[0]/0.002);
  if (nkbin>499) nkbin=499;
  Double_t Kc = calckcoulggg->GetY()[nkbin]; 
  Double_t gcpart = TMath::Exp(-qinv2*rad*rad);
  Double_t ecpart = (exps*exps*exps*exps)/((exps*exps+qinv2)*(exps*exps+qinv2));

  return norm * ((1-lam-expf) + lam*Kc*(1+gcpart) + expf*Kc*(1+ecpart) + q2sl*qinv2 + qbgs*x[0]);
}

bool fit1dcould(const char *fileName, Double_t &Rinv, Double_t &RinvE) 
{
  TFile *inFile = new TFile(fileName);
  if(inFile->IsZombie())
    return kFALSE;

  TH1D *ratq = (TH1D*) inFile->Get("cnumn1da");
  TH1D *den = (TH1D*) inFile->Get("cdenn1da");
  ratq->Divide(den);

  TFile *ifk = new TFile("ffcomp.root");
  calckcoulggg = (TGraph *) ifk->Get("KCoulomb");

  TF1 *funq = new TF1("funq",fungek,0.0,1.0,7);
  
  funq->SetParName(0,"Normalization");
  funq->SetParName(1,"Lambda");
  funq->SetParName(2,"Radius [fm]");
  funq->SetParName(3,"Exp fraction");
  funq->SetParName(4,"Exp slope");
  funq->SetParName(5,"Slope val");
  funq->SetParName(6,"Shift val");

  funq->SetParameter(0,1.0);
  funq->SetParameter(1,0.4);
  funq->SetParameter(2,1.0);
  //funq->SetParameter(3,0.4);
  funq->FixParameter(3,0.0);
  //funq->SetParameter(4,0.01);
  funq->FixParameter(4,0.0);
  funq->FixParameter(5,0.0);
  //funq->FixParameter(6,0.0);

//   funq->FixParameter(5,-6.51437e+00);
//   funq->FixParameter(6,-1.01904e+00);

  TFitResultPtr result = ratq->Fit(funq,"NS","",0.001,0.5);

  Rinv = result->Value(2);
  RinvE = result->FitResult::Error(2);
  return kTRUE;
}

int main(int argc, char **argv)
{
  Double_t Rinv, dRinv;

  if(argc < 5)
  {
    std::cout << "Example:" << std::endl
      << argv[0] << " file ktmin ktmax type" << std::endl
      << argv[0] << " outfilekk51a.root 0.12 0.4 kk" << std::endl;
      return 1;
  }
  std::string pairType(argv[3]);

  fit1dcould(argv[1], Rinv, dRinv);

  std::ofstream RinvFile((std::string(argv[4]) + std::string("_Rinv.out")).c_str(), std::ifstream::app);

  RinvFile << argv[2] << "\t" << argv[3] << "\t" << Rinv << "\t" << dRinv << std::endl;
  RinvFile.close();
}
