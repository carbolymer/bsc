#include "fitshanalyticaaabackshdircovcoulpars.hxx"

using namespace std;

double gLmL, gRoL, gRsL, gRlL;
double gLmS, gRoS, gRsS, gRlS;

Int_t dofix[20];

Double_t parsg[20];
Double_t parmin[20];
Double_t parmax[20];

double getavk(TGraph *grk, Double_t minb, Double_t maxb)
{
  Int_t bins = 0; 
  Int_t bine = 0;
  
  while ((minb>grk->GetX()[bins]) && (bins<grk->GetN())) bins++;
  if (bins == grk->GetN()) return 1.0;
  while ((maxb>grk->GetX()[bine]) && (bine<grk->GetN())) bine++;
  
  Double_t wsum=0.0;
  Double_t ksum=0.0;

  for (int iter=bins; iter<bine; iter++) {
    ksum += grk->GetY()[iter]*grk->GetX()[iter]*grk->GetX()[iter];
    wsum += grk->GetX()[iter]*grk->GetX()[iter];
  }
  
  return ksum/wsum;
}

double erfi(double ax)
{
  double oneoversqpi = 1.0/TMath::Sqrt(TMath::Pi());

  if (ax > 0)
    return ( 2*ax +
	     2*ax*ax*ax/3.0 +
	     ax*ax*ax*ax*ax/5.0 +
	     ax*ax*ax*ax*ax*ax*ax/21.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax/108.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/780.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/5400.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/42840.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/383040.0
	     )*oneoversqpi;
  else
    return ( 2*ax -
	     2*ax*ax*ax/3.0 +
	     ax*ax*ax*ax*ax/5.0 -
	     ax*ax*ax*ax*ax*ax*ax/21.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax/108.0 -
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/780.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/5400.0 -
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/42840.0 +
	     ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax*ax/383040.0
	     )*oneoversqpi;
}

TGraph *grc00;
TGraph *grc00f;
TGraph *grc20f;
TGraph *grc22f;

TGraph *grkcoul;

Int_t binb, bine;
Double_t fitb, fite;

Int_t isprf;
Int_t funout, funside, funlong;

Double_t rqval(Double_t Rpar, Double_t qval, int fsw)
{
  switch (fsw) {
  case 0:
    return exp(-Rpar * qval);
    break;
  case 1:
    return exp(-TMath::Sqrt(Rpar * qval));
    break;
  case 2:
    return 1.0/(qval*Rpar+1.0);
    break;
  }
}

Double_t funfull(Double_t *x, Double_t *par, Double_t *cv00, Double_t *cv20, Double_t *cv22)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Lm = par[0];
  Double_t Ro = par[1];
  Double_t Rs = par[2];
  Double_t Rl = par[3];
  Double_t Bt = par[4];
  Double_t Rbo = par[5];
  Double_t Rbs = par[6];
  Double_t Rbl = par[7];
  Double_t Cf20 = par[8];
  Double_t Cf22 = par[9];
  Double_t Dv20 = par[10];
  Double_t Dv22 = par[11];
  Double_t Qbeg = par[12];
  Double_t Qslp = par[13];
  
  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double Rbo2 = Rbo*Rbo;
  double Rbs2 = Rbs*Rbs;
  double Rbl2 = Rbl*Rbl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  double fv = 0.0;
  double fv00s=0.0, fv20s=0.0, fv22s=0.0;

  double um = 1.0;
  if (xx>Qbeg) um = 1 + Qslp*(xx-Qbeg)*(xx-Qbeg); 

  double kavval = getavk(grkcoul, xx-0.005, xx+0.005);

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);

      fv = ((1-Lm) + kavval*Lm*(1+rqval(Ro2,qo2,funout) * rqval(Rs2,qs2,funside) * rqval(Rl2,ql2,funlong)))*(1+Bt*rqval(Rbo2,qo2,0) * rqval(Rbs2,qs2,0) * rqval(Rbl2,ql2,0));
        
      //      cout << "Got q cth phi fv " << xx << " " << cth << " " << phi << "   " << fv << endl;

      fv00s += fv;
      fv20s += fv*(3*cth*cth - 1);
      fv22s += fv*sth2*TMath::Cos(2*phi);
    
      //      fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      //      fv += rqval(Ro2,qo2,funout) * rqval(Rs2,qs2,funside) * rqval(Rl2,ql2,funlong);
      // 	fv20 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(3*cth*cth-1);
      // 	fv22 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(1-cth*cth)*TMath::Cos(2*phi);
    }
  }
  fv00s /= BINSC*BINSP;
  fv20s /= BINSC*BINSP;
  fv22s /= BINSC*BINSP;
  
  (*cv00) = fv00s * um;
  (*cv20) = fv20s * Cf20 - Dv20;
  (*cv22) = fv22s * Cf22 - Dv22;
}

Double_t myfunctionegg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];
  Double_t qb = par[5];
  Double_t qsl = par[6];
  Double_t qbm;

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.0000001) {
//     if (xx<qb) qbm = 1.0;
//     else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;
    
//     return val * qbm;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  double fv = 0.0;
//     double fv20 = 0.0;
//     double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //      fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      fv += rqval(Ro2,qo2,funout) * rqval(Rs2,qs2,funside) * rqval(Rl2,ql2,funlong);
      // 	fv20 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(3*cth*cth-1);
      // 	fv22 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(1-cth*cth)*TMath::Cos(2*phi);
    }
  }
  fv /= BINSC*BINSP;
//     fv20 /= 100*200;
//     fv22 /= 100*200;

  if (xx<qb) qbm = 1.0;
  else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

  return (val + norm * fv)*qbm;

}

Double_t myfun20egg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.0000001) {    
//     return val;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  //  double fv = 0.0;
  double fv20 = 0.0;
//     double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      //      fv20 += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2)*(3*cth2-1);
      fv20 += rqval(Ro2,qo2,funout) * rqval(Rs2,qs2,funside) * rqval(Rl2,ql2,funlong)*(3*cth2-1);
      // 	fv22 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(1-cth*cth)*TMath::Cos(2*phi);
    }
  }
  //  fv /= 100*100;
  fv20 /= BINSC*BINSP;
//     fv22 /= 100*200;

//   if (xx<qb) qbm = 1.0;
//   else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

  return (val + norm * fv20);//*qbm;

}

Double_t myfun22egg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.0000001) {    
//     return val;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  //  double fv = 0.0;
//     double fv20 = 0.0;
  double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //      fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      // 	fv20 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(3*cth*cth-1);
      //      fv22 += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2)*sth2*TMath::Cos(2*phi);
      fv22 += rqval(Ro2,qo2,funout) * rqval(Rs2,qs2,funside) * rqval(Rl2,ql2,funlong)*sth2*TMath::Cos(2*phi);
    }
  }
  //  fv /= 100*100;
//     fv20 /= 100*200;
  fv22 /= BINSC*BINSP;

//   if (xx<qb) qbm = 1.0;
//   else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

  return (val + norm * fv22);//*qbm;

}

Double_t myfunctionggg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];
  Double_t qb = par[5];
  Double_t qsl = par[6];
  Double_t qbm;

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.00001) {
//     if (xx<qb) qbm = 1.0;
//     else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;
    
//     return val * qbm;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  double fv = 0.0;
//     double fv20 = 0.0;
//     double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //      fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      fv += rqval(Ro2,qo2,0) * rqval(Rs2,qs2,0) * rqval(Rl2,ql2,0);
      // 	fv20 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(3*cth*cth-1);
      // 	fv22 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(1-cth*cth)*TMath::Cos(2*phi);
    }
  }
  fv /= BINSC*BINSP;
//     fv20 /= 100*200;
//     fv22 /= 100*200;

//   if (xx<qb) qbm = 1.0;
//   else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

//  return (val + norm * fv)*qbm;
  return (val + norm * fv);

}

Double_t myfun20ggg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.00001) {    
//     return val;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  //  double fv = 0.0;
  double fv20 = 0.0;
//     double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      fv20 += rqval(Ro2,qo2,0) * rqval(Rs2,qs2,0) * rqval(Rl2,ql2,0)*(3*cth2-1);
      // 	fv22 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(1-cth*cth)*TMath::Cos(2*phi);
    }
  }
  //  fv /= 100*100;
  fv20 /= BINSC*BINSP;
//     fv22 /= 100*200;

//   if (xx<qb) qbm = 1.0;
//   else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

  return (val + norm * fv20);//*qbm;
}

Double_t myfun22ggg(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  if (isprf) xx *= 2;
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  Double_t val = par[3];
  Double_t norm = par[4];

//   double Rmin = Ro;
//   if (Rs < Rmin) Rmin = Rs;
//   if (Rl < Rmin) Rmin = Rl;
//   if (exp(-xx*xx*Rmin*Rmin) < 0.00001) {    
//     return val;
//   }

  double Ro2 = Ro*Ro;
  double Rs2 = Rs*Rs;
  double Rl2 = Rl*Rl;

  double qv = xx;
  double qv2 = xx*xx;
  double phi, cth, sth2, cth2;
  double qo2, qs2, ql2;
  //  double fv = 0.0;
//     double fv20 = 0.0;
  double fv22 = 0.0;

  for (int icth=0; icth<BINSC; icth++) {
    cth = (icth+0.5)*2.0/BINSC - 1.0;
    cth2 = cth*cth;
    sth2 = (1-cth2);
    ql2 = qv2*cth2;

    for (int iphi=0; iphi<BINSP; iphi++) {
      phi = (iphi+0.5)*TMath::Pi()*2/BINSP;

      qo2 = qv2*sth2*TMath::Cos(phi)*TMath::Cos(phi);
      qs2 = qv2*sth2*TMath::Sin(phi)*TMath::Sin(phi);
      
      //      fv += exp(-Ro2*qo2 - Rs2*qs2 - Rl2*ql2);
      // 	fv20 += exp(-Ro*Ro*qo - Rs*Rs*qs - Rl*Rl*ql)*(3*cth*cth-1);
      fv22 += rqval(Ro2,qo2,0) * rqval(Rs2,qs2,0) * rqval(Rl2,ql2,0)*sth2*TMath::Cos(2*phi);
    }
  }
  //  fv /= 100*100;
//     fv20 /= 100*200;
  fv22 /= BINSC*BINSP;

//   if (xx<qb) qbm = 1.0;
//   else qbm = 1.0 + (xx-qb)*(xx-qb)*qsl;

  return (val + norm * fv22);//*qbm;
}

TF1 *f1;
TF1 *f20;
TF1 *fpod20;
TF1 *f22;

TF1 *f10S;
TF1 *f20S;
TF1 *f22S;

TH1D *c00;
TH1D *c20;
TH1D *c22;
TH3D *ccov;

void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t Ro = par[0];
  Double_t Rs = par[1];
  Double_t Rl = par[2];
  
  Double_t n0 = par[3];
  Double_t n1 = par[4];
  Double_t n2 = par[5];
  Double_t ng = par[6];
  Double_t nl = par[7];
  Double_t nt = par[8];
  
  Double_t qb = par[9];
  Double_t qs = par[10];

  Double_t RoS = par[11];
  Double_t RsS = par[12];
  Double_t RlS = par[13];
  Double_t n0S = par[14];

  Double_t c20lam = par[15];
  Double_t c20rad = par[16];
  Double_t c20hval = par[17];
  Double_t c20mval = par[18];
  Double_t c20wval = par[19];

  Double_t chi2 = 0.0;
  Double_t tv00, tv20, tv22, qv;
  Double_t cv00, cv20, cv22;
  Double_t ce00, ce20, ce22;
  Double_t cov0020, cov0022, cov2022;
  Double_t dl00, dl20, dl22;

//   f1 ->SetParameters(Ro, Rs, Rl, ng , n0, qb, qs);
//   f20->SetParameters(Ro, Rs, Rl, nl, n0*n1);
//   f22->SetParameters(Ro, Rs, Rl, nt, n0*n2);
  
//   f10S->SetParameters(RoS, RsS, RlS, ng , n0S, qb, qs);
//   f20S->SetParameters(RoS, RsS, RlS, nl, n0S*n1);
//   f22S->SetParameters(RoS, RsS, RlS, nt, n0S*n2);
  
//   Double_t binmin = 0; 
//   Double_t binmax = 0;
//   Double_t kavval = 0;
//   Double_t c00val = 0;

  Double_t xx[1];
  Double_t part[14];
  part[0] = n0;
  part[1] = Ro;
  part[2] = Rs;
  part[3] = Rl;
  part[4] = n0S;
  part[5] = RoS;
  part[6] = RsS;
  part[7] = RlS;
  part[8] = n1;
  part[9] = n2;
  part[10] = nl;
  part[11] = nt;
  part[12] = qb;
  part[13] = qs;

  for (int iter=binb; iter<bine; iter++) {
    //    cout << "iter " << c00->GetNbinsX() << endl;
    //     binmin = c00->GetXaxis()->GetBinLowEdge(iter);
    //     binmax = c00->GetXaxis()->GetBinUpEdge(iter);
    //     kavval = getavk(grkcoul, binmin/2.0, binmax/2.0);

    qv = c00->GetXaxis()->GetBinCenter(iter);
    xx[0] = qv;

    funfull(xx, part, &tv00, &tv20, &tv22);
    cv00 = c00->GetBinContent(iter);
    ce00 = c00->GetBinError(iter);
    dl00 = tv00*ng -cv00;

    cv20 = c20->GetBinContent(iter);
    ce20 = c20->GetBinError(iter);
    dl20 = (tv20*ng+c20lam*exp(-c20rad*c20rad*qv*qv)+c20hval*exp(-(qv-c20mval)*(qv-c20mval)*c20wval*c20wval)) - cv20;

    cv22 = c22->GetBinContent(iter);
    ce22 = c22->GetBinError(iter);
    dl22 = tv22*ng - cv22;

    cov0020 = ccov->GetBinContent(iter, 1, 13);
    cov0022 = ccov->GetBinContent(iter, 1, 17);
    cov2022 = ccov->GetBinContent(iter, 13, 17);
    
    //    printf("Covs are %lf %lf %lf\n", cov0020, cov0022, cov2022);

    chi2 += (dl00*dl00/(ce00*ce00) + 
	     dl22*dl22/(ce22*ce22) +
	     dl20*dl20/(ce20*ce20) -
	     (dl00/ce00)*(dl20/ce20)*(cov0020/(ce00*ce20)) -
	     (dl00/ce00)*(dl22/ce22)*(cov0022/(ce00*ce22)) -
	     (dl20/ce20)*(dl22/ce22)*(cov2022/(ce20*ce22)));

  }

  f = chi2;
}

void GetPar(ifstream *inf, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax)
{
  char linebuf[1000];
  char buf[100];
  memset(buf, 0, 100);

  inf->getline(linebuf, 1000);
  istringstream *istr = new istringstream(linebuf);
  (*istr) >> (*parval);
  (*istr) >> buf;

  cout << "Read value |" << *parval << "| |" << buf << "|";
  if (strstr(buf, "F")) { 
    *isfixed = 1; 
    cout << " Fixed" << endl; 
  }
  else if (strstr(buf, "L")) {
    *isfixed = 2;
    cout << " Limits ";
    (*istr) >> (*parmin) >> (*parmax);
    cout << "|" << *parmin << "| |" << (*parmax) << "|" << endl;
  }
  else { 
    *isfixed = 0; 
    cout << endl; 
  }
  
  delete istr;
}

void fitshanalyticreal( char *pref,
   Double_t &Rout, Double_t &Rside, Double_t &Rlong, Double_t &Rinv, Double_t &lambda,
   Double_t &RoutE, Double_t &RsideE, Double_t &RlongE, Double_t &RinvE, Double_t &dlambda)
{
  TVirtualFitter *fitter;


  double xy[100];
  double ykc00r[100];

  double ykc00rf[100];
  double ykc20rf[100];
  double ykc22rf[100];

  double Ro = gRoL/0.197327;
  double Rs = gRsL/0.197327;
  double Rl = gRlL/0.197327;

  f1 = new TF1("myf00",myfunctionegg,0,0.6,7);
  f1->SetParameters(Ro, Rs, Rl,0.0,1.0,0.5,0.0);

  f20 = new TF1("myf20",myfun20egg,0,0.6,5);
  f20->SetParameters(Ro, Rs, Rl,0.0,1.0);

  f22 = new TF1("myf22",myfun22egg,0,0.6,5);
  f22->SetParameters(Ro, Rs, Rl,0.0,1.0);

  f10S = new TF1("myf00S",myfunctionggg,0,0.6,7);
  f10S->SetParameters(Ro, Rs, Rl,0.0,1.0,0.5,0.0);

  f20S = new TF1("myf20S",myfun20ggg,0,0.6,5);
  f20S->SetParameters(Ro, Rs, Rl,0.0,1.0);


  f22S = new TF1("myf22S",myfun22ggg,0,0.6,5);
  f22S->SetParameters(Ro, Rs, Rl,0.0,1.0);

  char hname[200];
  sprintf(hname, "CfnReYlm00%s", pref);
  c00 = (TH1D *) gDirectory->Get(hname);
  if (!c00) { printf("Could not get histogram %s\n", hname); exit(0); }
  sprintf(hname, "CfnReYlm20%s", pref);
  c20 = (TH1D *) gDirectory->Get(hname);
  sprintf(hname, "CfnReYlm22%s", pref);
  c22 = (TH1D *) gDirectory->Get(hname);
  sprintf(hname, "CovCfc%s", pref);
  ccov = (TH3D *) gDirectory->Get(hname);

  if (TMath::IsNaN(c00->GetBinContent(1))) {
    c00->SetBinContent(1,0.0);
    c00->SetBinError(1,0.0);

    c20->SetBinContent(1,0.0);
    c20->SetBinError(1,0.0);
    c22->SetBinContent(1,0.0);
    c22->SetBinError(1,0.0);
  }

  binb = c00->GetXaxis()->FindFixBin(fitb);
  bine = c00->GetXaxis()->FindFixBin(fite);

  fitter=TVirtualFitter::Fitter(0, 20);
  fitter->SetFCN(myfuncf);
//   fitter->SetParameter(0, "Rout"     ,Ro        ,0.001, 0.5,   100.0);
//   fitter->SetParameter(1, "Rside"    ,Rs        ,0.001, 0.5,   100.0);
//   fitter->SetParameter(2, "Rlong"    ,Rl        ,0.001, 0.5,   100.0);
//   fitter->SetParameter(3, "lambda"   ,gLmL       ,0.0001,0.0,   40.01);
//   fitter->SetParameter(4, "20coef"   ,1.14       ,0.0001,-10.0, 10.0);
//   fitter->SetParameter(5, "22coef"   ,1.25   ,0.0001,-10.0, 10.0);
//   fitter->SetParameter(6, "norm"     ,1.0       ,0.0001,0.0,   2.0);
//   fitter->SetParameter(7, "norm20"   ,0.0       ,0.0001,-1.0,   1.0);
//   fitter->SetParameter(8, "norm22"   ,0.0       ,0.0001,-1.0,   1.0);
//   fitter->SetParameter(9, "qbeg"     ,1.0       ,0.001,  0.5,   1.5);
//   fitter->SetParameter(10,"qslope"   ,0.01      ,0.0001, 0.0,  10.0);
//   fitter->SetParameter(11,"RoutS"    ,gRoS/0.197327    ,0.0001, 0.01,  10.0);
//   fitter->SetParameter(12,"RSideS"   ,gRsS/0.197327    ,0.0001, 0.01,  10.0);
//   fitter->SetParameter(13,"RlongS"   ,gRlS/0.197327    ,0.0001, 0.01,  10.0);
//   fitter->SetParameter(14,"lambdaS"  ,gLmS       ,0.0001, 0.0,  40.01); 

  fitter->SetParameter(0, "Rout"     ,parsg[0]/0.197327,  0.001, parmin[0]/0.197327,  parmax[0]/0.197327);
  fitter->SetParameter(1, "Rside"    ,parsg[1]/0.197327,  0.001, parmin[1]/0.197327,  parmax[1]/0.197327);
  fitter->SetParameter(2, "Rlong"    ,parsg[2]/0.197327,  0.001, parmin[2]/0.197327,  parmax[2]/0.197327);
  fitter->SetParameter(3, "lambda"   ,parsg[3],  0.001, parmin[3],  parmax[3]);
  fitter->SetParameter(4, "20coef"   ,parsg[4],  0.001, parmin[4],  parmax[4]);
  fitter->SetParameter(5, "22coef"   ,parsg[5],  0.001, parmin[5],  parmax[5]);
  fitter->SetParameter(6, "norm"     ,parsg[6],  0.001, parmin[6],  parmax[6]);
  fitter->SetParameter(7, "norm20"   ,parsg[7],  0.001, parmin[7],  parmax[7]);
  fitter->SetParameter(8, "norm22"   ,parsg[8],  0.001, parmin[8],  parmax[8]);
  fitter->SetParameter(9, "qbeg"     ,parsg[9],  0.001, parmin[9],  parmax[9]);
  fitter->SetParameter(10,"qslope"   ,parsg[10], 0.001, parmin[10], parmax[10]);
  fitter->SetParameter(11,"RoutS"    ,parsg[11]/0.197327, 0.001, parmin[11]/0.197327, parmax[11]/0.197327);
  fitter->SetParameter(12,"RSideS"   ,parsg[12]/0.197327, 0.001, parmin[12]/0.197327, parmax[12]/0.197327);
  fitter->SetParameter(13,"RlongS"   ,parsg[13]/0.197327, 0.001, parmin[13]/0.197327, parmax[13]/0.197327);
  fitter->SetParameter(14,"lambdaS"  ,parsg[14], 0.001, parmin[14], parmax[14]); 
  fitter->SetParameter(15,"c20lam"   ,parsg[15], 0.001, parmin[15], parmax[15]); 
  fitter->SetParameter(16,"c20rad"   ,parsg[16], 0.001, parmin[16], parmax[16]); 
  fitter->SetParameter(17,"c20gval"  ,parsg[17], 0.001, parmin[17], parmax[17]); 
  fitter->SetParameter(18,"c20mval"  ,parsg[18], 0.001, parmin[18], parmax[18]); 
  fitter->SetParameter(19,"c20wval"  ,parsg[19], 0.001, parmin[19], parmax[19]); 

  for (int ipar=0; ipar<20; ipar++) {
    if (dofix[ipar] == 1)
      fitter->FixParameter(ipar);
  }

  Double_t arglist[100];
  arglist[0] = 1;
  fitter->ExecuteCommand("CALL FCN", arglist, 1);
  //    fitter->FixParameter(0);
  //    fitter->FixParameter(1);
  //    fitter->FixParameter(2);
  //   fitter->FixParameter(3);
//    fitter->FixParameter(4);
//    fitter->FixParameter(5);
   //   fitter->FixParameter(6);
   //     fitter->FixParameter(11);
//   fitter->FixParameter(14);

//   fitter->FixParameter(3);
//   fitter->FixParameter(4);
//   fitter->FixParameter(15);
  
//   fitter->FixParameter(7);
//   fitter->FixParameter(8);
//   fitter->FixParameter(9);
//    fitter->FixParameter(11);
//    fitter->FixParameter(12);
//    fitter->FixParameter(13);
//    fitter->FixParameter(14);
  

  //fitter->FixParameter(15);
  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT", arglist, 1);
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
  
//   f1 ->SetParameters(fitter->GetParameter(0), fitter->GetParameter(1), fitter->GetParameter(2), 1.0, n0);
//   f20->SetParameters(Ro, Rs, Rl, 1.0, n0*n1);
//   f22->SetParameters(Ro, Rs, Rl, 1.0, n0*n2);


  double pars[20];
  pars[0] = fitter->GetParameter(0);
  pars[1] = fitter->GetParameter(1);
  pars[2] = fitter->GetParameter(2);
  pars[3] = fitter->GetParameter(3);
  pars[4] = fitter->GetParameter(4);
  pars[5] = fitter->GetParameter(5);
  pars[6] = fitter->GetParameter(6);
  pars[7] = fitter->GetParameter(7);
  pars[8] = fitter->GetParameter(8);
  pars[9] = fitter->GetParameter(9);
  pars[10] = fitter->GetParameter(10);
  pars[11] = fitter->GetParameter(11);
  pars[12] = fitter->GetParameter(12);
  pars[13] = fitter->GetParameter(13);
  pars[14] = fitter->GetParameter(14);
  pars[15] = fitter->GetParameter(15);
  pars[16] = fitter->GetParameter(16);
  pars[17] = fitter->GetParameter(17);
  pars[18] = fitter->GetParameter(18);
  pars[19] = fitter->GetParameter(19);

  double  chimult;
  int iter;
  double xyv[2];
  myfuncf(iter, xyv, chimult, pars, 1);
  cout << "chi2 is " << chimult;
  chimult /= (3*(bine-binb)-9);
  chimult = TMath::Sqrt(chimult);
  cout << "   " << chimult << endl;

  cout << "Norm   " << fitter->GetParameter(6) << " +/- " << fitter->GetParError(6)*chimult << endl;
  cout << "Lambda " << fitter->GetParameter(3) << " +/- " << fitter->GetParError(3)*chimult << endl;
  lambda = fitter->GetParError(3);
  dlambda = lambda*chimult;
  cout << "Rout   " << fitter->GetParameter(0)*0.197327 << " +/- " << fitter->GetParError(0)*0.197327*chimult << endl;
  cout << "Rside  " << fitter->GetParameter(1)*0.197327 << " +/- " << fitter->GetParError(1)*0.197327*chimult << endl;
  cout << "Rlong  " << fitter->GetParameter(2)*0.197327 << " +/- " << fitter->GetParError(2)*0.197327*chimult << endl;

  cout << "Dev20  " << fitter->GetParameter(7) << " +/- " << fitter->GetParError(7)*chimult << endl;
  cout << "Dev22  " << fitter->GetParameter(8) << " +/- " << fitter->GetParError(8)*chimult << endl;

  cout << "QBeg   " << fitter->GetParameter(9) << " +/- " << fitter->GetParError(9)*chimult << endl;
  cout << "QSlope " << fitter->GetParameter(10) << " +/- " << fitter->GetParError(10)*chimult << endl;

  cout << "LambdaS " << fitter->GetParameter(14) << " +/- " << fitter->GetParError(14)*chimult << endl;
  cout << "RoutS  " << fitter->GetParameter(11)*0.197327 << " +/- " << fitter->GetParError(11)*0.197327*chimult << endl;
  cout << "RsideS " << fitter->GetParameter(12)*0.197327 << " +/- " << fitter->GetParError(12)*0.197327*chimult << endl;
  cout << "RlongS " << fitter->GetParameter(13)*0.197327 << " +/- " << fitter->GetParError(13)*0.197327*chimult << endl;

  cout << "c20lam  " << fitter->GetParameter(15) << " +/- " << fitter->GetParError(15)*chimult << endl;
  cout << "c20rad  " << fitter->GetParameter(16) << " +/- " << fitter->GetParError(16)*chimult << endl;
  cout << "c20hval " << fitter->GetParameter(17) << " +/- " << fitter->GetParError(17)*chimult << endl;
  cout << "c20mval " << fitter->GetParameter(18) << " +/- " << fitter->GetParError(18)*chimult << endl;
  cout << "c20wval " << fitter->GetParameter(19) << " +/- " << fitter->GetParError(19)*chimult << endl;

  Ro = Rout = fitter->GetParameter(0)*0.197327;
  Double_t RoE = RoutE = fitter->GetParError(0)*0.197327*chimult;
  Rs = Rside = fitter->GetParameter(1)*0.197327;
  Double_t RsE = RsideE = fitter->GetParError(1)*0.197327*chimult;
  Rl = Rlong = fitter->GetParameter(2)*0.197327;
  Double_t RlE = RlongE = fitter->GetParError(2)*0.197327*chimult;
  
  Double_t Ri = Rinv = TMath::Sqrt((Ro*Ro + Rs*Rs + Rl*Rl)/3.0);
  Double_t RiE = RinvE = (RoE*Ro + RsE*Rs + RlE*Rl)/(3*Ri);

  cout << endl;
  cout << "Rinv " << Ri << " +/- " << RiE << endl;

  Ro = fitter->GetParameter(11)*0.197327;
  RoE = fitter->GetParError(11)*0.197327*chimult;
  Rs = fitter->GetParameter(12)*0.197327;
  RsE = fitter->GetParError(12)*0.197327*chimult;
  Rl = fitter->GetParameter(13)*0.197327;
  RlE = fitter->GetParError(13)*0.197327*chimult;

  Ri = TMath::Sqrt((Ro*Ro + Rs*Rs + Rl*Rl)/3.0);
  RiE = (RoE*Ro + RsE*Rs + RlE*Rl)/(3*Ri);

  cout << endl;
  cout << "RinvS " << Ri << " +/- " << RiE << endl;

  Double_t *covmatrix;
  covmatrix = fitter->GetCovarianceMatrix();
  printf("\tRout\t\tRside\t\tRlong\t\tLambda\t\tNorm\t\tDev20\t\tDev22\n");
  for (int itery=0; itery<9; itery++) {
    if ((itery == 4) || (itery == 5)) continue;
    if (itery == 0) printf("Rout\t");
    if (itery == 1) printf("Rside\t");
    if (itery == 2) printf("Rlong\t");
    if (itery == 3) printf("Lambda\t");
    if (itery == 6) printf("Norm\t");
    if (itery == 7) printf("Dev20\t");
    if (itery == 8) printf("Dev22\t");
    for (int iterx=0; iterx<9; iterx++) {
      if ((iterx == 4) || (iterx == 5)) continue;
      printf("%.6lf\t", fitter->GetCovarianceMatrixElement((itery>3)?itery-2:itery, (iterx>3)?iterx-2:iterx));
    }
    printf("\n");
  }
  printf("\n");
  printf("\tRout\t\tRside\t\tRlong\t\tLambda\t\tNorm\t\tDev20\t\tDev22\n");
  for (int itery=0; itery<9; itery++) {
    if ((itery == 4) || (itery == 5)) continue;
    if (itery == 0) printf("Rout\t");
    if (itery == 1) printf("Rside\t");
    if (itery == 2) printf("Rlong\t");
    if (itery == 3) printf("Lambda\t");
    if (itery == 6) printf("Norm\t");
    if (itery == 7) printf("Dev20\t");
    if (itery == 8) printf("Dev22\t");
    for (int iterx=0; iterx<9; iterx++) {
      if (iterx > itery) continue;
      if ((iterx == 4) || (iterx == 5)) continue;
      printf("%.6lf\t", 
	     fitter->GetCovarianceMatrixElement((itery>3)?itery-2:itery, (iterx>3)?iterx-2:iterx)/
	     TMath::Sqrt(fitter->GetCovarianceMatrixElement((itery>3)?itery-2:itery, (itery>3)?itery-2:itery)*
			 fitter->GetCovarianceMatrixElement((iterx>3)?iterx-2:iterx, (iterx>3)?iterx-2:iterx)));
    }
    printf("\n");
  }

#ifndef _NOSAVE_
  Double_t xyf[201];
  Double_t yk00L[201];
  Double_t yk20L[201];
  Double_t yk22L[201];

  Double_t yk00S[201];
  Double_t yk20S[201];
  Double_t yk22S[201];

  Double_t yk00A[201];
  Double_t yk20A[201];
  Double_t yk22A[201];
  
  int ip;
  double ix;
  double fend = c00->GetXaxis()->GetBinCenter(bine);

  Double_t ng = fitter->GetParameter(6);
  Double_t n0 = fitter->GetParameter(3);
  Double_t nl = fitter->GetParameter(7);
  Double_t nt = fitter->GetParameter(8);
  Double_t kavval, tv00, tv20, tv22, c00val;
  
  Double_t xx[1];
  Double_t part[14];
  part[0] = fitter->GetParameter(3);
  part[1] = fitter->GetParameter(0);
  part[2] = fitter->GetParameter(1);
  part[3] = fitter->GetParameter(2);
  part[4] = fitter->GetParameter(14);
  part[5] = fitter->GetParameter(11);
  part[6] = fitter->GetParameter(12);
  part[7] = fitter->GetParameter(13);
  part[8] = fitter->GetParameter(4);
  part[9] = fitter->GetParameter(5);
  part[10] = fitter->GetParameter(7);
  part[11] = fitter->GetParameter(8);
  part[12] = fitter->GetParameter(9);
  part[13] = fitter->GetParameter(10);
//   part[14] = fitter->GetParameter(15);
//   part[15] = fitter->GetParameter(16);
//   part[16] = fitter->GetParameter(17);
//   part[17] = fitter->GetParameter(18);
//   part[18] = fitter->GetParameter(19);


  for (ix=0.001,  ip=0; ix<fend; ix+=0.02, ip++) {
    //    cout << ix << " " << f1->Eval(ix) << endl;
    xx[0] = ix;

    xyf[ip] = ix;

    funfull(xx, part, &tv00, &tv20, &tv22);

    yk00A[ip] = tv00*fitter->GetParameter(6);
    yk20A[ip] = (tv20+fitter->GetParameter(15)*exp(-fitter->GetParameter(16)*fitter->GetParameter(16)*ix*ix) + fitter->GetParameter(17)*exp(-(ix-fitter->GetParameter(18))*(ix-fitter->GetParameter(18))*fitter->GetParameter(19)*fitter->GetParameter(19)))*fitter->GetParameter(6);

    //    cout << ix << "   " << tv20 << "  " << (fitter->GetParameter(15)*exp(-fitter->GetParameter(16)*fitter->GetParameter(16)*ix*ix) + fitter->GetParameter(17)*exp(-(ix-fitter->GetParameter(18))*(ix-fitter->GetParameter(18))*fitter->GetParameter(19)*fitter->GetParameter(19))) << "     " << yk20A[ip] << endl;

    yk22A[ip] = tv22*fitter->GetParameter(6);
  }

  part[0] = 0.0;

  for (ix=0.001,  ip=0; ix<fend; ix+=0.02, ip++) {
    //    cout << ix << " " << f1->Eval(ix) << endl;
    xx[0] = ix;

    funfull(xx, part, &tv00, &tv20, &tv22);

    yk00S[ip] = tv00*fitter->GetParameter(6);
    yk20S[ip] = (tv20+(fitter->GetParameter(15)*exp(-fitter->GetParameter(16)*fitter->GetParameter(16)*ix*ix) + fitter->GetParameter(17)*exp(-(ix-fitter->GetParameter(18))*(ix-fitter->GetParameter(18))*fitter->GetParameter(19)*fitter->GetParameter(19))))*fitter->GetParameter(6);
    yk22S[ip] = tv22*fitter->GetParameter(6);

  }

  part[0] = fitter->GetParameter(3);
  part[4] = 0.0;

  for (ix=0.001,  ip=0; ix<fend; ix+=0.02, ip++) {
    //    cout << ix << " " << f1->Eval(ix) << endl;
    xx[0] = ix;

    funfull(xx, part, &tv00, &tv20, &tv22);

    yk00L[ip] = tv00*fitter->GetParameter(6);
    yk20L[ip] = tv20*fitter->GetParameter(6);
    yk22L[ip] = tv22*fitter->GetParameter(6);
  }

  TGraph *gr00L = new TGraph(ip,xyf, yk00L);
  TGraph *gr20L = new TGraph(ip,xyf, yk20L);
  TGraph *gr22L = new TGraph(ip,xyf, yk22L);

  TGraph *gr00S = new TGraph(ip,xyf, yk00S);
  TGraph *gr20S = new TGraph(ip,xyf, yk20S);
  TGraph *gr22S = new TGraph(ip,xyf, yk22S);

  TGraph *gr00A = new TGraph(ip,xyf, yk00A);
  TGraph *gr20A = new TGraph(ip,xyf, yk20A);
  TGraph *gr22A = new TGraph(ip,xyf, yk22A);

  gr00L->SetLineColor(kBlue);
  gr20L->SetLineColor(kBlue);
  gr22L->SetLineColor(kBlue);

  gr00S->SetLineColor(kRed);
  gr20S->SetLineColor(kRed);
  gr22S->SetLineColor(kRed);


  TCanvas *canfit = new TCanvas ("canfit","canfit",600,1200);
  canfit->SetFillColor(0);
  canfit->SetFillStyle(40000);
  canfit->Divide(1,3,0.0001,0.0001);
  canfit->cd(1);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(40000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.13);
  c00->SetTitle(";q_{LCMS} [GeV/c];C_{0}^{0}");
  c00->GetXaxis()->SetTitleSize(0.06);
  c00->GetXaxis()->SetLabelSize(0.06);
  c00->GetYaxis()->SetTitleSize(0.06);
  c00->GetYaxis()->SetLabelSize(0.06);
  c00->SetMaximum(2.13);
  //  c00->GetXaxis()->SetTitle("q_{LCMS} [GeV/c]");
  c00->Draw();
  //  f1->Draw("");
  gr00L->Draw("CP");
  gr00S->Draw("CP");
  gr00A->Draw("CP");

  canfit->cd(2);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(40000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.13);
  c20->SetTitle(";q_{LCMS} [GeV/c];C_{2}^{0}");
  c20->GetXaxis()->SetTitleSize(0.06);
  c20->GetXaxis()->SetLabelSize(0.06);
  c20->GetYaxis()->SetTitleSize(0.06);
  c20->GetYaxis()->SetLabelSize(0.06);
  c20->SetMinimum(-0.149);
  c20->SetMaximum(0.149);
  c20->Draw();
  //  f20->Draw("");
  gr20L->Draw("CP");
  gr20S->Draw("CP");
  gr20A->Draw("CP");

  canfit->cd(3);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(40000);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.13);
  c22->SetTitle(";q_{LCMS} [GeV/c];C_{2}^{2}");
  c22->GetXaxis()->SetTitleSize(0.06);
  c22->GetXaxis()->SetLabelSize(0.06);
  c22->GetYaxis()->SetTitleSize(0.06);
  c22->GetYaxis()->SetLabelSize(0.06);
  c22->SetMinimum(-0.149);
  c22->SetMaximum(0.149);
  c22->Draw();
  //  f22->Draw("");
  gr22L->Draw("CP");
  gr22S->Draw("CP");
  gr22A->Draw("CP");

  canfit->SaveAs("canfit.root");
  canfit->SaveAs("canfit.png");

  TFile *ofit = new TFile("ofit.root","RECREATE");
  ofit->cd();
  c00->Write();
  gr00L->SetName("gr00L");
  gr00L->Write();
  gr00S->SetName("gr00S");
  gr00S->Write();
  gr00A->SetName("gr00A");
  gr00A->Write();
  c20->Write();
  gr20L->SetName("gr20L");
  gr20L->Write();
  gr20S->SetName("gr20S");
  gr20S->Write();
  gr20A->SetName("gr20A");
  gr20A->Write();
  c22->Write();
  gr22L->SetName("gr22L");
  gr22L->Write();
  gr22S->SetName("gr22S");
  gr22S->Write();
  gr22A->SetName("gr22A");
  gr22A->Write();
#endif
}

bool fitshanalyticaaabackshdircovcoulpars(const char *filname, 
  Double_t &Rout, Double_t &Rside, Double_t &Rlong, Double_t &Rinv, Double_t &lambda,
   Double_t &dRout, Double_t &dRside, Double_t &dRlong, Double_t &dRinv, Double_t &dlambda) 
{
  TFile *infilek = new TFile("ffcomp.root");
  if(infilek->IsZombie())
  {
    cout << "could not load ffcomp.root" << endl;
    return false;
  }
  grkcoul = (TGraph *) infilek->Get("KCoulomb");

  ifstream *inf;

  inf = new ifstream("fitpar.in");
  if(inf->fail())
  {
    cout << "loading parameter file failed" << endl;
    return false;
  }

  char numname[100];
  char denname[100];

  for (int iter=0; iter<20; iter++) {
    GetPar(inf, &parsg[iter], &dofix[iter], &parmin[iter], &parmax[iter]);
    cout << "Got fix " << iter << " " << dofix[iter] << endl;
  }
  
  (*inf) >> numname;
  (*inf) >> fitb;
  (*inf) >> fite;
  (*inf) >> isprf;
  (*inf) >> funout;
  (*inf) >> funside;
  (*inf) >> funlong;

  TFile *infile = new TFile(filname);
  infile->cd();

  fitshanalyticreal(numname, Rout, Rside, Rlong, Rinv, lambda, dRout, dRside, dRlong, dRinv, dlambda);
  return true;
}