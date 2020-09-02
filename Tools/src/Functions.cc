#include "../interface/Functions.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

TMatrixDSym compMomentumTensor(double r, std::vector<TLorentzVector> & inputVectors_)
{
  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();

  if ( inputVectors_.size() < 2 ){
    return momentumTensor;
  }

  // fill momentumTensor from inputVectors
  double norm = 1.;
  for ( int i = 0; i < (int)inputVectors_.size(); ++i ){
    double p2 = inputVectors_[i].Dot(inputVectors_[i]);
    double pR = ( r == 2. ) ? p2 : TMath::Power(p2, 0.5*r);
    norm += pR;
    double pRminus2 = ( r == 2. ) ? 1. : TMath::Power(p2, 0.5*r - 1.);
    momentumTensor(0,0) += pRminus2*inputVectors_[i].X()*inputVectors_[i].X();
    momentumTensor(0,1) += pRminus2*inputVectors_[i].X()*inputVectors_[i].Y();
    momentumTensor(0,2) += pRminus2*inputVectors_[i].X()*inputVectors_[i].Z();
    momentumTensor(1,0) += pRminus2*inputVectors_[i].Y()*inputVectors_[i].X();
    momentumTensor(1,1) += pRminus2*inputVectors_[i].Y()*inputVectors_[i].Y();
    momentumTensor(1,2) += pRminus2*inputVectors_[i].Y()*inputVectors_[i].Z();
    momentumTensor(2,0) += pRminus2*inputVectors_[i].Z()*inputVectors_[i].X();
    momentumTensor(2,1) += pRminus2*inputVectors_[i].Z()*inputVectors_[i].Y();
    momentumTensor(2,2) += pRminus2*inputVectors_[i].Z()*inputVectors_[i].Z();
  }


  // return momentumTensor normalized to determinant 1
  return (1./norm)*momentumTensor;
}

/// helper function to fill the 3 dimensional vector of eigen-values;
/// the largest (smallest) eigen-value is stored at index position 0 (2)
TVectorD compEigenValues(double r, std::vector<TLorentzVector> & inputVectors_)
{
  TVectorD eigenValues(3);
  TMatrixDSym myTensor = compMomentumTensor(r, inputVectors_);
  if( myTensor.IsSymmetric() ){
    if( myTensor.NonZeros() != 0 ) myTensor.EigenVectors(eigenValues);
  }

  // CV: TMatrixDSym::EigenVectors returns eigen-values and eigen-vectors
  //     ordered by descending eigen-values, so no need to do any sorting here...
  //std::cout << "eigenValues(0) = " << eigenValues(0) << ","
  //      << " eigenValues(1) = " << eigenValues(1) << ","
  //      << " eigenValues(2) = " << eigenValues(2) << std::endl;

  return eigenValues;
}

/// 1.5*(q1+q2) where 0<=q1<=q2<=q3 are the eigenvalues of the momentum tensor sum{p_j[a]*p_j[b]}/sum{p_j**2} 
/// normalized to 1. Return values are 1 for spherical, 3/4 for plane and 0 for linear events


float sphericity(double r,  std::vector<TLorentzVector> & inputVectors_)
{
  TVectorD eigenValues = compEigenValues(r, inputVectors_);
  double tmp1 = eigenValues(1);
  double tmp2 = eigenValues(2);
  return float(1.5*(eigenValues(1) + eigenValues(2)));
}


double clean(double x, double d) {
  return TMath::Finite(x) ? x : d;
}

double bound(double val, double low, double high) {
  return TMath::Max(low,TMath::Min(high,val));
}

int dsign(double x) {
    return sign(x);
}

double Mxx(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2) {
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
    return (v1+v2).M();
}

double MT(double pt1, double phi1, double pt2, double phi2)
{
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1,0,phi1,0);
    v2.SetPtEtaPhiM(pt2,0,phi2,0);
    return (v1+v2).M();
}

double SignedDeltaPhi(double phi1, double phi2) {
    double dPhi = phi1-phi2;
    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;
    return dPhi;
}

double DeltaPhi(double phi1, double phi2) {
    double dPhi = phi1-phi2;
    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;
    return fabs(dPhi);
}

double DeltaR2(double eta1, double phi1, double eta2, double phi2) {
    float dEta2 = (eta1-eta2); dEta2 *= dEta2;
    float dPhi = SignedDeltaPhi(phi1,phi2);
    return dEta2 + dPhi*dPhi;
}

double ExpErf(double x, double a, double b, double c) {
    double exp_ = TMath::Exp(c*x);
    double erf_ = TMath::Erf((x-a)/b);
    return exp_*(1+erf_)/2;
}

double Expnom(double x, double a, double b, double c){
   double exp_ = a * TMath::Exp(b*x)+c;
   return exp_;
}

double tightIDsf(int bit, float sf1, float sf2){
   double sf = 1;
   if((bit>>3)%2) sf = sf1;
   else sf = sf2;
   return sf;
}

double sfphonew(double eta){
  if(fabs(eta)<0.8) return 1.019;
  else if(fabs(eta)<1.4) return 1.015;
  else return 1;
}


double sfpho(double eta){
    double x[11]={-2.5,-2,-1.566,-1.444,-0.8,0,0.8,1.444,1.566,2,2.5};
    double y[10]={0.995122,1.05853,1,0.983432,0.967391,0.954933,0.929929,1,0.860262,1.00601};

    int idx=-1;
    for(int i=0; i<10; i++){
      if(eta>=x[i]&&eta<x[i+1]) idx=i;
    }

    if(idx>=0) return y[idx];
    else return 1;
}

double sfpho2018(double eta){
    double x[11]={-2.5,-2,-1.566,-1.444,-0.8,0,0.8,1.444,1.566,2,2.5};
    double y[10]={1.061,1.037,1,1.027,1.035,1.052,0.977,1,1.090,1.016};

    int idx=-1;
    for(int i=0; i<10; i++){
      if(eta>=x[i]&&eta<x[i+1]) idx=i;
    }

    if(idx>=0) return y[idx];
    else return 1;
}

double fr_pho(double pt){
  double x[18]={230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,100000};
  double y[17]={0.00495011,0.00471173,0.00348358,0.00411499,0.00424737,0.00466982,0.00517474,0.00533911,0.00800606,0.00811088,0.0079036,0.00839311,0.0113206,0.0115227,0.00903586,0.0103002,0.00774355};
   if(pt>230){
      double yy=1;
      for(int i=0;i<17; i++){
       if(pt>=x[i]&&pt<x[i+1]) yy=y[i];
      }
      return yy;
   }
  else return 1;

}

double photrigsf(double pt){
  double a1=0.335,b1=217.91,c1=0.065,d1=0.996;
  double a2=0.244,b2=212.34,c2=0.050,d2=1;
  double y1=c1+(d1-c1)/(1+exp(-a1*(pt-b1)));
  double y2=c2+(d2-c2)/(1+exp(-a2*(pt-b2)));
  return y1/y2;
}

double photrigsf2018(double pt){
  double a1=1.022,b1=218.39,c1=0.086,d1=0.999;
  double a2=0.301,b2=212.83,c2=0.062,d2=1.000;
  double y1=c1+(d1-c1)/(1+exp(-a1*(pt-b1)));
  double y2=c2+(d2-c2)/(1+exp(-a2*(pt-b2)));
  return y1/y2;
}


double metsf(double recoil){
  double x[27]={0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,500,600,700,800,900,1000};
  int bin=-1;

  double y1[26]={5.05247,2.44119,1.75749,0.625013,0.539698,0.526736,0.555387,0.642364,0.750472,0.84822,0.916332,0.959062,0.97798,0.988869,0.994051,0.994889,0.996896,0.99699,0.997758,0.995825,0.99643,0.997156,0.99562,0.999024,0.994479,0.996249};

  for(int i=0; i<26; i++){
    if(recoil>=x[i]&&recoil<=x[i+1]) bin=i;
  }
  
  if(bin>=0)  return y1[bin];
  else return 1;
}

double metsf2018(double recoil){
  double x[27]={0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,500,600,700,800,900,1000};
  int bin=-1;

  double y1[26]={1.91756,1.90197,2.28993,0.78825,0.52439,0.4962,0.540381,0.632407,0.74145,0.840749,0.911223,0.953317,0.97695,0.987622,0.993508,0.996202,0.996901,0.9963,0.997581,0.999614,0.999398,0.999141,0.998145,0.998598,0.999982,0.99824};

  for(int i=0; i<26; i++){
    if(recoil>=x[i]&&recoil<=x[i+1]) bin=i;
  }

  if(bin>=0)  return y1[bin];
  else return 1;
}

double sfwjets(double pt){
  double x[25]={150,170,200,230,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,10000};
  double y[24]={0.966481,0.958351,0.948513,0.938304,0.929394,0.919688,0.910425,0.900288,0.889122,0.877255,0.86653,0.856032,0.84536,0.834675,0.822288,0.810779,0.799689,0.788566,0.777413,0.765405,0.753711,0.742262,0.72925,0.716423};
   double kqcd=1.053*exp(-0.003163*pt)+0.746;
   if(pt>150){ 
      double yy=1;
      for(int i=0;i<24; i++){
       if(pt>=x[i]&&pt<x[i+1]) yy=y[i]*kqcd;
      }
      return yy;
   }
   else return 1;

}

double sfzjets(double pt){
  double x[25]={150,170,200,230,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,10000};
  double y[24]={0.970592,0.964424,0.956695,0.948747,0.941761,0.934246,0.927089,0.919181,0.909926,0.900911,0.892561,0.884353,0.8761,0.867687,0.858047,0.849014,0.840317,0.832017,0.823545,0.814596,0.806229,0.798038,0.789694,0.781163};
   double kqcd=1.434*exp(-0.00221*pt)+0.443;


   if(pt>150){ 
      double yy=1;
      for(int i=0;i<24; i++){
       if(pt>=x[i]&&pt<x[i+1]) yy=y[i]*kqcd;
      }
      return yy;
   }
   else return kqcd;

}

double sfzvv(double pt){
  double xqcd[48]={120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320,1360,1400,1440,1480,1520,1560,1600,1640,1680,1720,1760,1800,1840,1880,1920,1960,100000};
  double yqcd[47]={1.51607,1.42685,1.34526,1.27085,1.20313,1.14163,1.08587,1.03539,0.989713,0.948359,0.910398,0.875937,0.844316,0.815357,0.788256,0.761675,0.736038,0.713414,0.69041,0.668093,0.652968,0.636087,0.627602,0.612289,0.601274,0.587895,0.578841,0.568343,0.547456,0.547034,0.525202,0.503478,0.505757,0.490673,0.460639,0.466632,0.46423,0.446041,0.443809,0.443691,0.445554,0.449265,0.45469,0.461696,0.470149,0.479917,0.490865};

  double x[25]={150,170,200,230,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,10000};
  double y[24]={0.970592,0.964424,0.956695,0.948747,0.941761,0.934246,0.927089,0.919181,0.909926,0.900911,0.892561,0.884353,0.8761,0.867687,0.858047,0.849014,0.840317,0.832017,0.823545,0.814596,0.806229,0.798038,0.789694,0.781163};

   double kqcd=1;
   if(pt>120){
      for(int i=0;i<47; i++){
       if(pt>=xqcd[i]&&pt<xqcd[i+1]) kqcd=yqcd[i];
      }
   }


   if(pt>150){
      double yy=1;
      for(int i=0;i<24; i++){
       if(pt>=x[i]&&pt<x[i+1]) yy=y[i]*kqcd;
      }
      return yy;
   }
   else return yqcd[0]*y[0];
}

double sfgjets(double pt){

  double xqcd[25]={150,170,200,230,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,100000};
  double yqcd[24]={1.36708,1.39204,1.43241,1.39898,1.3592,1.32905,1.31791,1.29202,1.23844,1.21516,1.19233,1.17271,1.17239,1.1889,1.1316,1.10003,1.12668,1.09366,1.05034,1.01872,0.97818,0.965138,1.0759,1.17312};

  double x[25]={150,170,200,230,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,10000};
  double y[24]={0.993448,0.990594,0.987151,0.983324,0.980437,0.977336,0.974328,0.970768,0.968094,0.963285,0.959664,0.955693,0.951956,0.948179,0.94396,0.940431,0.936879,0.933292,0.929644,0.925767,0.922028,0.918239,0.914257,0.910161};
//   double kqcd=1.159*exp(-0.001944*pt)+1;
   if(pt>150){
      double yy=1;
      for(int i=0;i<24; i++){
       if(pt>=x[i]&&pt<x[i+1]) yy=y[i]*yqcd[i];
      }
      return yy;
   }
   else return yqcd[0]*y[0];
}

double combinedPhi(double pt1, double phi1, double pt2, double phi2)
{
  TVector2 p1, p2;
  p1.SetMagPhi(pt1,phi1); p2.SetMagPhi(pt2,phi2);
  TVector2 combined=p1+p2;
  return combined.Phi();
}

double combinedPt(double pt1, double phi1, double pt2, double phi2)
{
  TVector2 p1, p2;
  p1.SetMagPhi(pt1,phi1); p2.SetMagPhi(pt2,phi2);
  TVector2 combined=p1+p2;
  return combined.Mod();
}

double tauSF(double x, int i){
  double sf=1;
//2017
  if(i==1) sf=(x<=20)*0+ ( x > 20 && x <=25)*0.8890281+ ( x > 25 && x <=30)*0.9469252+ ( x > 30 && x <=35)*0.9380708+ ( x > 35 && x <=40)*0.9044377+ (x > 40 && x <= 500)*0.949084749706+ (x > 500 && x <= 1000)*(0.908802967687 + 0.0402817820196*(x/500.))+ (x > 1000)*(0.908802967687 + 0.0805635640393);
  else if(i==2) sf=(x<=20)*0+ ( x > 20 && x <=25)*0.6469781+ ( x > 25 && x <=30)*0.8071072+ ( x > 30 && x <=35)*0.8347448+ ( x > 35 && x <=40)*0.8107377+ (x > 40 && x <= 500)*0.861070663387+ (x > 500 && x <= 1000)*(0.908802967687 - 0.0477323042994*(x/500.))+ (x > 1000)*(0.908802967687 - 0.0954646085989);
  else sf=(x<=20)*0+ ( x > 20 && x <=25)*0.7680031+ ( x > 25 && x <=30)*0.8770162+ ( x > 30 && x <=35)*0.8864078+ ( x > 35 && x <=40)*0.8575877+ (x > 40)*0.908802967687;
//2018
/*
  if(i==1) sf=(x<=20)*0+ ( x > 20 && x <=25)*1.0660781+ ( x > 25 && x <=30)*1.0602211+ ( x > 30 && x <=35)*0.9593934+ ( x > 35 && x <=40)*0.9217598+ (x > 40 && x <= 500)*1.01133836263+ (x > 500 && x <= 1000)*(0.977607152447 + 0.0337312101811*(x/500.))+ (x > 1000)*(0.977607152447 + 0.0674624203623);
  else if(i==2) sf=(x<=20)*0+ ( x > 20 && x <=25)*0.8664641+ ( x > 25 && x <=30)*0.9080931+ ( x > 30 && x <=35)*0.8391794+ ( x > 35 && x <=40)*0.8256058+ (x > 40 && x <= 500)*0.915656669816+ (x > 500 && x <= 1000)*(0.977607152447 - 0.0619504826307*(x/500.))+ (x > 1000)*(0.977607152447 - 0.123900965261);
  else sf=(x<=20)*0+ ( x > 20 && x <=25)*0.9662711+ ( x > 25 && x <=30)*0.9841571+ ( x > 30 && x <=35)*0.8992864+ ( x > 35 && x <=40)*0.8736828+ (x > 40)*0.977607152447;
*/
  return sf;
}

double eleunc(double pt, double eta,int id){
  double x_r[6]={10,20,45,75,100,10000};
  double y_r[13]={-2.5,-2,-1.566,-1.444,-1,-0.5,0,0.5,1,1.444,1.566,2,2.5};
  double x_id[7]={10,20,35,50,100,200,10000};
  double y_id[11]={-2.5,-2,-1.566,-1.444,-0.8,0,0.8,1.444,1.566,2,2.5};
//2017
  double uncreco[5][12]={{0.0606322,0.0203413,0.0307087,0.0846001,0.0562095,0.0562095,0.0562095,0.0562095,0.0846001,0.0307087,0.0203413,0.0606322},
                         {0.00775375,0.00265162,0.0146941,0.00518345,0.00490105,0.00702552,0.00702574,0.00492188,0.00517995,0.0145467,0.00265709,0.00773705},
                         {0.0037683,0.00393991,0.0143485,0.00448452,0.00269231,0.0026979,0.00270351,0.0026951,0.00447612,0.0144763,0.00393576,0.00376805},
                         {0.0150853,0.00782194,0.0406877,0.00790241,0.00250203,0.00395202,0.00395202,0.00250203,0.00790241,0.0406877,0.00782194,0.0150853},
                         {0.00972415,0.0048987,0.0486547,0.00637478,0.00500406,0.001762,0.001762,0.00500406,0.00637478,0.0486547,0.0048987,0.00972415}};
  double uncveto[6][10]={{0.0148581,0.027122,1,0.0716767,0.0423557,0.0428655,0.0726424,1,0.0276708,0.0148115},
                         {0.0201298,0.0198773,1,0.0364958,0.0235296,0.023504,0.0363334,1,0.0198554,0.020174},
                         {0.0059409,0.00632067,1,0.00669325,0.00589688,0.00589053,0.00667918,1,0.00633405,0.00596},
                         {0.0100874,0.00618614,1,0.00545523,0.00483512,0.00482492,0.00542683,1,0.00620557,0.0101054},
                         {0.018279,0.042443,1,0.0297531,0.0224298,0.0222684,0.0298456,1,0.0423548,0.0184289},
                         {0.0492527,0.0842523,1,0.0622674,0.0135541,0.0134638,0.0613385,1,0.0867114,0.0480777}};

  double unctight[6][10]={{0.0879713,0.0606584,1,0.0646497,0.049341,0.0501492,0.0662741,1,0.0601476,0.0869647},
                          {0.0288558,0.0316227,1,0.0368784,0.0233786,0.0235205,0.036817,1,0.0316268,0.0283302},
                          {0.0155378,0.0164442,1,0.0180374,0.0162448,0.0162191,0.0179363,1,0.0165303,0.0154793},
                          {0.0205271,0.0111227,1,0.00994803,0.0162164,0.0161315,0.00997226,1,0.0112307,0.020281},
                          {0.0209164,0.0585176,1,0.0261115,0.0426125,0.0425823,0.0264055,1,0.0574792,0.0215128},
                          {0.112117,0.121532,1,0.0598762,0.0343347,0.0343132,0.0596292,1,0.124078,0.097496}};

//2018
/*
  double uncreco[5][12]={{0.0147606,0.0180825,0.0507138,0.020546,0.00736356,0.0215549,0.0215549,0.00736356,0.020546,0.0507138,0.0180825,0.0147606},
                         {0.00147782,0.00147238,0.00263521,0.00149019,0.00146929,0.00147696,0.00146627,0.00145871,0.00149023,0.00267015,0.00147083,0.00148233},
                         {0.00146706,0.00146022,0.0033705,0.00147466,0.00146396,0.001467,0.00145646,0.00145197,0.00147469,0.00353802,0.00145719,0.00146699},
                         {0.0113314,0.0104506,0.0283021,0.00682089,0.00522452,0.00455881,0.00447623,0.00493653,0.00680673,0.0275973,0.00957309,0.011167},
                         {0.0121418,0.00870061,0.0274722,0.00683116,0.00452188,0.00452188,0.00452188,0.00452188,0.00683116,0.0274722,0.00870061,0.0121418}};
  double uncveto[6][10]={{0.0123986,0.0207671,1,0.0181053,0.0222039,0.0215699,0.0183423,1,0.021196,0.0126006},
                         {0.0201435,0.0219587,1,0.0332007,0.0210758,0.0209334,0.0335497,1,0.0218629,0.0203266},
                         {0.00666115,0.00246603,1,0.00231481,0.00216437,0.00216208,0.00230757,1,0.00246346,0.00667544},
                         {0.0135708,0.00300606,1,0.00342604,0.0040573,0.00405725,0.00343321,1,0.00300604,0.0136142},
                         {0.0142779,0.0157512,1,0.0116066,0.00649975,0.00651315,0.0116544,1,0.0156806,0.014561},
                         {0.0612724,0.0283079,1,0.0184826,0.0193021,0.0192937,0.0180257,1,0.0272889,0.0588555}};
  double unctight[6][10]={{0.0383151,0.0466478,1,0.0420617,0.024811,0.0235889,0.0419376,1,0.0482084,0.0363606},
                          {0.0281661,0.0301148,1,0.0255931,0.0176143,0.0173567,0.0253841,1,0.0299706,0.0286798},
                          {0.016482,0.00394214,1,0.00619856,0.00345462,0.00343137,0.0061571,1,0.00398869,0.0165889},
                          {0.025059,0.00756358,1,0.00699055,0.0120398,0.0119654,0.0069297,1,0.00769257,0.0256944},
                          {0.0303811,0.0190396,1,0.0297827,0.0168526,0.0163624,0.0289321,1,0.0192371,0.0307287},
                          {0.0843876,0.0515224,1,0.0463371,0.0307746,0.0294709,0.0445142,1,0.0505062,0.0773356}};
*/
  double unc1=0,unc2=0;
  int indx,indy,indx2,indy2;
  if(!(pt>=10&&pt<10000&&fabs(eta)<2.5)) return 0;
  for(int i=0;i<5;i++){
   if(pt>=x_r[i]&&pt<x_r[i+1]) indx=i;
  }
  for(int i=0;i<6;i++){
   if(pt>=x_id[i]&&pt<x_id[i+1]) indx2=i;
  }
  for(int i=0;i<12;i++){
   if(eta>=y_r[i]&&eta<y_r[i+1]) indy=i;
  }
  for(int i=0;i<10;i++){
   if(eta>=y_id[i]&&eta<y_id[i+1]) indy2=i;
  }
  unc1=uncreco[indx][indy];
  if(id==0) unc2=uncveto[indx2][indy2];
  else unc2=unctight[indx2][indy2];
  return sqrt(unc1*unc1+unc2*unc2);
}

