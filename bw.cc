#include <complex>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TDecompBK.h"

using namespace std;

//#include "complexMatrix.h"
#include "bw.h"

double
blattWeisskopf(int L, double p)
{
  double z = pow(p / hbarc, 2);  // 1fm interaction radius, VES uses 5.2 GeV^-1
  double result;
  switch (L) {
  case 0:
    result = 1; break;
  case 1:
    result = sqrt(2*z/(z+1)); break;
  case 2:
    result = sqrt(13*z*z/(pow(z-3, 2) + 9*z)); break;
  case 3:
    result = sqrt(277*z*z*z/(z*pow(z-15, 2) + 9*pow(2*z-5, 2))); break;
  case 4:
    result = sqrt(12746*pow(z,4)/(pow(z*z-45*z+105,2) + 25*z*pow(2*z-21,2))); break;
  default:
    result = 0; break;
  }
  return result;
}


// This returns the square and takes the square as argument to avoid
// unnecessary complex variables.
double
blattWeisskopf2(int L, double p2)
{
  double z = p2 / hbarc / hbarc;  // 1fm interaction radius, VES uses 5.2 GeV^-1
  double result;
  switch (L) {
  case 0:
    result = 1; break;
  case 1:
    result = 2.*z/(z+1.); break;
  case 2:
    result = 13.*z*z/(pow(z-3., 2) + 9.*z); break;
  case 3:
    result = 277.*z*z*z/(z*pow(z-15., 2) + 9.*pow(2.*z-5., 2)); break;
  case 4:
    result = 12746.*pow(z,4)/(pow(z*z-45.*z+105.,2) + 25.*z*pow(2.*z-21.,2)); break;
  default:
    result = 0; break;
  }
  return result;
}

double
breakupMomentum(double s, double m1, double m2)
{
  if (s < pow(m1 + m2, 2))
    return 0;
  double result = .5*sqrt((s - pow(m1+m2, 2))*(s - pow(m1-m2,2)) / s);
  return result;
}

// Square of breakup momentum, can be negative.
double
breakupMomentum2(double s, double m1, double m2)
{
  return .25*(s - pow(m1+m2, 2))*(s - pow(m1-m2,2)) / s;
}

// complex version with analytic continuation below threshold
complex<double>
breakupMomentumComplex(double s, double m1, double m2) 
{
  double q2 = breakupMomentum2(s, m1, m2);
  double q  = sqrt(fabs(q2));
  if (q2 < 0)
    return complex<double>(0, q);
  else
    return complex<double>(q, 0);
} 

complex<double>
AMP_M(double s, double m1, double m2, double alpha1, double alpha2)
{    
  // parameter
  //matrix<std::complex<double> > T(2, 2);
  vector<TMatrixD> a(2, TMatrixD(2,2));
  vector<TMatrixD> c(5, TMatrixD(2,2));
  vector<TMatrixD> pol(3, TMatrixD(1,2));
  TMatrixD sP(1,2);
  int vesSheet(1);
  int kachaev(1);
  
  const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1

  a[0](0, 0) =  0.1131;  // AMP Table 1, M solution: f_2^2 // schmarrn, steht eigentlich fuer a11
  a[0](0, 1) =  0.0150;  // AMP Table 1, M solution: f_1^3 // a12
  a[0](1, 0) =  0.0150;  // AMP Table 1, M solution: f_1^3 // a12
  a[0](1, 1) = -0.3216;  // AMP Table 1, M solution: f_2^3 // a22
  a[1](0, 0) = f[0] * f[0];
  a[1](0, 1) = f[0] * f[1];
  a[1](1, 0) = f[1] * f[0];
  a[1](1, 1) = f[1] * f[1];

  c[0](0, 0) =  0.0337;                // AMP Table 1, M solution: c_11^0
  c[1](0, 0) = -0.3185;                // AMP Table 1, M solution: c_11^1
  c[2](0, 0) = -0.0942;                // AMP Table 1, M solution: c_11^2
  c[3](0, 0) = -0.5927;                // AMP Table 1, M solution: c_11^3
  c[4](0, 0) =  0.1957;                // AMP Table 1, M solution: c_11^4
  c[0](0, 1) = c[0](1, 0) = -0.2826;  // AMP Table 1, M solution: c_12^0
  c[1](0, 1) = c[1](1, 0) =  0.0918;  // AMP Table 1, M solution: c_12^1
  c[2](0, 1) = c[2](1, 0) =  0.1669;  // AMP Table 1, M solution: c_12^2
  c[3](0, 1) = c[3](1, 0) = -0.2082;  // AMP Table 1, M solution: c_12^3
  c[4](0, 1) = c[4](1, 0) = -0.1386;  // AMP Table 1, M solution: c_12^4
  c[0](1, 1) =  0.3010;                // AMP Table 1, M solution: c_22^0
  c[1](1, 1) = -0.5140;                // AMP Table 1, M solution: c_22^1
  c[2](1, 1) =  0.1176;                // AMP Table 1, M solution: c_22^2
  c[3](1, 1) =  0.5204;                // AMP Table 1, M solution: c_22^3
  c[4](1, 1) = -0.3977;                // AMP Table 1, M solution: c_22^4

  pol[0](0, 0) =  0.1393;              // AMP Table 1, M solution: a_1^0
  pol[1](0, 0) = -0.02775;             // AMP Table 1, M solution: a_1^1
  pol[2](0, 0) =  0.3952;              // AMP Table 1, M solution: a_1^2
  pol[0](0, 1) =  3.241;               // AMP Table 1, M solution: a_2^0
  pol[1](0, 1) = -3.432;               // AMP Table 1, M solution: a_2^1
  pol[2](0, 1) =  1.141;               // AMP Table 1, M solution: a_2^2
  
  sP(0, 0) = -0.0074;  // AMP Table 1, M solution: s_0
  sP(0, 1) =  0.9828;  // AMP Table 1, M solution: s_1
  
  if (kachaev){
    // change parameters according to Kachaev's prescription                                                                                            
    c[4](0, 0) = 0; // was 0.1957;                                                                                                                     
    c[4](1, 1) = 0; // was -0.3977;                                                                                                                    
    a[0](0, 1) = 0; // was 0.0150                                                                                                                      
    a[0](1, 0) = 0; // was 0.0150                                                                                                                      
    // a[1] are the f's from the AMP paper                                                                                                             
    a[1](0, 0) = 0;
    a[1](0, 1) = 0;
    a[1](1, 0) = 0;
    a[1](1, 1) = 0;
  }
  
  const complex<double> imag(0, 1);
  
  double m = sqrt(s);
  if (fabs(s - sP(0, 1)) < 1e-6) {
    m += 1e-6;
    s = m * m;
  }

  double mKm = (mK+mKs)/2.;

  const complex<double> qPiPi   = breakupMomentumComplex(s, mPi, mPi);
  const complex<double> qPi0Pi0 = breakupMomentumComplex(s, mPi0, mPi0);
  const complex<double> qKK     = breakupMomentumComplex(s, mK, mK);
  const complex<double> qK0K0   = breakupMomentumComplex(s, mKs, mKs);
  complex<double>       qKmKm   = breakupMomentumComplex(s, mKm, mKm );
  
  vector<complex<double> > rho(2);
  if ( vesSheet ) {
    if (qKmKm.imag() > 0)
      qKmKm *= -1;
    rho[0] = (2. * qPiPi) / m;
    rho[1] = (2. * qKmKm) / m;
  } else {
    rho[0] = ((2. * qPiPi) / m + (2. * qPi0Pi0) / m) / 2.;
    rho[1] = ((2. * qKK)   / m + (2. * qK0K0)   / m) / 2.;
  }
  
  const double scale = (s / (4 * mKm * mKm)) - 1;
  
  vector<complex<double> > M(4);
  for (unsigned int i = 0; i < 4; ++i) {
    // translate 1x4 into 2x2
    int x = i/2;
    int y = i%2;
    
    for (unsigned int j = 0; j < 2; ++j) {
      const complex<double> fa = 1. / (s - sP(0, j));
      M[i] += fa * a[j](x,y);
    }
    
    for (unsigned int j = 0; j < 5; ++j) {
      const complex<double> sc = pow(scale, (int)j);
      M[i] += sc * c[j](x,y);
    }
  }
  
  for (unsigned int i = 0; i < 4; ++i) {
    int index = i/2;
    if (i == 0 || i ==3)
      M[i] -= imag*rho[index];
  }
  
  //cout << M[0] << " " << M[1] << " " << M[2] << " " << M[3] << endl;
  
  // modification: off-diagonal terms set to 0
  M[1] = 0;
  M[2] = 0;

  //vector<complex<double> > T(4);
  complex<double> det = M[0]*M[3] - M[1]*M[2];
  //T[0] =    M[3]/det;
  //T[1] = -1*M[2]/det;
  //T[2] = -1*M[1]/det;
  //T[3] =    M[0]/det;
  complex<double> alpha(alpha1, alpha2);
  const complex<double> amp = M[3]/det - alpha * M[2]/det;
  
  //if (_debug)
  //printDebug << name() << "(m = " << maxPrecision(mass) << " GeV) = "
  //<< maxPrecisionDouble(amp) << endl;
  
  return amp;
} 

complex<double>
AMP_K1_my(double s, double m1, double m2)
{    
  // parameter
  //matrix<std::complex<double> > T(2, 2);
  vector<TMatrixD> a(2, TMatrixD(2,2));
  vector<TMatrixD> c(5, TMatrixD(2,2));
  TMatrixD sP(1,2);
  
  sP(0, 0) = -0.0110;  // AMP Table 1, K1 solution: s_0
  sP(0, 1) =  0.9247;  // AMP Table 1, K1 solution: s_1

  const double f[2] = {-0.2242, 0.5829};  // AMP Table 1, K1 solution: f_1^1 and f_2^1

  c[0](0, 0) =  0.7347;                // AMP Table 1, K1 solution: c_11^0
  c[1](0, 0) = -0.5266;                // AMP Table 1, K1 solution: c_11^1
  c[2](0, 0) =  2.6151;                // AMP Table 1, K1 solution: c_11^2
  c[3](0, 0) = -1.7747;                // AMP Table 1, K1 solution: c_11^3
  c[4](0, 0) =  0.8031;                // AMP Table 1, K1 solution: c_11^4
  c[0](0, 1) = c[0](1, 0) = -3.2762;  // AMP Table 1, K1 solution: c_12^0
  c[1](0, 1) = c[1](1, 0) = -0.6662;  // AMP Table 1, K1 solution: c_12^1
  c[2](0, 1) = c[2](1, 0) =  0.8778;  // AMP Table 1, K1 solution: c_12^2
  c[3](0, 1) = c[3](1, 0) = -2.1190;  // AMP Table 1, K1 solution: c_12^3
  c[4](0, 1) = c[4](1, 0) =  0.2319;  // AMP Table 1, K1 solution: c_12^4
  c[0](1, 1) = -2.6785;                // AMP Table 1, K1 solution: c_22^0
  c[1](1, 1) =  7.9951;                // AMP Table 1, K1 solution: c_22^1
  c[2](1, 1) =  5.5763;                // AMP Table 1, K1 solution: c_22^2
  c[3](1, 1) = -1.4956;                // AMP Table 1, K1 solution: c_22^3
  c[4](1, 1) =  0     ;                // AMP Table 1, K1 solution: c_22^4
  
  const complex<double> imag(0, 1);
  
  double m = sqrt(s);

  double mKm = (mK+mKs)/2.;
  
  const complex<double> qPiPi   = breakupMomentumComplex(s, mPi, mPi);
  const complex<double> qPi0Pi0 = breakupMomentumComplex(s, mPi0, mPi0);
  const complex<double> qKK     = breakupMomentumComplex(s, mK, mK);
  const complex<double> qK0K0   = breakupMomentumComplex(s, mKs, mKs);
  //complex<double>       qKmKm   = breakupMomentumComplex(s, mKm, mKm );
  
  vector<complex<double> > rho(2);
  rho[0] = ((2. * qPiPi) / m + (2. * qPi0Pi0) / m) / 2.;
  rho[1] = ((2. * qKK)   / m + (2. * qK0K0)   / m) / 2.;

  const double scale = (s / (4 * mKm * mKm)) - 1;
  
  vector<complex<double> > K(4);
  //vector<complex<double> > Khat(4);
  for (unsigned int i = 0; i < 4; ++i) {
    // translate 1x4 into 2x2
    int x = i/2;
    int y = i%2;
    
    K[i] = 1;//( s - sP(0,0) ) / (4*mK*mK); //wrong in paper, comes later
    K[i] /= ( sP(0,1) - s ) * ( sP(0,1) - sP(0,0) );
    K[i] = f[x] * f[y];
    
    //cout << K[0] << " " << K[1] << " " << K[2] << " " << K[3] << endl;
    
    for (unsigned int j = 0; j < 5; ++j) {
      const complex<double> sc = pow(scale, (int)j);
      K[i] += sc * c[j](x,y);
    }

    K[i] *= ( s - sP(0,0) ) / (4*mK*mK);

    //Khat[i] = K[i]/( s - sP(0,0) );
  }

  // invert
  vector<complex<double> > M(4);
  complex<double> det_K = (K[0]*K[3] - K[1]*K[2]);
  det_K = 1. / det_K;
  M[0] += det_K * K[3];
  M[1] -= det_K * K[1];
  M[2] -= det_K * K[2];
  M[3] += det_K * K[0];

  for (unsigned int i = 0; i < 4; ++i) {
    int index = i/2;
    if (i == 0 || i ==3)
      M[i] -= imag*rho[index];
  }
  
  // invert
  vector<complex<double> > T(4);
  complex<double> det_M = (M[0]*M[3] - M[1]*M[2]);
  det_M = 1. / det_M;
  T[0] += det_M * M[3];
  T[1] -= det_M * M[1];
  T[2] -= det_M * M[2];
  T[3] += det_M * M[0];

  // modification: off-diagonal terms set to 0
  //T[1] = 0;
  //T[2] = 0;
  
  const complex<double> amp = T[0];
  
  return amp;
} 

complex<double>
AMP_Kred(double s, double m1, double m2, double alpha1, double alpha2)
{    
  const complex<double> zi(0, 1);
  const double pi = 3.1415926;
  //double pie = pi/180.;
  //complex<double> zie = zi*pie;
  complex<double> zieps = zi*0.00001;
  //double pip = atan(1.)*4;
  double empic2 = mPi*mPi*4.;
  //double empi02 = mPi0*mPi0*4.;
  //double emp = 0.137;        // mPim???
  //double emrho2 = 0.77*0.77; // mRho???
  double emK02 = mKs*mKs*4.;
  double emKc2 = mK*mK*4.;
  //double emK2 = 0.4956*0.4956*4.; // mKm???
  //double ehadbot = 0.28;
  //double ehadtop = 1.6;

  // parameter
  vector<TMatrixD> c(3, TMatrixD(2,2));
  TMatrixD f(2,2);
  
  double s0 = 0.41*empic2/4.;             // S0
  const double S[3] = {s0, 0.26261, 1.0811};  // S0, S1, S2
  //const double S[3] = {s0, 0.26261, alpha1};  // S0, S1, S2

  f(0,0) =  0.38949; // F11
  f(0,1) =  0.33961; // F12
  f(1,0) =  0.24150; // F21
  f(1,1) = -0.78538; // F22
    
  c[0](0, 0) =  0.14760;                 // C110
  c[1](0, 0) =  0.062181;                // C111
  c[2](0, 0) =  0.029465;                // C112
  c[0](0, 1) = c[0](1, 0) =  0.10914;    // C120
  c[1](0, 1) = c[1](1, 0) = -0.17912;    // C121
  c[2](0, 1) = c[2](1, 0) =  0.10758;    // C122
  c[0](1, 1) = -0.27253;                 // C220
  c[1](1, 1) =  0.79442;                 // C221
  c[2](1, 1) = -0.49529;                 // C222
  
  
  complex<double> zs = s + zieps;
  complex<double> zrho1 = sqrt(1. - empic2/zs); // does that work??
  //double rho1 = sqrt(1. - empic2/s);
  //double rho2;
  complex<double> zrho2 = sqrt(1. - emK02/zs)*0.5 + sqrt(1. - emKc2/zs)*0.5;
  complex<double> zf1 = zrho1/pi * log((zrho1+1.)/(zrho1-1.));
  complex<double> zf2 = zrho2/pi * log((zrho2+1.)/(zrho2-1.));
  /*if (s > emK02) 
    rho2 = abs(zrho2);
  else
  rho2 = 0;*/
  double q2 = s/emKc2 - 1;
  double q4 = q2*q2;
  double q[3] = {1,q2,q4};

  vector<double> AQL(4);
  for (unsigned int i = 0; i < 4; ++i) {
    // translate 1x4 into 2x2
    int x = i/2;
    int y = i%2;
    
    AQL[i] = 0;
    for (unsigned int j = 0; j < 3; ++j)
      AQL[i] += q[j] * c[j](x,y);
  }

  vector<double> APL1(4);
  APL1[0] = f(0,0)*f(0,0);
  APL1[1] = f(0,0)*f(1,0);
  APL1[2] = f(0,0)*f(1,0);
  APL1[3] = f(1,0)*f(1,0);
  for (unsigned int i = 0; i < 4; ++i) 
    APL1[i] /= (S[1] - s)*(S[1] - S[0]);
  
  vector<double> APL2(4);
  APL2[0] = f(0,1)*f(0,1);
  APL2[1] = f(0,1)*f(1,1);
  APL2[2] = f(0,1)*f(1,1);
  APL2[3] = f(1,1)*f(1,1);
  for (unsigned int i = 0; i < 4; ++i) 
    APL2[i] /= (S[2] - s)*(S[2] - S[0]);

  vector<double> ALL(4);
  for (unsigned int i = 0; i < 4; ++i){
    ALL[i] = APL1[i] + APL2[i] + AQL[i];
    ALL[i] *= (s - S[0])/emKc2;
  }
  double detLL = ALL[0]*ALL[3] - ALL[1]*ALL[2];
  
  complex<double> denominator = 1. + zf1 * ALL[0] + zf2 * ALL[3] + zf1 * zf2 * detLL;

  vector<complex<double> > T(4);
  for (unsigned int i = 0; i < 4; ++i)
    T[i] = ALL[i];
  T[0] += zf2 * detLL;
  T[3] += zf1 * detLL;
  
  for (unsigned int i = 0; i < 4; ++i){
    T[i] /= denominator;
    T[i] *= emKc2/(s - S[0]);
  }
  
  complex<double> alpha(alpha1, alpha2);
  //complex<double> alpha(0.09, -0.5);
  const complex<double> amp = T[0] + alpha * T[1];
    
  //if (sqrt(s) < 1.6) 
  return amp;
  //else return zi;
} 

complex<double>
Flatte(double s, double m1, double m2, double m0, double g1, double g2)
{
  double m = sqrt(s);

  const complex<double> imag(0, 1);
  
  double rho_pi = 2. * breakupMomentum(s, m1, m2) / m;
  complex<double> rho_K   = 2. * breakupMomentumComplex(s, mK, mK) / m;

  double numerator = m0 * g1*rho_pi; 

  complex<double> denominator = m0*m0 - s - imag * m0 * ( g1*rho_pi + g2*rho_K);

  return numerator / denominator;
}

complex<double>
BW(double s, double m1, double m2,
   int J, double m0, double Gamma0)
{
  double m = sqrt(s);
  double q0 = breakupMomentum(m0*m0, m1, m2);
  double q = breakupMomentum(s, m1, m2);

  double Gamma = Gamma0 * m0/m * q/q0 * pow(blattWeisskopf(J, q) / blattWeisskopf(J, q0), 2);

  double numerator = sqrt(Gamma);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*Gamma);

  double norm = 1/(m0*sqrt(Gamma0));

  return numerator / denominator / norm;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b,
	  double m0, double Gamma0, int Ja, int Jb, double BR)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(Ja, qa) / blattWeisskopf(Ja, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(Jb, qb) / blattWeisskopf(Jb, q0b), 2);

  double numerator = sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*(BR*Gammaa + (1-BR)*Gammab));

  return numerator / denominator;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b, double m1c, double m2c,
	  double m0, double Gamma0, int Ja, int Jb, int Jc, double BRb, double BRc)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(Ja, qa) / blattWeisskopf(Ja, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(Jb, qb) / blattWeisskopf(Jb, q0b), 2);

  double q0c = breakupMomentum(m0*m0, m1c, m2c);
  double qc = breakupMomentum(s, m1c, m2c);

  double Gammac = Gamma0 * m0/m * qc/q0c * pow(blattWeisskopf(Jc, qc) / blattWeisskopf(Jc, q0c), 2);

  double numerator = sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*((1 - BRb - BRc)*Gammaa + BRb*Gammab + BRc*Gammac));

  double norm = 1/(m0*sqrt(Gamma0));

  return numerator / denominator / norm;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b, double m1c, double m2c, double m1d, double m2d,
	  double m0, double Gamma0, int Ja, int Jb, int Jc, int Jd, double BRb, double BRc, double BRd)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(Ja, qa) / blattWeisskopf(Ja, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(Jb, qb) / blattWeisskopf(Jb, q0b), 2);

  double q0c = breakupMomentum(m0*m0, m1c, m2c);
  double qc = breakupMomentum(s, m1c, m2c);

  double Gammac = Gamma0 * m0/m * qc/q0c * pow(blattWeisskopf(Jc, qc) / blattWeisskopf(Jc, q0c), 2);

  double q0d = breakupMomentum(m0*m0, m1d, m2d);
  double qd = breakupMomentum(s, m1d, m2d);

  double Gammad = Gamma0 * m0/m * qd/q0d * pow(blattWeisskopf(Jd, qd) / blattWeisskopf(Jd, q0d), 2);

  double numerator = sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*((1 - BRb - BRc -BRd)*Gammaa + BRb*Gammab + BRc*Gammac + BRd*Gammad));

  double norm = 1/(m0*sqrt(Gamma0));

  return numerator / denominator / norm;
}

complex<double>
BW_a2_pietap(double s)
{
  return BW(s, mPi, mEtaP, 2, 1.3183, 0.107);
}

complex<double>
BW_a2_pieta(double s)
{
  return BW(s, mPi, mEta, 2, 1.3183, 0.107);
}

complex<double>
BW_a2_pieta_coupled(double s)
{
  return BWcoupled(s, mPi, mEta, mPi, 0.77, 1.3183, 0.107, 2, 1, 0.2);
}

complex<double>
BW_a2_pietap_coupled(double s)
{
  return BWcoupled(s, mPi, mEtaP, mPi, 0.77, mPi, mEta, 1.3183, 0.107, 2, 1, 2, 0.70/0.85, 0.145/0.85);
  //return BWcoupled(s, mPi, mEtaP, mPi, 0.77, 2, 1.3183, 0.107, 0.99) / abs(BWcoupled(1.3183*1.3183, mPi, mEtaP, mPi, 0.77, 2, 1.3183, 0.107, 0.99));
}

complex<double>
BW_a2_pirho(double s)
{
  return BW(s, mPi, .77, 1, 1.3183, 0.107);
}

complex<double>
BW_a4_pieta(double s)
{
  return BW(s, mPi, mEta, 4, 2.001, 0.235);
}

complex<double>
BW_a4_pietap(double s)
{
  return BW(s, mPi, mEtaP, 4, 2.001, 0.235) / abs(BW(2.001*2.001, mPi, mEtaP, 4, 2.001, 0.235));
}

complex<double>
BW_exotic_pietap(double s)
{
  return BW(s, mPi, mEtaP, 1, 1.579, 0.340);
}

complex<double>
BW_exotic_pieta(double s)
{
  return BW(s, mPi, mEta, 1, 1.597, 0.340);
}

#if 0
complex<double>
Kmatrix_2_to_1(double s, double m1, double m2, int J,
	       double mA, double GammaA0, double mB, double GammaB0)
{
  // Gamma1 and Gamma2 also contain the relative coupling strengths.
  double m = sqrt(s);
  double q = breakupMomentum(s, m1, m2);
  double q0A = breakupMomentum(mA*mA, m1, m2);
  double q0B = breakupMomentum(mB*mB, m1, m2);

  double GammaA = GammaA0 * mA/m * q/q0A * pow(blattWeisskopf(J, q) / blattWeisskopf(J, q0A), 2);
  double GammaB = GammaB0 * mB/m * q/q0B * pow(blattWeisskopf(J, q) / blattWeisskopf(J, q0B), 2);


  double numerator = (mA * GammaA + mB * GammaB);
  double realDenom = (mA*mA - s)*(mB*mB - s);
  double imagDenom = -mA*GammaA*(mB*mB - s) - mB*GammaB*(mA*mA - s);

  return numerator / complex<double>(realDenom, imagDenom);  
}

complex<double>
Kmatrix_2_to_2(size_t iChannel,
	       double s, double m11, double m12, double m21, double m22,
	       int J,
	       double mA, double GammaA0, complex<double> phase, double mB, double GammaB0)
{
  double m = sqrt(s);
  double q12 = breakupMomentum2(s, m11, m12);
  double q01A = breakupMomentum(mA*mA, m11, m12);
  double q01B = breakupMomentum(mB*mB, m11, m12);
  double q22 = breakupMomentum2(s, m21, m22);
  double q02A = breakupMomentum(mA*mA, m21, m22);
  double q02B = breakupMomentum(mB*mB, m21, m22);

  double rho1real = 2*breakupMomentum(s, m11, m12) / m;
  double rho2real = 2*breakupMomentum(s, m21, m22) / m;
  double rhoreal[2] = { rho1real, rho2real };

  complex<double> rho1;
  if (q12 > 0)
    rho1 = 2.*sqrt(q12)/m;
  else
    rho1 = complex<double>(0,1)*2.*sqrt(-q12)/m;
  complex<double> rho2;
  if (q22 > 0)
    rho2 = 2.*sqrt(q22)/m;
  else
    rho2 = complex<double>(0,1)*2.*sqrt(-q22)/m;

  //complex<double> rho[2] = { rho1, rho2 };

  double g1 = 0.15;
  double g2 = 0.85;

  double Gamma1A = g1*GammaA0 * blattWeisskopf2(J,q12) / pow(blattWeisskopf(J,q01A),2);
  double Gamma1B = g1*GammaB0 * blattWeisskopf2(J,q12) / pow(blattWeisskopf(J,q01B),2);

  double Gamma2A = g2*GammaA0 * blattWeisskopf2(J,q22) / pow(blattWeisskopf(J,q02A),2);
  double Gamma2B = g2*GammaB0 * blattWeisskopf2(J,q22) / pow(blattWeisskopf(J,q02B),2);

  // The K matrix is multiplied with the product of its poles to avoid
  // singularities.
  double K[4];
  K[0] = (mA*sqrt(Gamma1A*Gamma1A) * (mB*mB - s)
	  + mB*sqrt(Gamma1B*Gamma1B) * (mA*mA - s));
  K[1] = (mA*sqrt(Gamma1A*Gamma2A) * (mB*mB - s)
	  + mB*sqrt(Gamma1B*Gamma2B) * (mA*mA - s));
  K[2] = (mA*sqrt(Gamma2A*Gamma1A) * (mB*mB - s)
	  + mB*sqrt(Gamma2B*Gamma1B) * (mA*mA - s));
  K[3] = (mA*sqrt(Gamma2A*Gamma2A) * (mB*mB - s)
	  + mB*sqrt(Gamma2B*Gamma2B) * (mA*mA - s));

  double poles = (mA*mA - s)*(mB*mB - s);

  // 1 - i K rho multiplied by poles
  complex<double> KK[4];
  const complex<double> ii(0,1);
  KK[0] = poles - ii*K[0]*rho1;
  KK[1] =       - ii*K[1]*rho2;
  KK[2] =       - ii*K[2]*rho1;
  KK[3] = poles - ii*K[3]*rho2;

  complex<double> invKK[4];
  if (!invertComplex2x2Matrix(KK, invKK))
    {
      cerr << "couldn't invert at m = " << m << endl;
      return 0;
    }

  // Production vector, singularities multiplied away
  complex<double> P[2];
  P[0] = (mA*sqrt(Gamma1A) * (mB*mB - s)
	  + phase*mB*sqrt(Gamma1B) * (mA*mA - s));
  P[1] = (mA*sqrt(Gamma2A) * (mB*mB - s)
	  +  phase*mB*sqrt(Gamma2B) * (mA*mA - s));
  complex<double> T[2];
  for (int i = 0; i < 2; i++)
    {
      T[i] = 0;
      for (int j = 0; j < 2; j++)
	T[i] += invKK[2*i+j]*P[j]*rhoreal[i];
    }
  return T[iChannel];
}

complex<double>
Kmatrix_2_to_3(size_t iChannel, double s,
	       double m11, double m12, double m21, double m22, double m31, double m32,
	       int J,
	       double mA, double GammaA0, complex<double> phase, double mB, double GammaB0)
{
  double m = sqrt(s);
  double q12 = breakupMomentum2(s, m11, m12);
  double q01A = breakupMomentum(mA*mA, m11, m12);
  double q01B = breakupMomentum(mB*mB, m11, m12);
  double q22 = breakupMomentum2(s, m21, m22);
  double q02A = breakupMomentum(mA*mA, m21, m22);
  double q02B = breakupMomentum(mB*mB, m21, m22);
  double q32 = breakupMomentum2(s, m31, m32);
  double q03A = breakupMomentum(mA*mA, m31, m32);
  double q03B = breakupMomentum(mB*mB, m31, m32);

  complex<double> rho1;
  if (q12 > 0)
    rho1 = 2.*sqrt(q12)/m;
  else
    rho1 = complex<double>(0,1)*2.*sqrt(-q12)/m;
  complex<double> rho2;
  if (q22 > 0)
    rho2 = 2.*sqrt(q22)/m;
  else
    rho2 = complex<double>(0,1)*2.*sqrt(-q22)/m;
  complex<double> rho3;
  if (q32 > 0)
    rho3 = 2.*sqrt(q32)/m;
  else
    rho3 = complex<double>(0,1)*2.*sqrt(-q32)/m;

  double rho1real = 2*breakupMomentum(s, m11, m12) / m;
  double rho2real = 2*breakupMomentum(s, m21, m22) / m;
  double rho3real = 2*breakupMomentum(s, m31, m32) / m;
  double rhoreal[3] = { rho1real, rho2real, rho3real };

  double g1 = 0.4;
  double g2 = 0.4;
  double g3 = 0.2;

  double Gamma1A = g1*GammaA0 * blattWeisskopf2(J,q12) / pow(blattWeisskopf(J,q01A),2);
  double Gamma1B = 0*g1*GammaB0 * blattWeisskopf2(J,q12) / pow(blattWeisskopf(J,q01B),2);

  double Gamma2A = g2*GammaA0 * blattWeisskopf2(J,q22) / pow(blattWeisskopf(J,q02A),2);
  double Gamma2B = 0*g2*GammaB0 * blattWeisskopf2(J,q22) / pow(blattWeisskopf(J,q02B),2);

  double Gamma3A = g3*GammaA0 * blattWeisskopf2(J,q32) / pow(blattWeisskopf(J,q03A),2);
  double Gamma3B = 0*g3*GammaB0 * blattWeisskopf2(J,q32) / pow(blattWeisskopf(J,q03B),2);

#if 1
  // The K matrix is multiplied with the product of its poles to avoid
  // singularities.
  complex<double> K[9];
  K[0] = (mA*sqrt(Gamma1A*Gamma1A) * (mB*mB - s)
	  + mB*sqrt(Gamma1B*Gamma1B) * (mA*mA - s));
  K[1] = (mA*sqrt(Gamma1A*Gamma2A) * (mB*mB - s)
	  + mB*sqrt(Gamma1B*Gamma2B) * (mA*mA - s));
  K[2] = (mA*sqrt(Gamma1A*Gamma3A) * (mB*mB - s)
	  + mB*sqrt(Gamma1B*Gamma3B) * (mA*mA - s));
  K[3] = (mA*sqrt(Gamma2A*Gamma1A) * (mB*mB - s)
	  + mB*sqrt(Gamma2B*Gamma1B) * (mA*mA - s));
  K[4] = (mA*sqrt(Gamma2A*Gamma2A) * (mB*mB - s)
	  + mB*sqrt(Gamma2B*Gamma2B) * (mA*mA - s));
  K[5] = (mA*sqrt(Gamma2A*Gamma3A) * (mB*mB - s)
	  + mB*sqrt(Gamma2B*Gamma3B) * (mA*mA - s));
  K[6] = (mA*sqrt(Gamma3A*Gamma1A) * (mB*mB - s)
	  + mB*sqrt(Gamma3B*Gamma1B) * (mA*mA - s));
  K[7] = (mA*sqrt(Gamma3A*Gamma2A) * (mB*mB - s)
	  + mB*sqrt(Gamma3B*Gamma2B) * (mA*mA - s));
  K[8] = (mA*sqrt(Gamma3A*Gamma3A) * (mB*mB - s)
	  + mB*sqrt(Gamma3B*Gamma3B) * (mA*mA - s));

  double poles = (mA*mA - s)*(mB*mB - s);

  // 1 - i K rho multiplied by poles
  complex<double> KK[9];
  const complex<double> ii(0,1);
  KK[0] = poles - ii*K[0]*rho1;
  KK[1] =       - ii*K[1]*rho2;
  KK[2] =       - ii*K[2]*rho3;
  KK[3] =       - ii*K[3]*rho1;
  KK[4] = poles - ii*K[4]*rho2;
  KK[5] =       - ii*K[5]*rho3;
  KK[6] =       - ii*K[6]*rho1;
  KK[7] =       - ii*K[7]*rho2;
  KK[8] = poles - ii*K[8]*rho3;

  complex<double> invKK[9];
  if (!invertComplex3x3Matrix(KK, invKK))
    {
      cerr << "couldn't invert at m = " << m << endl;
      return 0;
    }

  // Production vector, singularities multiplied away
  complex<double> P[3];
  P[0] = (mA*sqrt(Gamma1A) * (mB*mB - s)
	  + phase*mB*sqrt(Gamma1B) * (mA*mA - s));
  P[1] = (mA*sqrt(Gamma2A) * (mB*mB - s)
	  + phase*mB*sqrt(Gamma2B) * (mA*mA - s));
  P[3] = (mA*sqrt(Gamma3A) * (mB*mB - s)
	  + phase*mB*sqrt(Gamma3B) * (mA*mA - s));
#endif

#if 0
  complex<double> K[9];
  K[0] = (mA*sqrt(Gamma1A*Gamma1A / rho1 / rho1) / (mA*mA - s)
	  + mB*sqrt(Gamma1B*Gamma1B / rho1 / rho1) / (mB*mB - s));
  K[1] = (mA*sqrt(Gamma1A*Gamma2A / rho1 / rho2) / (mA*mA - s)
	  + mB*sqrt(Gamma1B*Gamma2B / rho1 / rho2) / (mB*mB - s));
  K[2] = (mA*sqrt(Gamma1A*Gamma3A / rho1 / rho3) / (mA*mA - s)
	  + mB*sqrt(Gamma1B*Gamma3B / rho1 / rho3) / (mB*mB - s));
  K[3] = (mA*sqrt(Gamma2A*Gamma1A / rho2 / rho1) / (mA*mA - s)
	  + mB*sqrt(Gamma2B*Gamma1B / rho2 / rho1) / (mB*mB - s));
  K[4] = (mA*sqrt(Gamma2A*Gamma2A / rho2 / rho2) / (mA*mA - s)
	  + mB*sqrt(Gamma2B*Gamma2B / rho2 / rho2) / (mB*mB - s));
  K[5] = (mA*sqrt(Gamma2A*Gamma3A / rho2 / rho3) / (mA*mA - s)
	  + mB*sqrt(Gamma2B*Gamma3B / rho2 / rho3) / (mB*mB - s));
  K[6] = (mA*sqrt(Gamma3A*Gamma1A / rho3 / rho1) / (mA*mA - s)
	  + mB*sqrt(Gamma3B*Gamma1B / rho3 / rho1) / (mB*mB - s));
  K[7] = (mA*sqrt(Gamma3A*Gamma2A / rho3 / rho2) / (mA*mA - s)
	  + mB*sqrt(Gamma3B*Gamma2B / rho3 / rho2) / (mB*mB - s));
  K[8] = (mA*sqrt(Gamma3A*Gamma3A / rho3 / rho3) / (mA*mA - s)
	  + mB*sqrt(Gamma3B*Gamma3B / rho3 / rho3) / (mB*mB - s));

  double poles = 1; //(mA*mA - s)*(mB*mB - s);

  // 1 - i K rho multiplied by poles
  complex<double> KK[9];
  const complex<double> ii(0,1);
  KK[0] = complex<double>(poles) - ii*K[0]*rho1;
  KK[1] =                        - ii*K[1]*rho2;
  KK[2] =                        - ii*K[2]*rho3;
  KK[3] =                        - ii*K[3]*rho1;
  KK[4] = complex<double>(poles) - ii*K[4]*rho2;
  KK[5] =                        - ii*K[5]*rho3;
  KK[6] =                        - ii*K[6]*rho1;
  KK[7] =                        - ii*K[7]*rho2;
  KK[8] = complex<double>(poles) - ii*K[8]*rho3;

  complex<double> invKK[9];
  if (!invertComplex3x3Matrix(KK, invKK))
    {
      cerr << "couldn't invert at m = " << m << endl;
      return 0;
    }

  // Production vector, singularities multiplied away
  complex<double> P[3];
  P[0] = (mA*sqrt(Gamma1A / rho1) / (mA*mA - s)
	  + phase*mB*sqrt(Gamma1B / rho1) / (mB*mB - s));
  P[1] = (mA*sqrt(Gamma2A / rho2) / (mA*mA - s)
	  + phase*mB*sqrt(Gamma2B / rho2) / (mB*mB - s));
  P[3] = (mA*sqrt(Gamma3A / rho3) / (mA*mA - s)
	  + phase*mB*sqrt(Gamma3B / rho3) / (mB*mB - s));
#endif

  /*
  complex<double> T[4];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      {
	T[2*i+j] = 0;
	for (int k = 0; k < 2; k++)
	  T[2*i+j] += invKK[2*i+k]*K[2*k+j];
	T[2*i+j] *= sqrt(rho[i]*rho[j]);
      }
  return T[iChannel * 2 + 1];
  */
  complex<double> T[3];
  for (int i = 0; i < 3; i++)
    {
      T[i] = 0;
      for (int j = 0; j < 3; j++)
	T[i] += invKK[3*i+j]*P[j];
      T[i] *= rhoreal[i];
    }
  return T[iChannel];
}


/*
double
fitFunc(double *x, double *p)
{
  double s = x[0]*x[0];

  double iK892 = p[0];
  double mK892 = p[1];
  double wK892 = p[2];
  double iK1430 = p[3];
  double mK1430 = p[4];
  double wK1430 = p[5];
  double p0 = p[6];
  double p1 = p[7];
  double p2 = p[8];
  double p3 = p[9];
  double p4 = p[10];

  complex<double> bwK892 = BW(s, mPi, mKs, 1, mK892, wK892);
  double intK892 = norm(bwK892)*breakupMomentum(s, mPi, mKs);

  complex<double> bwK1430 = BW(s, mPi, mKs, 2, mK1430, wK1430);
  double intK1430 = norm(bwK1430)*breakupMomentum(s, mPi, mKs);

  double poly = p0*(1 + p1*x[0] + p2*pow(x[0],2) + p3*pow(x[0],3) + p4*pow(x[0],4));
  return iK892*intK892 + iK1430*intK1430 + poly;
}


void bw()
{
  hmKpifit->Draw();
  f = new TF1("f", fitFunc, 0.65, 1.65, 11);
  f->FixParameter(1, .892);
  f->FixParameter(2, .030);
  f->FixParameter(4, 1.43);
  f->FixParameter(5, .1);
  f->FixParameter(6,0);
  f->FixParameter(7,0);
  f->FixParameter(8,0);
  f->FixParameter(9,0);
  f->FixParameter(10,0);
  hmKpifit->Fit("f");
}
*/


void
bw()
{
  size_t nPoints = 20000;

  for (int iStrength = 0; iStrength < 5; iStrength++)
    {
      double strength = iStrength * .4;
      for (int iPlot = 0; iPlot < 5; iPlot++)
	{
	  complex<double> phase(strength*cos(iPlot*M_PI/4), strength*sin(iPlot*M_PI/4));
	  //cout << phase << endl;

	  double x[nPoints];
	  double y1[nPoints];
	  double y2[nPoints];
	  double y3[nPoints];
	  double y4[nPoints];
	  double y5[nPoints];
	  double y6[nPoints];
	  double y7[nPoints];
	  double y8[nPoints];

	  double m1 = mEta;
	  double m2 = mPi;

	  for (size_t i = 1; i <= nPoints; i++)
	    {
	      double m = m1 + m2 + 2./nPoints*i;
	      complex<double> bw = (//Kmatrix_2_to_2(0, m*m,  mEta, mPi, mEtaP, mPi,
				    //		   2,
				    //		   1.32, .107, phase, 1.745, 0.3)
				    //Kmatrix_2_to_2(1, m*m,  mEta, mPi, mRho, mPi,
				    //		   2,
				    //		   1.3201, .107, phase, 1.745, 0.3)
				    Kmatrix_2_to_3(0, m*m,
				    		   mEta, mPi, mEtaP, mPi, 0.77, mPi,
				    		   2,
				    		   1.32, 0.107, phase, 1.745, 0.3)
				    );

	      x[i-1] = m;
	      y1[i-1] = real(bw);
	      y2[i-1] = imag(bw);
	      y3[i-1] = norm(bw);
	      y4[i-1] = arg(bw);

	      bw = (//Kmatrix_2_to_1(m*m, m1, m2, 2,
		    //		  1.32, .13, 1.8, .8)
		    Kmatrix_2_to_2(0, m*m, mEta, mPi, mRho, mPi,
				   2,
				   1.3201, 0.107, phase, 1.745, 0.3)
		    );

	      y5[i-1] = real(bw);
	      y6[i-1] = imag(bw);
	      y7[i-1] = norm(bw);
	      y8[i-1] = arg(bw);
	    }

	  TCanvas* c = new TCanvas();
	  c->Divide(2,3);

	  TGraph* g1 = new TGraph(nPoints, x, y3);
	  TGraph* g2 = new TGraph(nPoints, x, y4);
	  TGraph* g3 = new TGraph(nPoints, x, y1);
	  TGraph* g4 = new TGraph(nPoints, x, y2);
	  TGraph* gArgand1 = new TGraph(nPoints, y1, y2);
	  TGraph* g5 = new TGraph(nPoints, x, y7); g5->SetLineColor(kRed);
	  TGraph* g6 = new TGraph(nPoints, x, y8); g6->SetLineColor(kRed);
	  TGraph* g7 = new TGraph(nPoints, x, y5); g7->SetLineColor(kRed);
	  TGraph* g8 = new TGraph(nPoints, x, y6); g8->SetLineColor(kRed);
	  TGraph* gArgand2 = new TGraph(nPoints, y5, y6); gArgand2->SetLineColor(kRed);
	  c->cd(1);
	  g1->Draw("AL"); g5->Draw("L");
	  c->cd(2);
	  g2->Draw("AL"); g6->Draw("L");
	  c->cd(3);
	  g3->Draw("AL"); g7->Draw("L");
	  c->cd(4);
	  g4->Draw("AL"); g8->Draw("L");
	  c->cd(5);
	  gArgand1->Draw("AL");
	  c->cd(6);
	  gArgand2->Draw("AL");
	}
    }
}


#if 0
void
bw()
{
  ofstream of("bw.txt");

  double x[2000];
  double y1[2000];
  double y2[2000];
  double y3[2000];
  double y4[2000];
  double y5[2000];

  for (int i = 1; i <= 2000; i++)
    {
      double m = mPi + mEtaP + i*0.0005; //mPi + 0.77 + i*0.001;
      complex<double> z =BW_a2_pietap_coupled(m*m);
      complex<double> z2 = BW(m*m, mEtaP, mPi, 2, 1.3183, 0.107);
      of << m << " "
	   << real(z) << " " << imag(z) << " " << norm(z) << " " << arg(z)*180/3.142 << " "
	 << real(z2) << " " << imag(z2) << " " << norm(z2) << " " << arg(z2)*180/3.142 << " "
	 << arg(z / z2) << " " << breakupMomentum(m*m, mEtaP, mPi)
	 << " " << blattWeisskopf(4, breakupMomentum(m*m, mEta, mPi))
	 << endl;

      x[i-1] = m;
      y1[i-1] = norm(z);
      y2[i-1] = norm(z2);
      y3[i-1] = arg(z)*180/3.141592653589793;
      y4[i-1] = arg(z2)*180/3.141592653589793;
      y5[i-1] = y3[i-1] - y4[i-1];
    }

  TGraph* g1 = new TGraph(2000, x, y1);
  TGraph* g2 = new TGraph(2000, x, y2);
  g1->SetLineStyle(2);

  g1->SetTitle("a_{2} intensity with different BWs");
  g1->Draw("AL"); // coupled -> further right
  g2->Draw("L");

  new TCanvas;

  TGraph* g3 = new TGraph(2000, x, y3);
  TGraph* g4 = new TGraph(2000, x, y4);
  TGraph* g5 = new TGraph(2000, x, y5);
  g3->SetLineStyle(2);
  g5->SetLineStyle(3);

  g3->SetTitle("a_{2} phase with different BWs");
  g3->Draw("AL"); // coupled -> further right
  g4->Draw("L");
  g5->Draw("L");
}
#endif

#endif
