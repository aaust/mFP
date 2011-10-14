#ifndef BW_H__
#define BW_H__

#include <complex>

const double mPi = .13957018;
const double mPi0 = .1349766;
const double mK = .493677;
const double mKs = .497614;
const double mEta = .547853;
const double mEtaP = .95778;

const double hbarc = 0.1973269631; // GeV fm

double blattWeisskopf(int L, double p);
double breakupMomentum(double s, double m1, double m2);
std::complex<double> BW(double s, double m1, double m2, int J, double m0, double Gamma0);
std::complex<double> BW_a2_pietap_coupled(double s);
std::complex<double> BW_a2_pieta_coupled(double s);
std::complex<double> BW_a4_pietap(double s);
std::complex<double> BWcoupled(double s, double m1a, double m2a, double m1b, double m2b,
			       int J, double m0, double Gamma0, double BR);

#endif
