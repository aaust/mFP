#include <complex>
#include <iostream>
#include <fstream>

using namespace std;

const double mPi = .13957018;
const double mPi0 = .1349766;
const double mEta = .547853;
const double mEtaP = .95778;

const double hbarc = 0.1973269631; // MeV fm

double
blattWeisskopf(int L, double p)
{
  double z = pow(p / hbarc, 2);  // 1fm interaction radius
  double result;
  switch (L) {
  case 0:
    result = 1; break;
  case 1:
    result = sqrt(2*z/(z+1)); break;
  case 2:
    result = sqrt(13*z*z/(z*pow(z-3, 2) + 9*z)); break;
  case 3:
    result = sqrt(277*z*z*z/(z*pow(z-15, 2)+ 9*pow(2*z-5, 2))); break;
  case 4:
    result = sqrt(12746*pow(z,4)/(pow(z*z-45*z+105,2) + 25*z*pow(2*z-21,2))); break;
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

complex<double>
BW(double s, double m1, double m2,
   int J, double m0, double Gamma0)
{
  double m = sqrt(s);
  double q0 = breakupMomentum(m0*m0, m1, m2);
  double q = breakupMomentum(s, m1, m2);

  double Gamma = Gamma0 * m0/m * q/q0 * pow(blattWeisskopf(J, q) / blattWeisskopf(J, q0), 2);

  double numerator = blattWeisskopf(J, q) * sqrt(Gamma);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*Gamma);

  return numerator / denominator;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b,
	  int J, double m0, double Gamma0, double BR)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(J, qa) / blattWeisskopf(J, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(J, qb) / blattWeisskopf(J, q0b), 2);

  double numerator = blattWeisskopf(J, qa) * sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*(BR*Gammaa + (1-BR)*Gammab));

  return numerator / denominator;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b, double m1c, double m2c,
	  int J, double m0, double Gamma0, double BRb, double BRc)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(J, qa) / blattWeisskopf(J, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(J, qb) / blattWeisskopf(J, q0b), 2);

  double q0c = breakupMomentum(m0*m0, m1c, m2c);
  double qc = breakupMomentum(s, m1c, m2c);

  double Gammac = Gamma0 * m0/m * qc/q0c * pow(blattWeisskopf(J, qc) / blattWeisskopf(J, q0c), 2);

  double numerator = blattWeisskopf(J, qa) * sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*((1 - BRb - BRc)*Gammaa + BRb*Gammab + BRc*Gammac));

  return numerator / denominator;
}

complex<double>
BWcoupled(double s, double m1a, double m2a, double m1b, double m2b, double m1c, double m2c, double m1d, double m2d,
	  int J, double m0, double Gamma0, double BRb, double BRc, double BRd)
{
  double m = sqrt(s);
  double q0a = breakupMomentum(m0*m0, m1a, m2a);
  double qa = breakupMomentum(s, m1a, m2a);

  double Gammaa = Gamma0 * m0/m * qa/q0a * pow(blattWeisskopf(J, qa) / blattWeisskopf(J, q0a), 2);

  double q0b = breakupMomentum(m0*m0, m1b, m2b);
  double qb = breakupMomentum(s, m1b, m2b);

  double Gammab = Gamma0 * m0/m * qb/q0b * pow(blattWeisskopf(J, qb) / blattWeisskopf(J, q0b), 2);

  double q0c = breakupMomentum(m0*m0, m1c, m2c);
  double qc = breakupMomentum(s, m1c, m2c);

  double Gammac = Gamma0 * m0/m * qc/q0c * pow(blattWeisskopf(J, qc) / blattWeisskopf(J, q0c), 2);

  double q0d = breakupMomentum(m0*m0, m1d, m2d);
  double qd = breakupMomentum(s, m1d, m2d);

  double Gammad = Gamma0 * m0/m * qd/q0d * pow(blattWeisskopf(J, qd) / blattWeisskopf(J, q0d), 2);

  double numerator = blattWeisskopf(J, qa) * sqrt(Gammaa);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*((1 - BRb - BRc -BRd)*Gammaa + BRb*Gammab + BRc*Gammac + BRd*Gammad));

  return numerator / denominator;
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
  return BWcoupled(s, mPi, mEta, mPi, 0.77, 2, 1.3183, 0.107, 0.2);
}

complex<double>
BW_a2_pietap_coupled(double s)
{
  return BWcoupled(s, mPi, mEtaP, mPi, 0.77, mPi, mEta, 2, 1.3183, 0.107, 0.70/0.85, 0.15/0.85);
}

complex<double>
BW_a2_pirho(double s)
{
  return BW(s, mPi, .77, 2, 1.3183, 0.107);
}

complex<double>
BW_a4_pieta(double s)
{
  return BW(s, mPi, mEta, 4, 2.001, 0.235);
}

complex<double>
BW_a4_pietap(double s)
{
  return BW(s, mPi, mEtaP, 4, 2.001, 0.235);
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

#if 1
void
bw()
{
  ofstream of("bw.txt");

  for (int i = 1; i < 2000; i++)
    {
      double m = mPi + mEtaP + i*0.001; //mPi + 0.77 + i*0.001;
      complex<double> z = BW_a4_pietap(m*m);
      complex<double> z2 = BW_a4_pieta(m*m);
      of << m << " "
	   << real(z) << " " << imag(z) << " " << norm(z)*breakupMomentum(m*m, mEtaP, mPi) << " " << arg(z)*180/3.142 << " "
	 << real(z2) << " " << imag(z2) << " " << norm(z2)*breakupMomentum(m*m, mEta, mPi) << " " << arg(z2)*180/3.142 << " "
	 << arg(z / z2) << " " << breakupMomentum(m*m, mEta, mPi)
	 << " " << breakupMomentum(m*m, mEtaP, mPi)
	   << endl;
    }
}
#endif
