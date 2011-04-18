#include <complex>
#include <iostream>

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

  double numerator = blattWeisskopf(J, q);
  complex<double> denominator = complex<double>(m0*m0 - s, -m0*Gamma);

  return numerator / denominator;
}

complex<double>
BW_a2_pietap(double s)
{
  return BW(s, mPi, mEtaP, 2, 1.3183, 0.107);
}

complex<double>
BW_exotic_pietap(double s)
{
  return BW(s, mPi, mEtaP, 1, 1.579, 0.340);
}

complex<double>
BW_a2_pieta(double s)
{
  return BW(s, mPi, mEta, 2, 1.3183, 0.107);
}

complex<double>
BW_exotic_pieta(double s)
{
  return BW(s, mPi, mEta, 1, 1.597, 0.340);
}

#if 0
void
bw()
{
  for (int i = 0; i < 2000; i++)
    {
      double m = mPi + mEtaP + i*0.001;
      complex<double> z = BW_a2_pietap(m*m);
      complex<double> z2 = BW_exotic_pietap(m*m);
      cout << m << " "
	   << real(z) << " " << imag(z) << " " << norm(z)*breakupMomentum(m*m, mPi, mEtaP) << " " << arg(z) << " "
	   << real(z2) << " " << imag(z2) << " " << norm(z2) << " " << arg(z2) << " "
	   << arg(z / z2)
	   << endl;
    }
}
#endif
