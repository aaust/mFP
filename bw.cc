#include <complex>
#include <iostream>

using namespace std;

double
blattWeisskopf(int L, double p)
{
  double z = p / 0.1973;  // 1fm interaction radius
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
  return BW(s, .139, .958, 2, 1.3183, 0.107);
}

#if 0
void
bw()
{
  for (int i = 0; i < 2000; i++)
    {
      double m = .139 + .958 + i*0.001;
      complex<double> z = BW_a2_pietap(m*m);
      cout << m << " "
	   << real(z) << " " << imag(z) << " " << abs(z) << " " << arg(z)
	   << endl;
    }
}
#endif
