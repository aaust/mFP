#include <math.h>

#include "Math/SpecFunc.h"

#include "event.h"

// Returns the real / imaginary part of the amplitude, depending on which one
// is non-zero.

// The amplitude is:  Ylm(theta, phi) - epsilon (-)^m Yl-m(theta, phi)
//                  = Ylm(theta, phi) - epsilon Ylm(theta, phi)*
//                  = Ylm(theta, 0) (e^(i m phi) - epsilon e^(-i m phi))
//                  = 2 Ylm(theta, 0) {i sin, cos}(m phi)

// NOTE the phase is ignored, as different reflectivities don't interfere

double
event::decayAmplitude(int reflectivity, int l, int m) const
{
  double spherical = ROOT::Math::sph_legendre(l, m, this->theta);

  // This absorbs the factor 2 from e^i \pm e^-i = 2 {i sin, cos}
  double factor = 1; //sqrt(nrdevents);
  if (m != 0)
    factor *= sqrt(2.);    

  if (reflectivity == +1)
    return factor*spherical*sin(m*this->phi);
  else
    return factor*spherical*cos(m*this->phi);
}


double
event::MCweight(int reflectivity, int l1, int m1, int l2, int m2) const
{
  // This is real and no conjugate is employed because of the special
  // form of the two-pseudoscalar decay amplitudes.
  return (this->decayAmplitude(reflectivity,l1,m1)
	  * this->decayAmplitude(reflectivity,l2,m2));
}


// Contribution of this event to the momement H_x(LM)
// this is D^L_{M0}(theta, phi, 0) = conj(Y^L_M(theta, phi))
complex<double>
event::momentWeight(int L, int M) const
{
  double spherical = ROOT::Math::sph_legendre(L, M, this->theta);
  return spherical * complex<double>(cos(M*this->phi), -sin(M*this->phi));
}
