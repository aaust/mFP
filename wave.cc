#include "wave.h"
#include "event.h"


complex<double>
coherent_waves::sum(const double* x, const event& e) const
{
  complex<double> result = 0;
  size_t idx = 0;
  vector<wave>::const_iterator it;
  for (it = this->waves.begin(); it != this->waves.end(); it++)
    {
      complex<double> a(x[idx], x[idx+1]);
      result += a * e.decayAmplitude(this->reflectivity,*it);
      idx += 2;
    }

  return result;
}
