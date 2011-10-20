#include <complex>
#include <vector>
#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;

#include "bw.h"
#include "fitModel.h"

class fitModelEtaPi : public fitModel {
public:
  fitModelEtaPi() : mass(0), Dwave(0), Pwave(0), Gwave(0) {};
  void evaluateAt(double mass_, const std::vector<double>& x);
  std::complex<double> valueForWave(const char* name) const;
private:
  double mass;
  std::complex<double> Dwave;
  std::complex<double> Pwave;
  std::complex<double> Gwave;
};

void
fitModelEtaPi::evaluateAt(double mass_, const vector<double>& x)
{
  mass = mass_;

  double m1 = mEta;
  double m2 = mPi;
  double phaseSpace = breakupMomentum(mass*mass, m1, m2);

  const double *par = &x[1];
  Dwave = 0;
  double Dbg = pow(phaseSpace, 2*2+1)/mass*exp(-x[19]*phaseSpace*phaseSpace)*(x[20] + mass*x[21] + mass*mass*x[22]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0 && mass > 0.77 + mPi)
    Dwave = x[0]*(sqrt(phaseSpace)*BWcoupled(mass*mass, mPi, mEta, mPi, 0.77, par[0], par[1], 2, 2, 0.15/0.85)
		  + Dbg
		  //+ (complex<double>(par[4], par[5])
		  //	*BW(mass*mass, m1, m2, 2, par[2], par[3]))
		  );
  complex<double> phaseD = 1;
  if (abs(Dwave) != 0)
    phaseD = Dwave / abs(Dwave);
  Dwave /= phaseD;

  par = &x[9];
  Pwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[0] > mEta + mPi)
    Pwave = complex<double>(x[7],x[8])*sqrt(phaseSpace)*BW(mass*mass, m1, m2, 1, par[0], par[1]);
  Pwave /= phaseD;

  par = &x[13];
  Gwave = 0;
  double Gbg = pow(phaseSpace, 2*4+1)/mass*exp(-x[23]*phaseSpace*phaseSpace)*(x[24] + mass*x[25] + mass*mass*x[26]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Gwave = complex<double>(x[11],x[12])*sqrt(phaseSpace)*(BW(mass*mass, m1, m2, 4, par[0], par[1])
							   + Gbg
							   //+ (complex<double>(par[4], par[5])
							   //	*BW(mass*mass, m1, m2, 4, par[2], par[3]))
							   );
  Gwave /= phaseD;
}


complex<double>
fitModelEtaPi::valueForWave(const char* name) const
{
  if (!strcmp(name, "D+"))
    return Dwave;
  else if (!strcmp(name, "P+"))
    return Pwave;
  else if (!strcmp(name, "G+"))
    return Gwave;

  cerr << "unknown wave '" << name << "' requested" << endl;
  exit(1);
}


fitModel*
fitModel::getFitModelForName(const string& name)
{
  if (!strcmp(name.c_str(), "etaPi"))
    return new fitModelEtaPi();

  cerr << "No fit model with modelName == '" << name << "' known" << endl;
  abort();
}
