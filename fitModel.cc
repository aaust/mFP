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
  fitModelEtaPi() : mass(0), Dwave(0), DwaveBG(0), Pwave(0), Gwave(0), GwaveBG(0) {};
  void evaluateAt(double mass_, const std::vector<double>& x);
  std::complex<double> valueForWave(const char* name) const;
private:
  double mass;
  std::complex<double> phaseD;
  std::complex<double> Dwave;
  std::complex<double> DwaveBG;
  std::complex<double> Pwave;
  std::complex<double> Gwave;
  std::complex<double> GwaveBG;
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
  DwaveBG = x[0]*pow(phaseSpace, 2+.5)/mass*exp(-x[19]*phaseSpace*phaseSpace)*x[20]*(1 + mass*x[21] + mass*mass*x[22]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0 && mass > 0.77 + mPi)
    Dwave = (x[0]*sqrt(mass)*BWcoupled(mass*mass, mPi, mEta, mPi, 0.77, par[0], par[1], 2, 2, 0.15/0.85)
	     + DwaveBG
	     //+ (complex<double>(par[4], par[5])
	     //	*BW(mass*mass, m1, m2, 2, par[2], par[3]))
	     );
  phaseD = 1;
  if (abs(Dwave) != 0)
    phaseD = Dwave / abs(Dwave);
  Dwave /= phaseD;

  par = &x[9];
  Pwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[0] > mEta + mPi)
    Pwave = complex<double>(x[7],x[8])*sqrt(mass)*BW(mass*mass, m1, m2, 1, par[0], par[1]);
  Pwave /= phaseD;

  par = &x[13];
  Gwave = 0;
  GwaveBG = complex<double>(x[11],x[12])*pow(phaseSpace, 4+.5)/mass*exp(-x[23]*phaseSpace*phaseSpace)*(x[24] + mass*x[25] + mass*mass*x[26]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Gwave = (complex<double>(x[11],x[12])*sqrt(mass)*BW(mass*mass, m1, m2, 4, par[0], par[1])
	     + GwaveBG
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
  else if (!strcmp(name, "D+BG"))
    return DwaveBG;
  else if (!strcmp(name, "phaseD"))
    return phaseD;
  else if (!strcmp(name, "P+"))
    return Pwave;
  else if (!strcmp(name, "G+"))
    return Gwave;
  else if (!strcmp(name, "G+BG"))
    return GwaveBG;

  cerr << "unknown wave '" << name << "' requested" << endl;
  exit(1);
}

class fitModelEtaPpi : public fitModel {
public:
  fitModelEtaPpi() : mass(0), Dwave(0), DwaveBG(0), Pwave(0), Gwave(0), GwaveBG(0) {};
  void evaluateAt(double mass_, const std::vector<double>& x);
  std::complex<double> valueForWave(const char* name) const;
private:
  double mass;
  std::complex<double> Dwave;
  std::complex<double> DwaveBG;
  std::complex<double> phaseD;
  std::complex<double> Pwave;
  std::complex<double> Gwave;
  std::complex<double> GwaveBG;
};

void
fitModelEtaPpi::evaluateAt(double mass_, const vector<double>& x)
{
  mass = mass_;

  double m1 = mEtaP;
  double m2 = mPi;
  double phaseSpace = breakupMomentum(mass*mass, m1, m2);

  const double *par = &x[1];
  Dwave = 0;
  DwaveBG = x[0]*pow(phaseSpace, 2+.5)/mass*exp(-x[19]*phaseSpace*phaseSpace)*x[20]*(1 + mass*x[21] + mass*mass*x[22]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Dwave = (x[0]*sqrt(mass)*BWcoupled(mass*mass, mPi, mEtaP, mPi, mEta, mPi, 0.77, par[0], par[1], 2, 2, 2, 0.15/0.85, 0.70/0.85)
	     + DwaveBG
	     //+ (complex<double>(par[4], par[5])
	     //	*BW(mass*mass, m1, m2, 2, par[2], par[3]))
	     );
  phaseD = 1;
  if (abs(Dwave) != 0)
    phaseD = Dwave / abs(Dwave);
  Dwave /= phaseD;

  par = &x[9];
  Pwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[0] > mEtaP + mPi)
    Pwave = complex<double>(x[7],x[8])*sqrt(mass)*BW(mass*mass, m1, m2, 1, par[0], par[1]);
  Pwave /= phaseD;

  par = &x[13];
  Gwave = 0;
  GwaveBG = complex<double>(x[11],x[12])*pow(phaseSpace, 4+.5)/mass*exp(-x[23]*phaseSpace*phaseSpace)*(x[24] + mass*x[25] + mass*mass*x[26]);
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Gwave = (complex<double>(x[11],x[12])*sqrt(mass)*BW(mass*mass, m1, m2, 4, par[0], par[1])
	     + GwaveBG
	     //+ (complex<double>(par[4], par[5])
	     //	*BW(mass*mass, m1, m2, 4, par[2], par[3]))
	     );
  Gwave /= phaseD;
}


complex<double>
fitModelEtaPpi::valueForWave(const char* name) const
{
  if (!strcmp(name, "D+"))
    return Dwave;
  else if (!strcmp(name, "D+BG"))
    return DwaveBG;
  else if (!strcmp(name, "phaseD"))
    return phaseD;
  else if (!strcmp(name, "P+"))
    return Pwave;
  else if (!strcmp(name, "G+"))
    return Gwave;
  else if (!strcmp(name, "G+BG"))
    return GwaveBG;

  cerr << "unknown wave '" << name << "' requested" << endl;
  exit(1);
}

class fitModelKKbarSD : public fitModel {
public:
  fitModelKKbarSD() : mass(0), Swave(0), SwaveBG(0), Dwave(0), DwaveBG(0) {};
  void evaluateAt(double mass_, const std::vector<double>& x);
  std::complex<double> valueForWave(const char* name) const;
private:
  double mass;
  std::complex<double> Swave;
  std::complex<double> SwaveBG;
  std::complex<double> Dwave;
  std::complex<double> DwaveBG;
  std::complex<double> phaseS;

  std::complex<double> f2;
  std::complex<double> f2prime;

  std::complex<double> f01370;
  std::complex<double> f01500;
  std::complex<double> f01710;
};

void
fitModelKKbarSD::evaluateAt(double mass_, const vector<double>& x)
{
  mass = mass_;
  
  double m1 = mK;
  double m2 = mK;
  //double phaseSpace = breakupMomentum(mass*mass, m1, m2);
  
  //par = &x[13];
  
  DwaveBG = complex<double>(x[8],x[9]) /** pow(mass-2*mK, x[26])*/ * exp(-1*x[10]*mass - x[11]*mass*mass);

  f2 = complex<double>(x[0],x[1]) * BW(mass*mass, m1, m2, 2, x[2], x[3]);
  f2prime = complex<double>(x[4], x[5]) * BW(mass*mass, m1, m2, 2, x[6], x[7]);
  
  Dwave = f2 + f2prime + DwaveBG;
  
  //  SwaveBG = complex<double>(x[29],x[30]) * pow(mass - 2*mK,x[31]) + exp(-1*x[32]*mass - x[33]*mass*mass);    

  f01370 = complex<double>(x[16], x[17])*BW(mass*mass, m1, m2, 0, x[18], x[19]);
  f01500 = complex<double>(x[20], x[21])*BW(mass*mass, m1, m2, 0, x[22], x[23]);
  f01710 = complex<double>(x[24],x[25])*BW(mass*mass, m1, m2, 0, x[26], x[27]);
  
  Swave = complex<double>(x[12], x[13])*BW(mass*mass, m1, m2, 0, x[14], x[15])
    + f01500
    + f01710
    + f01370;
  //+ SwaveBG;
  
  phaseS = 1;
  if (abs(Swave) != 0)
    phaseS = Swave / abs(Swave);
  Swave /= phaseS;

  Dwave /= phaseS;
    
  //SwaveBG = complex<double>(x[21],x[22])*(x[23] + mass*x[24] + mass*mass*x[25]);
}


complex<double>
fitModelKKbarSD::valueForWave(const char* name) const
{
  if (!strcmp(name, "S0"))
    return Swave;
  else if (!strcmp(name, "S0BG"))
    return SwaveBG;
  else if (!strcmp(name, "D0"))
    return Dwave;
  else if (!strcmp(name, "D0BG"))
    return DwaveBG;
  else if (!strcmp(name, "f2"))
    return f2;
  else if (!strcmp(name, "f2prime"))
    return f2prime;
  else if (!strcmp(name, "f01370"))
    return f01370;
  else if (!strcmp(name, "f01500"))
    return f01500;
  else if (!strcmp(name, "f01710"))
    return f01710;

  cerr << "unknown wave '" << name << "' requested" << endl;
  exit(1);
}

fitModel*
fitModel::getFitModelForName(const string& name)
{
  if (!strcmp(name.c_str(), "etaPi"))
    return new fitModelEtaPi();

  if (!strcmp(name.c_str(), "etaPpi"))
    return new fitModelEtaPpi();

  if (!strcmp(name.c_str(), "KKbarSD"))
    return new fitModelKKbarSD();

  cerr << "No fit model with modelName == '" << name << "' known" << endl;
  abort();
}
