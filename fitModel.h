#ifndef FITMODEL_H__
#define FITMODEL_H__

#include <complex>
#include <vector>

class fitModel {
 public:
  fitModel(double mass_, const std::vector<double>& x) {};
  virtual std::complex<double> valueForWave(const char *name) const = 0;
};

class fitModelEtaPi : public fitModel {
public:
  fitModelEtaPi(double mass_, const std::vector<double>& x);
  std::complex<double> valueForWave(const char* name) const;
private:
  double mass;
  std::complex<double> Dwave;
  std::complex<double> Pwave;
  std::complex<double> Gwave;
};

#endif
