#ifndef FITMODEL_H__
#define FITMODEL_H__

#include <complex>
#include <vector>

class fitModel {
 public:
  fitModel() {};
  virtual void evaluateAt(double mass_, const std::vector<double>& x) = 0;
  virtual std::complex<double> valueForWave(const char *name) const = 0;

  static fitModel* getFitModelForName(const string& name);
};

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

#endif
