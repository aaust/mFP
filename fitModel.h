#ifndef FITMODEL_H__
#define FITMODEL_H__

#include <complex>
#include <vector>

class fitModel {
 public:
  fitModel() {};
  virtual ~fitModel() {};
  virtual void evaluateAt(double mass_, const std::vector<double>& x) = 0;
  virtual std::complex<double> valueForWave(const char *name) const = 0;

  static fitModel* getFitModelForName(const string& name);
};

#endif
