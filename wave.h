#ifndef WAVE_H__
#define WAVE_H__

#include <complex>
#include <vector>
#include <map>

#include "TH1.h"
#include "TFitterMinuit.h"

using namespace std;

struct wave {
  string name;
  size_t l, m;
  size_t idx;   // Index into the fit variables.
  TH1* histIntensity;
  map<string, TH1*> mHistPhase;
  bool phaseLocked;

  wave() { histIntensity = 0; }
  wave(int ll, int mm) { l = ll; m = mm; histIntensity = 0; }
  wave(const char* name_, int ll, int mm) : name(name_), l(ll), m(mm) { histIntensity = 0; }
  wave(const char* name_, int ll, int mm, int nBins, double lower, double upper, bool phaseLocked_ = false)
    : name(name_), l(ll), m(mm), phaseLocked(phaseLocked_)
  { buildHists(nBins, lower, upper); }
  wave(const wave& o);

  ~wave() {} // histograms are owned by ROOT

  void setIndex(int idx_) { idx = idx_; }
  size_t getIndex() const { return idx; }
  // Thanks to the brilliance of Minuit2, the covariance matrix doesn't
  // contain the fixed parameters, making the mapping of fit parameters
  // to covariance matrix elements awkward.  This should help.
  size_t idxInCovariance(const TFitterMinuit* minuit) const;

  void buildHists(int nBins, double lower, double upper);
  TH1* getHistIntensity() const { return histIntensity; }
  TH1* getHistPhase(const wave& other);
  void fillHistIntensity(int iBin, const TFitterMinuit* minuit);
  void fillHistPhase(int iBin, const wave& other, const TFitterMinuit* minuit);
};

class event;

struct coherent_waves {
  int reflectivity;
  int spinflip;
  std::vector<wave> waves;

  coherent_waves() {}
  coherent_waves(const coherent_waves& o) { reflectivity = o.reflectivity; spinflip = o.spinflip; waves = o.waves; }

  std::complex<double> sum(const std::vector<double>& x, const event& e) const;
  std::vector<wave>& getWaves() { return waves; }
  const std::vector<wave>& getWaves() const { return waves; }
  size_t getNwaves() const { return waves.size(); }
};

struct waveset : public std::vector<coherent_waves> {
public:
  size_t getNwaves() const {
    size_t count = 0;
    for (size_t i = 0; i < this->size(); i++)
      count += (*this)[i].waves.size();
    return count;
  }
};

#endif
