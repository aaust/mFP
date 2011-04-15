#ifndef LIKELIHOOD_H__
#define LIKELIHOOD_H__

#include <vector>
#include <map>

#include "Minuit2/FCNBase.h"

#include "wave.h"
#include "event.h"

class likelihood : public ROOT::Minuit2::FCNBase {
  std::vector<coherent_waves> ws;
  std::vector<event> RDevents;
  std::vector<event> MCevents;
  std::vector<event> MCallEvents;

  std::vector<std::vector<event> > binnedRDevents;
  std::vector<std::vector<event> > binnedMCevents;
  std::vector<double> binnedEtaAcc;

  size_t nBins;
  double threshold;
  double binWidth;

  size_t currentBin;

  mutable std::map<int, double> weights;

public:
  likelihood(waveset ws_,
	     std::vector<event>& RDevents_,
	     std::vector<event>& MCevents_,
	     std::vector<event>& MCallEvents_,
	     size_t nBins_, double threshold_, double binWidth_);

  double Up() const { return 0.5; }

public:
  double
  decay(int reflectivity, int l, int m, double theta, double phi) const;

  double
  probabilityDensity(const std::vector<double>& x, double theta, double phi) const;

  double
  probabilityDensity(const std::vector<double>& x, const event& e) const;

  double
  MCweight(int reflectivity, const wave& w1, const wave& w2) const;

  double
  calc_mc_likelihood(const std::vector<double>& x) const;

  double
  calc_rd_likelihood(const std::vector<double>& x) const;

  double
  calc_likelihood(const std::vector<double>& x) const;

  double
  operator()(const std::vector<double>& x) const;

  void
  setBin(size_t iBin) { currentBin = iBin;
    massLow = threshold + iBin*binWidth;
    massHigh = threshold + (iBin+1)*binWidth;
  }

  void
  clearWeights() { weights.clear(); }
};

#endif