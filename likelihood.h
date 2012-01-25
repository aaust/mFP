#ifndef LIKELIHOOD_H__
#define LIKELIHOOD_H__

#include <vector>
#include <map>

#include "Minuit2/FCNBase.h"

#include "wave.h"
#include "event.h"

class TH2;

class likelihood : public ROOT::Minuit2::FCNBase {
  std::vector<coherent_waves> ws;

  std::vector<std::vector<event> > binnedRDevents;
  std::vector<std::vector<event> > binnedMCevents;
  std::vector<event> flatMCevents;
  std::vector<double> binnedEtaAcc;

  size_t nBins;
  double threshold;
  double binWidth;

  size_t idxBranching;

  size_t currentBin;

  mutable std::map<int, double> weights;

public:
  likelihood(waveset ws_,
	     std::vector<event>& RDevents,
	     std::vector<event>& MCevents,
	     std::vector<event>& MCallEvents,
	     size_t nBins_, double threshold_, double binWidth_,
	     size_t idxBranching);

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
  
  void
  fillPredict(const std::vector<double>& x, TH2* hth, TH2* hph) const;

  double
  operator()(const std::vector<double>& x) const;

  std::complex<double>
  calcMoment(int L, int M) const;

  void
  setBin(size_t iBin) { currentBin = iBin;
    massLow = threshold + iBin*binWidth;
    massHigh = threshold + (iBin+1)*binWidth;
  }

  size_t
  eventsInBin() const { return binnedRDevents[currentBin].size(); }

  void
  clearWeights() { weights.clear(); }

private:
  // Not implemented, used to check that likelihoods are not copied,
  // which previously happened unintendedly.
  likelihood& operator=(const likelihood&);
};

#endif
