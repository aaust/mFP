#ifndef LIKELIHOOD_H__
#define LIKELIHOOD_H__

#include <vector>
#include <map>

#include "Minuit2/FCNBase.h"
#include "TMatrixD.h"  // or any other two-index object

#include "wave.h"
#include "event.h"
#include "eventStream.h"

class likelihood : public ROOT::Minuit2::FCNBase {
  std::vector<coherent_waves> ws;
  eventStream* RDevents;
  eventStream* MCevents;

  size_t nBins, lastIdx;
  double threshold;
  double binWidth;

  size_t idxBranching;

  size_t currentBin;

  mutable std::map<int, double> weights;
  mutable size_t eventsAccepted;
  mutable size_t eventsInBin;

  TMatrixD MCweights;

public:
  likelihood(const waveset& ws_,
	     eventStream* RDevents_,
	     eventStream* MCevents_,
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

  void
  calcMCweights();

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

  std::complex<double>
  calcMoment(int L, int M) const;

  void
  setBin(size_t iBin) { currentBin = iBin;
    massLow = threshold + iBin*binWidth;
    massHigh = threshold + (iBin+1)*binWidth;
    RDevents->setBin(iBin);
    MCevents->setBin(iBin);
    cout << "setting bin" << endl;
    calcMCweights();
  }

  size_t
  countEventsInBin() const {
    return MCevents->countEvents(true, true);
  }

  void
  clearWeights() { weights.clear(); MCweights.Clear(); cout << "cleared weights." << endl; }
};

#endif
