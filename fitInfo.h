#ifndef FITINFO_H__
#define FITINFO_H__

//
// The data stored in the output tree's UserInfo need to be put in a
// TList, which in turn means they have to be put in data structures
// inheriting from TObject, so we take care of this here.
//
// Collect relevant global information concerning the fit.

#include <vector>
#include <string>

#include "TObject.h"
#include "TString.h"

#include "wave.h"
#include "startingValue.h"

class fitInfo : public TObject
{
public:
  fitInfo() {};
  fitInfo(const std::string& modelName_,
	  const waveset& ws_,
	  const std::vector<tStartingValue>& fitVars,
	  size_t nBins_, double threshold_,
	  double binWidth_, size_t nVars_);
  ~fitInfo();

  size_t getNvars() const { return nVars; }

  const std::string& getModelName() const { return modelName; }
  const waveset& getWaveSet() const { return ws; }
  const std::vector<std::string>& getParamNames() const { return paramNames; }
  const std::vector<bool>& getFixed() const { return fixed; }
  size_t getNbins() const { return nBins; }
  double getThreshold() const { return threshold; }
  double getBinWidth() const { return binWidth; }
  double getLower() const { return threshold; }
  double getUpper() const { return threshold + nBins*binWidth; }
private:
  std::string modelName;
  waveset ws;
  std::vector<std::string> paramNames;
  std::vector<bool> fixed;
  size_t nBins;
  double threshold, binWidth;
  size_t nVars;

  ClassDef(fitInfo, 2)
};

#endif
