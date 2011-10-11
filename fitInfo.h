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

#include "startingValue.h"

class fitInfo : public TObject
{
public:
  fitInfo() {};
  fitInfo(const std::vector<tStartingValue>& fitVars, size_t nBins_);

  size_t getNparams() const { return paramNames.size(); }

  const std::vector<std::string>& getParamNames() const { return paramNames; }
  const std::vector<bool>& getFixed() const { return fixed; }
  size_t getNbins() const { return nBins; }
private:
  std::vector<std::string> paramNames;
  std::vector<bool> fixed;
  size_t nBins;

  ClassDef(fitInfo, 1)
};

#endif
