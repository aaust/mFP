#include <vector>
#include <string>

using namespace std;

#include "fitInfo.h"

fitInfo::fitInfo(const string& modelName_,
		 const waveset& ws_,
		 const vector<tStartingValue>& fitVars,
		 size_t nBins_, double threshold_,
		 double binWidth_, size_t nVars_)
  : modelName(modelName_), ws(ws_), nBins(nBins_),
    threshold(threshold_), binWidth(binWidth_), nVars(nVars_)
{
  for (size_t i = 0; i < fitVars.size(); i++)
    {
      paramNames.push_back(fitVars[i].name);
      fixed.push_back(fitVars[i].fixed);
    }
}

#ifndef __CINT__
ClassImp(fitInfo)
#endif
