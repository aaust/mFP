#include <vector>

#include "fitInfo.h"

fitInfo::fitInfo(const std::vector<tStartingValue>& fitVars, size_t nBins_, size_t nVars_)
  : nBins(nBins_), nVars(nVars_)
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
