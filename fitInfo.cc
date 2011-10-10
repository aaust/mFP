#include <vector>

#include "fitInfo.h"

fitInfo::fitInfo(const std::vector<tStartingValue>& fitVars)
{
  for (size_t i = 0; i < fitVars.size(); i++)
    {
      paramNames.push_back(fitVars[i].name.c_str());
    }
}

#ifndef __CINT__
ClassImp(fitInfo)
#endif
