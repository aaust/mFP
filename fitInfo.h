#ifndef FITINFO_H__
#define FITINFO_H__

//
// The data stored in the output tree's UserInfo need to be put in a
// TList, which in turn means they have to be put in data structures
// inheriting from TObject, so we take care of this here.

#include <vector>
#include <string>

#include "TObject.h"
#include "TString.h"

#include "startingValue.h"

class fitInfo : public TObject
{
public:
  fitInfo() {};
  fitInfo(const std::vector<tStartingValue>& fitVars);
private:
  std::vector<std::string> paramNames;

  ClassDef(fitInfo, 1)
};

#endif
