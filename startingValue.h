#ifndef STARTINGVALUE_H__
#define STARTINGVALUE_H__

#include <string>

#include "Rtypes.h" // for Double_t

struct tStartingValue {
  std::string name;
  Double_t value;
  bool fixed;

  tStartingValue(const std::string& n, Double_t v, bool f)
    : name(n), value(v), fixed(f)
  {}
};

#endif
