#ifndef WAVE_H__
#define WAVE_H__

#include <complex>
#include <vector>

using namespace std;

struct wave {
  string name;
  int l, m;

  wave() { }
  wave(int ll, int mm) { l = ll; m = mm; }
  wave(const char* name_, int ll, int mm) : name(name_), l(ll), m(mm) {}
  wave(const wave& o) { l = o.l; m = o.m; name = o.name; }
};

class event;

struct coherent_waves {
  int reflectivity;
  int spinflip;
  std::vector<wave> waves;

  coherent_waves() {}
  coherent_waves(const coherent_waves& o) { reflectivity = o.reflectivity; spinflip = o.spinflip; waves = o.waves; }

  std::complex<double> sum(const double* x, const event& e) const;
  size_t getNwaves() { return waves.size(); }
};

typedef std::vector<coherent_waves> waveset;

#endif
