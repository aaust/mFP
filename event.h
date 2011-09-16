#ifndef EVENT_H__
#define EVENT_H__

#include <map>

#include "wave.h"

extern double massLow, massHigh;

class event {
public:
  double mass;
  double tPrime;
  double theta;
  double phi;

  event() { mass = tPrime = theta = phi = 0; }
  /*
  event(const event& other) { theta = other.theta; phi = other.phi; }
  */
  event(double th, double ph) { mass = tPrime = 0; theta = th; phi = ph; }
  event(double mass_, double tPrime_, double th, double ph)
  { mass = mass_; tPrime = tPrime_; theta = th; phi = ph; }
  event(const event& o) { mass = o.mass; tPrime = o.tPrime; theta = o.theta; phi = o.phi; }

  double decayAmplitude(int reflectivity, const wave& w) const
  { return decayAmplitude(reflectivity, w.l, w.m); };
  double decayAmplitude(int reflectivity, int l, int m) const;
  double MCweight(int reflectivity, const wave& w1, const wave& w2) const
  { return MCweight(reflectivity, w1.l, w1.m, w2.l, w2.m); }
  double MCweight(int reflectivity, int l1, int m1, int l2, int m2) const;
  std::complex<double> momentWeight(int L, int M) const;

  bool accepted() const {
    return (this->mass >= massLow && this->mass < massHigh
	    //&& this->tPrime > 0.1 && this->tPrime < 0.3
	    //&& this->tPrime > 0.1
	    //&& this->tPrime > 0.1
	    //&& cos(theta) > -0.7
	    && 1);
  }
};

#endif
