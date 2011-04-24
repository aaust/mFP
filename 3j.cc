// This is of course absurdely inefficient, giving each factorized
// number a vector of powers would immediately give the multiplication
// etc routines all they need to not scan through all nPrime elements
// everytime an operation is performed.  Also, the sign belongs in the
// fraction, the current handling is more than awkward.  Anyway, it
// works.

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

#include "wave.h"

namespace {
  long long
  fac(long long n)
  {
    long long result = 1;
    if (n > 1)
      for (int i = 1; i <= n; i++)
	result *= i;
    return result;
  }


  // Calculate the 3j-symbol from the explicit expression given in http://dlmf.nist.gov/34.2
  // First their delta (34.2.5):

  double
  delta(long j1, long j2, long j3)
  {
    // Valid argument range is assumed.  This is assured in the 3j function.
    double result = (sqrt(fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3))
		     / sqrt(fac(j1 + j2 + j3 + 1)));
    //cout << "delta = " << result << endl;
    return result;
  }


  // The complete prefactor in (34.2.4).  Again, rangechecks are left to the main function.
  double
  prefactor(long j1, long j2, long j3, long m1, long m2, long m3)
  {
    int sign = ((j1 - j2 - m3) & 0x1) ? -1 : 1;  // odd / even
    double result = sign*delta(j1, j2, j3)*sqrt(fac(j1+m1)*fac(j1-m1)*fac(j2+m2)*fac(j2-m2)*fac(j3+m3)*fac(j3-m3));
    //cout << "prefactor = " << result << endl;
    return result;
  }
}


double
threeJ(long j1, long j2, long j3, long m1, long m2, long m3)
{
  if (j1 < 0 || j2 < 0 || j3 < 0)
    return 0;
  if (j1 > j2 + j3)
    return 0;
  if (j2 > j3 + j1)
    return 0;
  if (j3 > j1 + j2)
    return 0;
  if (abs(j1 - j2) > j3)
    return 0;
  if (abs(j2 - j3) > j1)
    return 0;
  if (abs(j3 - j1) > j2)
    return 0;
  if (m1 + m2 + m3 != 0)
    return 0;
  if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
    return 0;

  // The 3j-symbol is != 0.  Find the valid range for the summation in loc.cit. (34.2.4)
  long minS = 0;
  minS = max(minS, -(j3 - j2 + m1));
  minS = max(minS, -(j3 - j1 - m2));
  long maxS = j1 + j2 - j3;  // needs refinement
  maxS = min(maxS, j1 - m1);
  maxS = min(maxS, j2 + m2);

  if (minS > maxS)
    return 0;

  double sum = 0;
  for (int s = minS; s <= maxS; s++)
    {
      double add = ((s & 0x1 ? -1. : 1.)
		    / fac(s) / fac(j1+j2-j3-s)
		    / fac(j1-m1-s) / fac(j2+m2-s)
		    / fac(j3-j2+m1+s) / fac(j3-j1-m2+s));
      //cout << add << endl;
      sum+=add;
    }

  return prefactor(j1,j2,j3,m1,m2,m3)*sum;
}


double
theta(int m)
{
  if (m == 0)
    return .5;
  else
    return sqrt(.5);
}


// Get decomposition of H(LM) in terms of the waveset ws.
void
decomposeMoment(int L, int M, const waveset& ws)
{
  for (size_t iWs = 0; iWs < ws.size(); iWs++)
    {
      int eps = ws[iWs].reflectivity;

      const vector<wave>& w = ws[iWs].waves;
      for (size_t iW1 = 0; iW1 < w.size(); iW1++)
	{
	  const wave& w1 = w[iW1];
	  for (size_t iW2 = 0; iW2 < w.size(); iW2++)
	    {
	      const wave& w2 = w[iW2];
	      double threeJ1 = threeJ(L, w1.l, w2.l, 0, 0, 0);
	      if (threeJ1 == 0)
		continue;

	      int sign = (w1.m + M & 0x1 ? -1 : 1);
	      int sign1 = (w1.m + 1 & 0x1 ? -1 : 1);
	      int sign2 = (w2.m + 1 & 0x1 ? -1 : 1);
	      int sign12 = (w1.m + w2.m & 0x1 ? -1 : 1);
	      double parentheses = (threeJ(L, w1.l, w2.l, -M, -w1.m, w2.m)
				      + threeJ(L, w1.l, w2.l, -M, w1.m, w2.m) * eps * sign1
				      + threeJ(L, w1.l, w2.l, -M, -w1.m, -w2.m) * eps * sign2
				      + threeJ(L, w1.l, w2.l, -M, w1.m, -w2.m) * sign12) * sign;
	      if (parentheses == 0)
		continue;

	      cout << " + "
		   << (theta(w1.m)*theta(w2.m)*sqrt((2*w1.l+1)*(2*w2.l+1))
		       *threeJ1*parentheses)
		   << "*rho(eps = " << eps << ", " << w1.l << ", " << w1.m << ", " << w2.l << ", " << w2.m
		   << ")";
	    }
	}
    }
  cout << endl;
}
