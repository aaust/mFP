// This is of course absurdely inefficient, giving each factorized
// number a vector of powers would immediately give the multiplication
// etc routines all they need to not scan through all nPrime elements
// everytime an operation is performed.  Also, the sign belongs in the
// fraction, the current handling is more than awkward.  Anyway, it
// works.

#include <iostream>
#include <math.h>

using namespace std;

// Prime number list taken from http://www.mathematical.com/primes0to1000k.html
static const int primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
			      103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
			      199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
			      313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
			      433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
			      563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
			      673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
			      811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
			      941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,
			      1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,
			      1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,
			      1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,
			      1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,
			      1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,
			      1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,
			      1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
			      1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,
			      1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,
			      2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,
			      2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,
			      2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,
			      2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,
			      2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,
			      2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,
			      2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,
			      2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,
			      3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,
			      3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,
			      3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,
			      3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,
			      3527,3529,3533,3539,3541,3547,3557,3559,3571,3581,3583,3593,3607,3613,3617,
			      3623,3631,3637,3643,3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,3733,
			      3739,3761,3767,3769,3779,3793,3797,3803,3821,3823,3833,3847,3851,3853,3863,
			      3877,3881,3889,3907,3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,4001,
			      4003,4007,4013,4019,4021,4027,4049,4051,4057,4073,4079,4091,4093,4099,4111,
			      4127,4129,4133,4139,4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,4241,
			      4243,4253,4259,4261,4271,4273,4283,4289,4297,4327,4337,4339,4349,4357,4363,
			      4373,4391,4397,4409,4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,4507,
			      4513,4517,4519,4523,4547,4549,4561,4567,4583,4591,4597,4603,4621,4637,4639,
			      4643,4649,4651,4657,4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,4759,
			      4783,4787,4789,4793,4799,4801,4813,4817,4831,4861,4871,4877,4889,4903,4909,
			      4919,4931,4933,4937,4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,5009,
			      5011,5021,5023,5039,5051,5059,5077,5081,5087,5099,5101,5107,5113,5119,5147,
			      5153,5167,5171,5179,5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,5281,
			      5297,5303,5309,5323,5333,5347,5351,5381,5387,5393,5399,5407,5413,5417,5419,
			      5431,5437,5441,5443,5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,5527,
			      5531,5557,5563,5569,5573,5581,5591,5623,5639,5641,5647,5651,5653,5657,5659, };
static const int nPrimes = sizeof(primes) / sizeof(int);

class factoredInteger {
  int value;
  int powerDecomposition[nPrimes];
  int sign;

public:
  factoredInteger(int number)
  {
    int i = value = number;

    memset (powerDecomposition, 0, sizeof(int)*nPrimes);

    sign = 0;
    if (i == 0)
      return;

    sign = i / abs(i);
    i = i*sign;

    for (int n = nPrimes - 1; n >= 0 && i != 1; n--)
      {
	while ((i % primes[n]) == 0)
	  {
	    powerDecomposition[n]++;
	    i /= primes[n];
	  }
      }
    if (i != 1)
      {
	cout << "can't decompose " << number << endl;
      }
  }

  factoredInteger(const int decomp[nPrimes], int sgn = 1)
    : sign(sgn)
  {
    value = sgn;
    for (int i = 0; i < nPrimes; i++)
      {
	powerDecomposition[i] = decomp[i];
	value *= pow(primes[i], powerDecomposition[i]);
      }
    if (sign == 0)
      value =0;
  }

  factoredInteger
  operator*(int i) const
  {
    return (*this)*factoredInteger(i);
  }

  factoredInteger
  operator*(const factoredInteger& o) const
  {
    int decomp[nPrimes];
    for (int i = 0; i < nPrimes; i++)
      {
	decomp[i] = o.powerDecomposition[i] + powerDecomposition[i];
      }
    return factoredInteger(decomp, o.sign * sign);
  }

  factoredInteger
  operator+(const factoredInteger& o) const
  {
    return factoredInteger(this->value + o.value);
  }

  // this modifies both the argument and this.
  void
  divideByGCD(factoredInteger& o)
  {
    if (o.sign == 0 || this->sign == 0)
      return;

    if (o.sign != this->sign)
      {
	this->sign = -1;
	o.sign = 1;
      }
    else
      {
	this->sign = o.sign = 1;
      }
    for (int i = 0; i < nPrimes; i++)
      {
	int pow = std::min(powerDecomposition[i], o.powerDecomposition[i]);
	powerDecomposition[i] -= pow;
	o.powerDecomposition[i] -= pow;
      }

    fixValue();
    o.fixValue();
  }

  void
  fixValue()
  {
    value = sign;
    for (int i = 0; i < nPrimes; i++)
      {
	value *= pow(primes[i], powerDecomposition[i]);
      }
  }

  int
  getValue() const
  {
    return value;
  }

  int
  getSign() const
  {
    return sign;
  }

  void
  write() const
  {
    if (sign == -1)
      cout << "-";
    if (sign == 0)
      {
	cout << 0 << endl;
	return;
      }
    bool first = true;
    for (int i = 0; i < nPrimes; i++)
      {
	if (powerDecomposition[i] != 0)
	  {
	    if (!first)
	      cout << " * ";
	    else
	      first = false;
	    cout << primes[i] << "^" << powerDecomposition[i];
	  }
      }
    if (first)
      cout << "1";
    cout << endl;
  }
};

class fraction {
  factoredInteger numerator, denominator;

public:
  fraction(int num, int denom)
    : numerator(num), denominator(denom)
  {
    this->reduce();
  }

  fraction(const factoredInteger& num, const factoredInteger& denom)
    : numerator(num), denominator(denom)
  {
    this->reduce();
  }

  const fraction&
  reduce()
  {
    numerator.divideByGCD(denominator);
    return *this;
  }

  fraction
  operator*(const factoredInteger& other)
  {
    return fraction(this->numerator * other, this->denominator);
  }

  fraction
  operator*(const fraction& other)
  {
    // reduce
    fraction f1(this->numerator, other.denominator);
    fraction f2(other.numerator, this->denominator);

    // calculate with reduced values.
    return fraction(f1.numerator*f2.numerator, f1.denominator*f2.denominator);
  }

  fraction
  operator+(const fraction& other)
  {
    factoredInteger newDenom(this->denominator * other.denominator);

    return fraction(this->numerator * other.denominator + other.numerator * this->denominator, newDenom);
  }

  int
  getSign() const
  {
    return numerator.getSign();
  }    

  void
  write() const
  {
    cout << numerator.getValue() << " / " << denominator.getValue() << endl;
  }
};


namespace {
  long
  fac(long n)
  {
    long result = 1;
    if (n > 1)
      for (long i = 1; i <= n; i++)
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

  fraction
  delta_square(int j1, int j2, int j3)
  {
    return fraction(fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3),
		    fac(j1 + j2 +j3 + 1));
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

  fraction
  prefactor_square(long j1, long j2, long j3, long m1, long m2, long m3)
  {
    int sign = ((j1 - j2 - m3) & 0x1) ? -1 : 1;  // odd / even
    fraction result = fraction(delta_square(j1,j2,j3) * factoredInteger(fac(j1+m1)*fac(j1-m1)*fac(j2+m2)
						   *fac(j2-m2)*fac(j3+m3)*fac(j3-m3))) * sign;

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
      double add = (pow(-1,s)
	      / fac(s) / fac(j1+j2-j3-s)
	      / fac(j1-m1-s) / fac(j2+m2-s)
	      / fac(j3-j2+m1+s) / fac(j3-j1-m2+s));
      sum+=add;
    }

  return prefactor(j1,j2,j3,m1,m2,m3)*sum;
}


// This returns the signed square of the 3j-symbol.
fraction
threeJalgebraically(long j1, long j2, long j3, long m1, long m2, long m3)
{
  fraction zero(0,1);

  if (j1 < 0 || j2 < 0 || j3 < 0)
    return zero;
  if (j1 > j2 + j3)
    return zero;
  if (j2 > j3 + j1)
    return zero;
  if (j3 > j1 + j2)
    return zero;
  if (abs(j1 - j2) > j3)
    return zero;
  if (abs(j2 - j3) > j1)
    return zero;
  if (abs(j3 - j1) > j2)
    return zero;
  if (m1 + m2 + m3 != 0)
    return zero;
  if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
    return zero;

  // The 3j-symbol is != 0.  Find the valid range for the summation in
  // http://dlmf.nist.gov/34.2 (34.2.4)
  long minS = 0;
  minS = max(minS, -(j3 - j2 + m1));
  minS = max(minS, -(j3 - j1 - m2));
  long maxS = j1 + j2 - j3;  // needs refinement
  maxS = min(maxS, j1 - m1);
  maxS = min(maxS, j2 + m2);

  if (minS > maxS)
    return zero;

  fraction sum(0, 1);
  for (int s = minS; s <= maxS; s++)
    {
      fraction add((fraction(1, fac(s)) * fraction(1, fac(j1+j2-j3-s))
		   * fraction(1, fac(j1-m1-s)) * fraction(1, fac(j2+m2-s))
		    * fraction(1, fac(j3-j2+m1+s)) * fraction(1, fac(j3-j1-m2+s)))*pow(-1,s));
      sum = sum + add;
    }

  return prefactor_square(j1,j2,j3,m1,m2,m3)*sum*sum*sum.getSign();
}
