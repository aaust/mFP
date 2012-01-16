//
// This writes the moments as function of the given waveset.
//

#include <complex>
#include <vector>
#include <map>

using namespace std;

#include "wave.h"
#include "3j.h"

int
main()
{
  vector<wave> positive;
  positive.push_back(wave("P+", 1, 1));
  positive.push_back(wave("D+", 2, 1));
  //positive.push_back(wave("F+", 3, 1));
  positive.push_back(wave("G+", 4, 1));
  //positive.push_back(wave("D++", 2, 2));

  vector<wave> negative;
  negative.push_back(wave("S0", 0, 0));
  negative.push_back(wave("P0", 1, 0));
  negative.push_back(wave("P-", 1, 1));
  negative.push_back(wave("D0", 2, 0));
  negative.push_back(wave("D-", 2, 1));

  coherent_waves wsPos, wsNeg;
  wsPos.reflectivity = +1;
  wsPos.spinflip = +1;
  wsPos.waves = positive;

  wsNeg.reflectivity = -1;
  wsNeg.spinflip = +1;
  wsNeg.waves = negative;

  waveset ws;
  //ws.push_back(wsPos);
  ws.push_back(wsNeg);

  size_t lastIdx = 0;
  for (size_t i = 0; i < ws.size(); i++)
    for (size_t j = 0; j < ws[i].waves.size(); j++, lastIdx += 2)
      ws[i].waves[j].setIndex(lastIdx);

  // Find the non-zero moments for the given waveset.
  std::vector<std::pair<size_t, size_t> > vecMom = listOfMoments(ws);

  // Prepare the moment histograms
  map<std::pair<size_t,size_t>, TH1*> mhMoments;
  for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin();
       it != vecMom.end(); it++)
    {
      cout << "H(" << it->first << ", " << it->second << ") = ";
      for (size_t iCoh = 0; iCoh < ws.size(); iCoh++)
	{
	  const coherent_waves& waves = ws[iCoh];
	  int eps = waves.reflectivity;

	  for (size_t i = 0; i < waves.waves.size(); i++)
	    {
	      const wave& w1 = waves.waves[i];
	      for (size_t j = 0; j <= i; j++)
		{
		  const wave& w2 = waves.waves[j];

		  double coeff = getCoefficient(eps, it->first, it->second, w1.getL(), w1.getM(), w2.getL(), w2.getM());
		  if (coeff == 0)
		    continue;

		  if (j == i)
		    {
		      if (coeff > 0)
			cout << "+ " << coeff;
		      else if (coeff < 0)
			cout << "- " << -coeff;
		      cout << " |" << w1.getName() << "|^2 ";
		    }
		  else
		    {
		      if (coeff > 0)
			cout << "+ " << 2*coeff;
		      else if (coeff < 0)
			cout << "- " << -2*coeff;
		      cout << " Re(" << w1.getName() << "* " << w2.getName() << ") ";
		    }
		}
	    }
	}
      cout << endl;
    }
}
