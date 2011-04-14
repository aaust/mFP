#include <iostream>
#include <stdio.h>
#include <string.h>

#include "control.h"

using namespace std;

string dataFile;
bool flatMC;
string MCFile;
double threshold;
int nBins;
double binWidth;
int nFits;


static bool
readString(FILE* fd, string& result)
{
  char text[999];
  if (fscanf(fd, "%s", text) != 1) // possible overflow
    return false;
  result = string(text);
  return true;
}


static bool
readInt(FILE *fd, int& result)
{
  return fscanf(fd, "%d", &result) == 1;
}


static bool
readDouble(FILE* fd, double& result)
{
  return fscanf(fd, "%lf", &result) == 1;
}


bool
readControlFile(string fileName)
{
  FILE* fd = fopen(fileName.c_str(), "r");
  if (!fd)
    {
      cerr << "Can't open input file '" << fileName << "'" << endl;
      return false;
    }

  bool haveDataFile = false;
  bool flatMC = false;
  bool haveMCFile = false;
  bool haveNFits = false;
  bool haveThreshold = false;
  bool haveNBins = false;
  bool haveBinWidth = false;
  while(!feof(fd))
    {
      char key[999];
      if (fscanf(fd, "%s", key) != 1)  // possible overflow here
	break;

      if (!strcasecmp(key, "datafile"))
	{
	  if (!readString(fd, dataFile))
	    break;
	  haveDataFile = true;
	}
      else if (!strcasecmp(key, "flatMC"))
	{
	  flatMC = true;
	}
      else if (!strcasecmp(key, "mcfile"))
	{
	  if (!readString(fd, MCFile))
	    break;
	  haveMCFile = true;
	}
      else if (!strcasecmp(key, "threshold"))
	{
	  if (!readDouble(fd, threshold))
	    break;
	  haveThreshold = true;
	}
      else if (!strcasecmp(key, "nfits"))
	{
	  if (!readInt(fd, nFits))
	    break;
	  haveNFits = true;
	}
      else if (!strcasecmp(key, "nbins"))
	{
	  if (!readInt(fd, nBins))
	    break;
	  haveNBins = true;
	}
      else if (!strcasecmp(key, "binwidth"))
	{
	  if (!readDouble(fd, binWidth))
	    break;
	  haveBinWidth = true;
	}
      else
	{
	  cerr << "unknown key '" << key << "' encountered" << endl;
	  return false;
	}
    }

  bool good = true;
  if (!haveDataFile)
    {
      cerr << "No data file given." << endl;
      good = false;
    }
  if (!flatMC && !haveMCFile)
    {
      cerr << "No MC file given." << endl;
      good = false;
    }
  if (!haveThreshold)
    {
      cerr << "No threshold given." << endl;
      good = false;
    }
  if (!haveNFits)
    {
      cerr << "No number of fits given." << endl;
      good = false;
    }
  if (!haveNBins)
    {
      cerr << "No number of bins given." << endl;
      good = false;
    }
  if (!haveBinWidth)
    {
      cerr << "No bin width given." << endl;
      good = false;
    }
  return good;
}

