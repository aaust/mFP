#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "control.h"

using namespace std;

vector<string> dataFiles;
bool flatMC;
vector<string> MCFiles;
double threshold;
int nBins;
double binWidth;
int nFits;
bool continuous;
bool ambiguous;
bool predict;
string modelName;


static bool
readString(FILE* fd, string& result)
{
  char *text;
  if (fscanf(fd, "%as", &text) != 1) // %as is a GNU extension
    return false;
  result = string(text);
  free(text);
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
  bool haveFlatMC = false;
  bool haveMCFile = false;
  bool haveNFits = false;
  bool haveThreshold = false;
  bool haveNBins = false;
  bool haveBinWidth = false;
  bool haveContinuous = false;
  bool haveAmbiguous = false;
  bool havePredict = false;
  bool haveModelName = false;
  while(!feof(fd))
    {
      char key[999];
      if (fscanf(fd, "%s", key) != 1)  // possible overflow here
	break;

      if (!strcasecmp(key, "datafile"))
	{
	  string s;
	  if (!readString(fd, s))
	    break;
	  dataFiles.push_back(s);
	  haveDataFile = true;
	}
      else if (!strcasecmp(key, "flatMC"))
	{
	  haveFlatMC = true;
	}
      else if (!strcasecmp(key, "mcfile"))
	{
	  string s;
	  if (!readString(fd, s))
	    break;
	  MCFiles.push_back(s);
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
      else if (!strcasecmp(key, "continuous"))
	{
	  haveContinuous = true;
	}
      else if (!strcasecmp(key, "ambiguous"))
	{
	  haveAmbiguous = true;
	}
      else if (!strcasecmp(key, "predict"))
	{
	  havePredict = true;
	}
      else if (!strcasecmp(key, "modelName"))
	{
	  if (!readString(fd, modelName))
	    break;
	  haveModelName = true;
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
  if (!haveFlatMC && !haveMCFile)
    {
      cerr << "No MC file given." << endl;
      good = false;
    }
  if (!haveFlatMC && MCFiles.size() != dataFiles.size())
    {
      cerr << "different number of MC files and data files." << endl;
      good = false;
    }
  flatMC = haveFlatMC;
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
  continuous = haveContinuous;
  ambiguous = haveAmbiguous;
  predict = havePredict;
  return good;
}

