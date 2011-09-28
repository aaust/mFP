#ifndef __GHIST_H__
#define __GHIST_H__

// Allows to do
//  gHist.Fill("histName", "title of histgram", 100, 0, 1,
//             value)
// instead of
//
//  somewhere
//    TH1* hmuell;
//  somewhere else
//    hmuell = new TH1D("histName", "title of histogram", 100, 0, 1);
//  yet somewhere else
//    hmuell->Fill(value);
//
// In order to save some typing when filling the same histogram
// several times, a pointer to the histogram can be obtained via
// getHist.
//
// Some trivialities remain:
//
// TODO: give Fill a return value?  Either return hist or whatever
// TH1::Fill returns
//
// TODO: filling with weights?

#include <map>
#include <stdint.h>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace {
  class gHistClass
  {
  private:
    std::map<uint64_t, TH1*> mHists;

    // Calculate an FNV hash in order to not have to search through
    // the map by comparing strings.  Collisions should be absurdely
    // unlikely.
    uint64_t hashString(const char* str) const {
      const uint64_t FNV_offset_basis = 14695981039346656037ULL;
      const uint64_t FNV_prime = 1099511628211ULL;

      uint64_t hash = FNV_offset_basis;
      const char* p = str;
      while (*p)
	{
	  hash = hash ^ *p;
	  hash = hash * FNV_prime;
	  p++;
	}
      return hash;
    }
  public:
    TH1* getHist(const char* name, const char* title,
		 Int_t nBinsX, Double_t lowX, Double_t highX)
    {
      uint64_t hash = hashString(name);
      if (mHists.find(hash) == mHists.end())
	{
	  mHists[hash] = new TH1D(name, title,
				  nBinsX, lowX, highX);
	}
      return dynamic_cast<TH1*>(mHists[hash]);
    }

    void Fill(const char* name, const char* title,
	      Int_t nBinsX, Double_t lowX, Double_t highX,
	      Double_t valueX, Double_t weight = 1)
    {
      getHist(name, title, nBinsX, lowX, highX)
	->Fill(valueX, weight);
    }

    TH2* getHist(const char* name, const char* title,
		 Int_t nBinsX, Double_t lowX, Double_t highX,
		 Int_t nBinsY, Double_t lowY, Double_t highY)
    {
      uint64_t hash = hashString(name);
      if (mHists.find(hash) == mHists.end())
	mHists[hash] = new TH2D(name, title,
				nBinsX, lowX, highX, nBinsY, lowY, highY);
      return dynamic_cast<TH2*>(mHists[hash]);
    }

    void Fill(const char* name, const char* title,
	      Int_t nBinsX, Double_t lowX, Double_t highX,
	      Int_t nBinsY, Double_t lowY, Double_t highY,
	      Double_t valueX, Double_t valueY, Double_t weight = 1)
    {
      getHist(name, title, nBinsX, lowX, highX, nBinsY, lowY, highY)
	->Fill(valueX, valueY, weight);
    }

    TH3* getHist(const char* name, const char* title,
		 Int_t nBinsX, Double_t lowX, Double_t highX,
		 Int_t nBinsY, Double_t lowY, Double_t highY,
		 Int_t nBinsZ, Double_t lowZ, Double_t highZ)
    {
      uint64_t hash = hashString(name);
      if (mHists.find(hash) == mHists.end())
	mHists[hash] = new TH3D(name, title,
				nBinsX, lowX, highX,
				nBinsY, lowY, highY,
				nBinsZ, lowZ, highZ);
      return dynamic_cast<TH3*>(mHists[hash]);
    }

    void Fill(const char* name, const char* title,
	      Int_t nBinsX, Double_t lowX, Double_t highX,
	      Int_t nBinsY, Double_t lowY, Double_t highY,
	      Int_t nBinsZ, Double_t lowZ, Double_t highZ,
	      Double_t valueX, Double_t valueY, Double_t valueZ, Double_t weight = 1)
    {
      getHist(name, title, nBinsX, lowX, highX, nBinsY, lowY, highY, nBinsZ, lowZ, highZ)
	->Fill(valueX, valueY, valueZ, weight);
    }

    /*
    // Syntactic sugar.
    gHistClass* operator->()
    {
      return this;
    }
    */

    void clear()
    {
      mHists.clear();
    }
  } gHist;
}
#endif
