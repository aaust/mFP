#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H

#include <string>

extern std::string dataFile;
extern bool flatMC;
extern std::string MCFile;
extern double threshold;
extern int nBins;
extern double binWidth;
extern int nFits;

bool
readControlFile(std::string fileName);

#endif
