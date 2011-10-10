#ifndef CONTROL_FILE_H
#define CONTROL_FILE_H

#include <vector>
#include <string>

extern std::vector<std::string> dataFiles;
extern bool flatMC;
extern std::vector<std::string> MCFiles;
extern double threshold;
extern int nBins;
extern double binWidth;
extern int nFits;
extern bool continuous;

bool
readControlFile(std::string fileName);

#endif
