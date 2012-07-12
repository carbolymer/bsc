#ifndef _MERGER_HXX_
#define _MERGER_HXX_

#include <fstream>
#include <vector>
#include <TFile.h>
#include <TFileMerger.h>
#include <TH1D.h>
#include <TMath.h>
#include <libgen.h>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#define TMP_FILE "./tmp/outfile_merged"

typedef boost::tuple<std::string,double,double> ktBinFile;

Bool_t checkFile(const char *fileName);
Bool_t mergeFiles(std::vector<ktBinFile> &inputFiles, ktBinFile &outputFile);
void decodeKtBin(const unsigned int binNumber, double &ktMin, double &ktMax);
void parseFileName(const char *fileName, double &ktMin, double &ktMax);
Bool_t loadFileList(std::ifstream &inputFilesList, std::vector<ktBinFile> &inputFiles);



#endif