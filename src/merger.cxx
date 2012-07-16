#include "merger.hxx"

Bool_t checkFile(const char *fileName)
{
	TFile* file = new TFile(fileName);
	if(file->IsZombie())
		return kFALSE;
	TH1D* cdenn1da = (TH1D*)  file->Get("cdenn1da");
	if(cdenn1da->Integral(1,TMath::Nint((double)cdenn1da->GetMaximumBin()/5)) > 0)
		return kTRUE;
	else
		return kFALSE;
}


Bool_t mergeFiles(std::vector<ktBinFile> &inputFiles, ktBinFile &outputFile)
{
	TFileMerger merger(kFALSE,kFALSE);
	merger.SetMsgPrefix("[merger]");
	merger.SetNotrees(kFALSE);

	if (!merger.OutputFile(boost::get<0>(outputFile).c_str(),kTRUE,1))
	{
		std::cerr << "Could not open output file: " << boost::get<0>(outputFile) << std::endl;
		return kFALSE;
	}
	for(std::vector<ktBinFile>::iterator it = inputFiles.begin() ; it != inputFiles.end() ; ++it)
		if(!merger.AddFile(boost::get<0>(*it).c_str()))
		{
			std::cerr << "could not add file: " << boost::get<0>(*it).c_str() << std::endl;
			return kFALSE;
		}
	boost::get<1>(outputFile) = boost::get<1>(inputFiles.front());
	boost::get<2>(outputFile) = boost::get<2>(inputFiles.back());
	return merger.Merge();
}

Bool_t decodeKtBin(const unsigned int binNumber, double &ktMin, double &ktMax)
{
	switch(binNumber)
	{
		case 51:
			ktMin = 0.12;
			ktMax = 0.2;
			break;
		case 52:
			ktMin = 0.2;
			ktMax = 0.3;
			break;
		case 53:
			ktMin = 0.3;
			ktMax = 0.4;
			break;
		case 54:
			ktMin = 0.4;
			ktMax = 0.5;
			break;
		case 55:
			ktMin = 0.5;
			ktMax = 0.6;
			break;
		case 56:
			ktMin = 0.6;
			ktMax = 0.7;
			break;
		case 57:
			ktMin = 0.7;
			ktMax = 0.8;
			break;
		case 58:
			ktMin = 0.8;
			ktMax = 1.0;
			break;
		case 59:
			ktMin = 1.0;
			ktMax = 1.2;
			break;
		default:
			std::cerr << "\tUnknown kt bin: " << binNumber << std::endl;
			return kFALSE;
			break;
	}
	return kTRUE;
}

Bool_t parseFileName(const char *fileName, double &ktMin, double &ktMax)
{
	unsigned int ktBin;
	ktMin = 0;
	ktMax = 0;
	boost::cmatch matches;
	boost::regex pattern("(.*?)(\\d+)(.*?)");
	boost::regex_match(basename((char*)fileName), matches, pattern);
	ktBin = boost::lexical_cast<unsigned int>(matches[2]);
	return decodeKtBin(ktBin, ktMin, ktMax);
}

Bool_t loadFileList(std::istream &inputFilesList, std::vector<ktBinFile> &inputFiles)
{
	char buffer[256];
	std::string tmpfile(TMP_FILE);
	ktBinFile outfile;
	static int counter = 0;
	std::vector<ktBinFile> filesToMerge;
	double ktBinMin, ktBinMax;
  	while (inputFilesList.good())
    {
    	std::fill(buffer, buffer+255, 0);
    	inputFilesList >> buffer;
    	if(buffer[0] == '\0')
    		continue;
   		std::cerr << buffer;
		if(checkFile(buffer))
    	{
			std::cerr << "\t\t[ OK ]" << std::endl;
    		if(filesToMerge.size() == 0)
			{	
				if(parseFileName(buffer, ktBinMin, ktBinMax))
					inputFiles.push_back(boost::make_tuple(std::string(buffer), ktBinMin, ktBinMax));
			}
			else
			{// merging
				std::cerr << "\t[ MERGING IN PROGRESS... ]" << std::endl;

				if(parseFileName(buffer, ktBinMin, ktBinMax))
					filesToMerge.push_back(boost::make_tuple(std::string(buffer), ktBinMin, ktBinMax));
				while(true)
				{
					boost::get<0>(outfile) = tmpfile;
					boost::get<0>(outfile) += boost::lexical_cast<std::string>(counter++);
					boost::get<0>(outfile) += ".root";
					ifstream testfile(boost::get<0>(outfile).c_str());
					if(testfile.good())
						continue;
					else
						break;
				}
				if(!mergeFiles(filesToMerge,outfile))
				{
					std::cerr << "Could not merge files: " << std::endl;
					for(std::vector<ktBinFile>::iterator it = filesToMerge.begin() ; it != filesToMerge.end() ; ++it)
						std::cerr << "\t" << boost::get<0>(*it) << std::endl;
					return kFALSE;
				}
				std::cerr << "Files merged into: " << boost::get<0>(outfile) << std::endl;
				inputFiles.push_back(outfile);
				filesToMerge.clear();
			}
		}
		else
		{	
			parseFileName(buffer, ktBinMin, ktBinMax);
			filesToMerge.push_back(boost::make_tuple(std::string(buffer), ktBinMin, ktBinMax));
			std::cerr << "\t\t[ REQUIRES MERGING ]" << std::endl;
		}
    }
    return kTRUE;
}


int main(int argc, char **argv)
{
	std::vector<ktBinFile> files;
	if(argc < 1)
	{	std::cerr << argv[0] << " " << "pp" << std::endl <<
				"Pass to the argument pair identifier and file list to stdin" << std::endl;
		return 1;
	}
	if(loadFileList(std::cin, files))
	{
		for(int i = 0, size = files.size(); i < size; ++i)
			std::cout << boost::get<0>(files[i]) << " " << boost::get<1>(files[i]) << " " << boost::get<2>(files[i]) << std::endl;
		return 0;
	}
	else
		std::cerr << "Could not load file list" << std::endl;
	return 1;
}