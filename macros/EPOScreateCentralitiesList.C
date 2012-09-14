#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <TChain.h>

using namespace std;

typedef struct
{
	string	name;
	float bMin;
	float bMax;
} Centrality;

const int nCentralities = 9;

// int main()
int EPOScreateCentralitiesList()
{
	char buffer[256];
	int iEvent, totalEventCount;
	Float_t bim;
	Centrality tmp;

	Centrality centralities[nCentralities];
	fstream centralitiesFiles[nCentralities];

	centralities[0].name = "0005";
	centralities[0].bMin = 0;
	centralities[0].bMax = 3.375;

	centralities[1].name = "0510";
	centralities[1].bMin = 3.375;
	centralities[1].bMax = 4.875;

	centralities[2].name = "1020";
	centralities[2].bMin = 4.875;
	centralities[2].bMax = 6.875;

	centralities[3].name = "2030";
	centralities[3].bMin = 6.875;
	centralities[3].bMax = 9.125;

	centralities[4].name = "3040";
	centralities[4].bMin = 9.125;
	centralities[4].bMax = 9.875;

	centralities[5].name = "4050";
	centralities[5].bMin = 9.875;
	centralities[5].bMax = 11.125;

	centralities[6].name = "5060";
	centralities[6].bMin = 11.125;
	centralities[6].bMax = 12.125;

	centralities[7].name = "7080";
	centralities[7].bMin = 12.125;
	centralities[7].bMax = 13.875;

	centralities[8].name = "8099";
	centralities[8].bMin = 13.875;
	centralities[8].bMax = 20;



	TChain *tEPOS = new TChain("teposevent");

	fstream fileList("eventFiles.list", fstream::out);
	while (cin.good())
    {
    	std::fill(buffer, buffer+255, 0);
    	cin >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	tEPOS->Add(buffer);
    	fileList << buffer << endl;
   	}
	fileList.close();
	cout << "Event list 'eventFiles.list' closed." << endl;

	for(int i = 0; i < nCentralities; ++i)
		centralitiesFiles[i].open((centralities[i].name + string(".centrality")).c_str(), fstream::out);

	tEPOS->SetBranchStatus("*",1);
	tEPOS->SetBranchAddress("bim", &bim);
	totalEventCount = tEPOS->GetEntries();

	cout << "Event files loaded. Total << " << totalEventCount << " events. Beginning classification." << endl;

	for(iEvent = 0; iEvent < totalEventCount; ++iEvent)
	{
		if(iEvent % 300 == 0)
			cout << iEvent << "/" << totalEventCount << endl;
		tEPOS->GetEntry(iEvent);
		for(int i = 0; i < nCentralities; ++i)
		{
			if(centralities[i].bMin <= bim && centralities[i].bMax > bim )
			{
				centralitiesFiles[i] << iEvent << endl;
				// cout << bim << " @ " << centralities[i].bMin << "-" << centralities[i].bMax << endl;
				break;
			}
		}
	}

	cout << "Closing centrality files" << endl;
	// closing files
	for(int i = 0; i < nCentralities; ++i)
		centralitiesFiles[i].close();

	return 0;
}
