#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <TH1D.h>
#include <TCanvas.h>
#include <TChain.h>

using namespace std;

const int maxParticleCount = 1e5;
vector<string> inFiles;
TH1D *bDist;

void* worker();

int EPOSimpactDistribution()
{
	char buffer[256];
	stringstream number;

	while (cin.good())
    {
    	std::fill(buffer, buffer+255, 0);
    	cin >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	// cout << "Loading: " << buffer << endl;
    	inFiles.push_back(string(buffer));
   	}

   	cout << "To analyze: " << inFiles.size() << " files" << endl;
   	
	bDist = new TH1D("bDIST", "Impact Parameter distribution;Impact parameter b;Entries", 1000, 0, 20);

	worker();

	TCanvas *canv = new TCanvas("bdist","b dist", 800, 600);
	bDist->Draw();
	canv->SaveAs("epos_b_dist.png");
	canv->SaveAs("epos_b_dist.root");
	delete canv;
	return 0;
}

void* worker()
{
	int iEvent, totalEventCount;
	Float_t bim;
	TChain *tEPOS;

	cout << "Worker has started" << endl;

	tEPOS = new TChain("teposevent");

	while(!inFiles.empty())
    {
    	// cout << "Loading " << inFiles.back() << endl;
    	tEPOS->Add(inFiles.back().c_str());
    	inFiles.pop_back();
    }

	tEPOS->SetBranchStatus("*",1);
	tEPOS->SetBranchAddress("bim", &bim);

	totalEventCount = tEPOS->GetEntries();

	for(iEvent = 0; iEvent < totalEventCount; ++iEvent)
	{
		tEPOS->GetEntry(iEvent);
		bDist->Fill(bim);
		if(iEvent % 200 == 0)
			cout << "\t" << iEvent << " /\t" << totalEventCount << endl;
	}
	// delete tEPOS;
	return 0;
}