#include <iostream>
#include <algorithm>
#include <TChain.h>

using namespace std;

const int maxParticleCount = 1e5;

int analyzeEPOS()
{
	int iEvent, totalEventCount, iParticle;
	Float_t bim;
	char buffer[256];

	TH1D *bDist = new TH1D("bDIST", "Impact Parameter distribution;Impact parameter b;Entries", 80, 0, 20);

	TChain *tEPOS = new TChain("teposevent");
	while (cin.good())
    {
    	std::fill(buffer, buffer+255, 0);
    	cin >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	cout << "Loading: " << buffer << endl;
    	tEPOS->Add(buffer);
   	}

	tEPOS->SetBranchStatus("*",1);
	tEPOS->SetBranchAddress("bim", &bim);
	// tEPOS->SetBranchAddress("zus", event.zus);
	// tEPOS->SetBranchAddress("px", event.px);
	// tEPOS->SetBranchAddress("py", event.py);
	// tEPOS->SetBranchAddress("pz", event.pz);
	// tEPOS->SetBranchAddress("e", event.e);
	// tEPOS->SetBranchAddress("x", event.x);
	// tEPOS->SetBranchAddress("y", event.y);
	// tEPOS->SetBranchAddress("z", event.z);
	// tEPOS->SetBranchAddress("t", event.t);
	// tEPOS->SetBranchAddress("id", event.id);
	// tEPOS->SetBranchAddress("ist", event.ist);
	// tEPOS->SetBranchAddress("ity", event.ity);
	// tEPOS->SetBranchAddress("ior", event.ior);
	// tEPOS->SetBranchAddress("jor", event.jor);

	totalEventCount = tEPOS->GetEntries();

	for(iEvent =0; iEvent < totalEventCount; ++iEvent)
	{
		tEPOS->GetEntry(iEvent);
		if(iEvent % 200 == 0)
			cout << "\t" << iEvent << " /\t" << totalEventCount << endl;
		bDist->Fill(bim);
	}
	TCanvas *canv = new TCanvas("bdist","b dist", 800, 600);
	bDist->Draw();
	canv->SaveAs("epos_b_dist.png");
	canv->SaveAs("epos_b_dist.root");
	delete canv;
	return 0;
}