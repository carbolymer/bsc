#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <TH1D.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TThread.h>

using namespace std;

const int NTHREADS = 1; // number of simultaneously running workers

const int maxParticleCount = 1e5;
vector<string> inFiles;
TH1D *bDist;
TThread *threads[NTHREADS];

void* worker(void*);

int EPOSimpactDistribution()
{
#ifdef __CINT__
   printf("This script can only be executed via ACliC: .x %s++\n",__FILE__);
   return;
#endif
	char buffer[256];
	int toLoad[NTHREADS];
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
   	
   	for(int i = 0; i < NTHREADS; ++i)
   	{
   		number << i;
		toLoad[i] = (int) ((double)inFiles.size()/NTHREADS);
   		if(i == NTHREADS-1)
   			toLoad[i] += inFiles.size()%NTHREADS;
   		cout << toLoad[i] << endl;
   		threads[i] = new TThread((string("thread_")+number.str()).c_str(), worker, (void*) (toLoad+i) );
   	}


	bDist = new TH1D("bDIST", "Impact Parameter distribution;Impact parameter b;Entries", 80, 0, 20);

	// launch worker
	cout << "launching workers" << endl;
   	for(int i = 0; i < NTHREADS; ++i)
   	{
   		threads[i]->Run();
   	}

   	TThread::Ps();

	cout << "waiting to join" << endl;
   	for(int i = 0; i < NTHREADS; ++i)
   	{
   		threads[i]->Join();
   	}

	TCanvas *canv = new TCanvas("bdist","b dist", 800, 600);
	bDist->Draw();
	canv->SaveAs("epos_b_dist.png");
	// canv->SaveAs("epos_b_dist.root");
	delete canv;
	return 0;
}

void* worker(void* numberOfFiles)
{
	int iEvent, totalEventCount, nFiles;
	Float_t bim;
	TChain *tEPOS;

	cout << "Worker has started" << endl;

	nFiles = *(int*) numberOfFiles;
	cout << "# Loading " << nFiles << " files." << endl;

	tEPOS = new TChain("teposevent");

	TThread::Lock();
	while(!inFiles.empty() && nFiles-- > 0)
    {
    	cout << "Loading " << inFiles.back() << endl;
    	tEPOS->Add(inFiles.back().c_str());
    	inFiles.pop_back();
    }
	TThread::UnLock();

	tEPOS->SetBranchStatus("*",1);
	tEPOS->SetBranchAddress("bim", &bim);

	totalEventCount = tEPOS->GetEntries();

	for(iEvent = 0; iEvent < totalEventCount; ++iEvent)
	{
		tEPOS->GetEntry(iEvent);
		if(iEvent % 200 == 0)
			cout << "\t" << iEvent << " /\t" << totalEventCount << endl;
		TThread::Lock();
		// cout << bim << endl;
		bDist->Fill(bim);
		TThread::UnLock();
	}
	delete tEPOS;
	return 0;
}