#include <iostream>
#include <string>
#include <map>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include "/home/mgalazyn/workspace/Therminator2/build/include/StructEvent.h"
#include "/home/mgalazyn/workspace/Therminator2/build/include/ParticleCoor.h"

const unsigned int maxParticleCount = 1e5;

// tree->Draw("sqrt(px*px+py*py)","pid=211 && t<100")
// "x*x+y*y"
// 211 - pi
// 321 - K
// 2212 - p

typedef struct
{
	Int_t		npart;
	Int_t		id[maxParticleCount];
	Int_t		mid[maxParticleCount];
	Float_t		x[maxParticleCount];
	Float_t		y[maxParticleCount];
	Float_t		z[maxParticleCount];
	Float_t		t[maxParticleCount];
	Float_t		px[maxParticleCount];
	Float_t		py[maxParticleCount];
	Float_t		pz[maxParticleCount];
	Float_t		E[maxParticleCount];
}Event;

void setBranchesAddresses(TChain* chain, Event &event);

int main(int argc, char **argv)
{
	int fileNumber = 1;
	char buffer[256];
	std::vector<std::string> inputFiles;
	if(argc < 1)
	{
		std::cout << "You have to specify output path in argument. Paths to input files should be given to stdin." << std::endl;
		return 1;
	}

	while (std::cin.good())
    {
    	std::fill(buffer, buffer+255, 0);
    	std::cin >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	inputFiles.push_back(buffer);
    }
    if(inputFiles.size() == 0)
    {
    	std::cout << "No input files were provided" << std::cin;
    	return 1;
    }

    // removing "/" at the end of path
    int i = 0;
   	while(argv[1][i++] != '\0');
   	if(argv[1][i-2] == '/')
   		argv[1][i-2] = '\0';


	int totalEventCount, previousCount = 0;
	Event eInitial, eFinal;
	StructEvent event;
	ParticleCoor particleCoordinates;

	TFile *outFile;
	TTree *particleTree;
	TTree *eventTree;
	TChain *chInitial = new TChain("ti");
	TChain *chFinal = new TChain("tf");

	// reading
	for(unsigned int i = 0; i < inputFiles.size(); ++i)
	{
		if( chInitial->Add(inputFiles[i].c_str()) == 0 
			|| chFinal->Add(inputFiles[i].c_str()) == 0 )
			std::cout << ">>> Could not open file: " << inputFiles[i] << std::endl;
		else
			std::cout << "Adding file:  " << inputFiles[i] << std::endl;
			
	}
	
	chInitial->SetBranchStatus("*",1);
	chFinal->SetBranchStatus("*",1);

	setBranchesAddresses(chInitial, eInitial);
	setBranchesAddresses(chFinal, eFinal);
	totalEventCount = chInitial->GetEntries();

	for(unsigned int iEvent = 0; iEvent < totalEventCount; ++iEvent)
	{
		if(iEvent % 50 == 0)
		{
			if(iEvent != 0)
			{
				outFile->Write();
				outFile->Close();
				delete outFile;
				delete particleTree;
				delete eventTree;
			}
			particleTree = new TTree(_PARTICLES_TREE_, "particle tree");    
			eventTree = new TTree(_EVENTS_TREE_,    "event tree");

			outFile = new TFile(TString::Format("%s/event%d.root",argv[1],fileNumber++),"RECREATE");
			if(outFile->IsZombie())
			{
				std::cout << "Output file opening failed." << std::endl;
				return 1;
			}

			outFile->cd();

		    particleTree->Branch(_PARTICLE_BRANCH_, &particleCoordinates, _PARTICLE_FORMAT_ );
		    eventTree->Branch(_EVENTS_BRANCH_, &event, _EVENTS_FORMAT_ );
			std::cout << "Writing to: " << TString::Format("%s/event%d.root",argv[1],fileNumber-1) << std::endl;
		}
		chInitial->GetEntry(iEvent);
		chFinal->GetEntry(iEvent);
		int eid = 0;

		for (int iParticle = 0; iParticle < eInitial.npart; ++iParticle)
		{
			particleCoordinates.mass = TMath::Sqrt(
				TMath::Power(eInitial.E[iParticle]/TMath::C(),2)
				- TMath::Power(eInitial.px[iParticle],2)
				- TMath::Power(eInitial.py[iParticle],2)
				- TMath::Power(eInitial.pz[iParticle],2)
			);
			particleCoordinates.t = eInitial.t[iParticle];
			particleCoordinates.x = eInitial.x[iParticle];
			particleCoordinates.y = eInitial.y[iParticle];
			particleCoordinates.z = eInitial.z[iParticle];
			particleCoordinates.e = eInitial.E[iParticle];
			particleCoordinates.px = eInitial.px[iParticle];
			particleCoordinates.py = eInitial.py[iParticle];
			particleCoordinates.pz = eInitial.pz[iParticle];
			particleCoordinates.decayed = 0;
			particleCoordinates.pid = eInitial.id[iParticle];
			particleCoordinates.fatherpid = 0;
			particleCoordinates.rootpid = 0;
			particleCoordinates.eid = ++eid;
			particleCoordinates.fathereid = 0;
			particleCoordinates.eventid = iEvent;

			particleTree->Fill();
		}

		for (int iParticle = 0; iParticle < eFinal.npart; ++iParticle)
		{
			particleCoordinates.mass = TMath::Sqrt(
				TMath::Power(eFinal.E[iParticle]*TMath::C(),2)
				- TMath::Power(eFinal.px[iParticle],2)
				- TMath::Power(eFinal.py[iParticle],2)
				- TMath::Power(eFinal.pz[iParticle],2)
			);
			particleCoordinates.t = eFinal.t[iParticle];
			particleCoordinates.x = eFinal.x[iParticle];
			particleCoordinates.y = eFinal.y[iParticle];
			particleCoordinates.z = eFinal.z[iParticle];
			particleCoordinates.e = eFinal.E[iParticle];
			particleCoordinates.px = eFinal.px[iParticle];
			particleCoordinates.py = eFinal.py[iParticle];
			particleCoordinates.pz = eFinal.pz[iParticle];
			particleCoordinates.decayed = 0;
			particleCoordinates.pid = eFinal.id[iParticle];
			particleCoordinates.fatherpid = eFinal.mid[iParticle];
			particleCoordinates.rootpid = 0;
			particleCoordinates.eid = ++eid;
			particleCoordinates.fathereid = 0;
			particleCoordinates.eventid = iEvent;

			particleTree->Fill();
		}
		
		event.eventID = 0;
		event.entries = eInitial.npart + eFinal.npart;
		event.entriesprev = previousCount;
		eventTree->Fill();
		previousCount = eInitial.npart + eFinal.npart;
	}

	outFile->Write();
	outFile->Close();
	delete outFile;

	return 0;
}

void setBranchesAddresses(TChain* chain, Event &event)
{
	chain->SetBranchAddress("npart",&event.npart);
	chain->SetBranchAddress("id", event.id);
	chain->SetBranchAddress("mid", event.mid );
	chain->SetBranchAddress("x", event.x);
	chain->SetBranchAddress("y", event.y);
	chain->SetBranchAddress("z", event.z);
	chain->SetBranchAddress("t", event.t);
	chain->SetBranchAddress("px", event.px);
	chain->SetBranchAddress("py", event.py);
	chain->SetBranchAddress("pz", event.pz);
	chain->SetBranchAddress("E", event.E);
}
