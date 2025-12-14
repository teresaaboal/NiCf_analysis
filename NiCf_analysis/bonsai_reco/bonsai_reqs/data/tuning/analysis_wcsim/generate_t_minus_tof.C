#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

#include "WCSimRootOptions.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootEvent.hh"

// Simple example of reading a generated Root file
int generate_t_minus_tof(const char *filename, bool verbose=false)
{
  // Clear global scope
  //gROOT->Reset();

  // Open the file
  TFile * file = new TFile(filename,"read");
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }

  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");

  // Get the number of events
  long nevent = tree->GetEntries();
  if(verbose) printf("nevent %ld\n",nevent);

  // Create a WCSimRootEvent to put stuff from the tree in
  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);

  // Force deletion to prevent memory leak
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0;
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
    exit(9);
  }
  geotree->GetEntry(0);

  // Options tree - only need 1 "event"
  TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
  WCSimRootOptions *opt = 0;
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;

  //output file
  TString foutname(filename);
  foutname.ReplaceAll(".root", "_tof.root");
  TFile fout(foutname, "RECREATE");

  //output histograms
  TH1D *hndigi = new TH1D("hndigi", "PMT Hits;N digitised hits;Frequency in bin", 8000, 0, 8000);
  TH1D *hvtxX = new TH1D("hvtxX", "Event VTX0;True vertex X / cm;Frequency in bin", 200, -1500, 1500);
  TH1D *hvtxY = new TH1D("hvtxY", "Event VTX1;True vertex Y / cm;Frequency in bin", 200, -1500, 1500);
  TH1D *hvtxZ = new TH1D("hvtxZ", "Event VTX2;True vertex Z / cm;Frequency in bin", 200, -1500, 1500);
  TH1D *httof = new TH1D("likelihood", "Digitised hit time - time of flight;Digitised hit time - time of flight / ns;Frequency in bin", 750, -50, 250);

  int num_trig=0;
  const double c = 21.58333; //cm/ns. Group speed in pure water.

  // Now loop over events
  for(long iev = 0; iev < nevent; iev++) {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(iev);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %ld\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());

      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    }
    hvtxX->Fill(wcsimrootevent->GetVtx(0));
    hvtxY->Fill(wcsimrootevent->GetVtx(1));
    hvtxZ->Fill(wcsimrootevent->GetVtx(2));

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();

    hndigi->Fill(ncherenkovdigihits);
    if(verbose){
      printf("Event number: %ld\n", iev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
    }

    // Look at digitized hit info

    // Get the number of digitized hits
    // Loop over sub events

    if(verbose) cout << "DIGITIZED HITS:" << endl;
    const int ntriggers = wcsimrootsuperevent->GetNumberOfEvents();
    if(ntriggers > 1)
      cerr << "WARN: There are " << ntriggers << " in this event. Just using first for t-tof histogram" << endl;
    for(int index = 0 ; index < 1; index++) {
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
      if(verbose) cout << "Sub event number = " << index << "\n";

      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
      if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      int ncherenkovdigihits_slots = wcsimrootevent->GetNcherenkovdigihits_slots();

      if(ncherenkovdigihits>0)
	num_trig++;
      int idigi = 0;
      for(int i=0; i < ncherenkovdigihits_slots; i++) {
	// Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
	TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	if(!element) continue;
	idigi++;

	WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit =
	  dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

	double time = wcsimrootcherenkovdigihit->GetT();
	const WCSimRootPMT * pmt = geo->GetPMTPtr(wcsimrootcherenkovdigihit->GetTubeId() - 1); //Tube ID runs from 1 to N. This method directly accesses TClonesArray, with contents running from 0 to N-1
	// distance_of_flight in cm
	double dof = TMath::Sqrt(
	  TMath::Power(pmt->GetPosition(0) - wcsimrootevent->GetVtx(0), 2) +
	  TMath::Power(pmt->GetPosition(1) - wcsimrootevent->GetVtx(1), 2) + 
	  TMath::Power(pmt->GetPosition(2) - wcsimrootevent->GetVtx(2), 2));
	double tof = dof / c;
	httof->Fill(time - tof);
	if(verbose){
	  if(i < 10) // Only print first XX=10 tubes
	    printf("q, t, tubeid: %f %f %d \n\tFlight distance (cm), tof (ns): %f %f\n",
		   wcsimrootcherenkovdigihit->GetQ(),
		   wcsimrootcherenkovdigihit->GetT(),
		   wcsimrootcherenkovdigihit->GetTubeId(),
		   dof, tof);
	}
      } // End of loop over Cherenkov digihits
      if(verbose)
	cout << idigi << " digits found; expected " << ncherenkovdigihits << endl;
    } // End of loop over trigger

    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();

  } // End of loop over events

  hvtxX->Write();
  hvtxY->Write();
  hvtxZ->Write();
  hndigi->Write();
  httof->Write();

  TCanvas can;
  can.SetLogy(1);
  httof->Draw();
  TString plotname(filename);
  plotname.ReplaceAll(".root", "_tof.pdf");
  can.SaveAs(plotname);

  std::cout<<"num_trig "<<num_trig<<"\n";

  return 0;
}
