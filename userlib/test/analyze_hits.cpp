#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSSamplingSection.hh"
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

void myHoneycomb(TH2Poly* map, Double_t xstart,
		Double_t ystart, Double_t a,  // side length
		Int_t k,     // # hexagons in a column
		Int_t s)     // # columns
		{
	// Add the bins
	Double_t numberOfHexagonsInAColumn;
	Double_t x[6], y[6];
	Double_t xloop, yloop, ytemp;
	xloop = xstart;
	yloop = ystart + a * TMath::Sqrt(3) / 2.0;
	for (int sCounter = 0; sCounter < s; sCounter++) {

		ytemp = yloop; // Resets the temp variable

		// Determine the number of hexagons in that column
		if (sCounter % 2 == 0) {
			numberOfHexagonsInAColumn = k;
		} else {
			numberOfHexagonsInAColumn = k - 1;
		}

		for (int kCounter = 0; kCounter < numberOfHexagonsInAColumn;
				kCounter++) {

			// Go around the hexagon
			x[0] = xloop;
			y[0] = ytemp;
			x[1] = x[0] + a / 2.0;
			y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
			x[2] = x[1] + a;
			y[2] = y[1];
			x[3] = x[2] + a / 2.0;
			y[3] = y[1] - a * TMath::Sqrt(3) / 2.0;
			;
			x[4] = x[2];
			y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
			;
			x[5] = x[1];
			y[5] = y[4];

			map->AddBin(6, x, y);

			// Go up
			ytemp += a * TMath::Sqrt(3);
		}

		// Increment the starting position
		if (sCounter % 2 == 0)
			yloop += a * TMath::Sqrt(3) / 2.0;
		else
			yloop -= a * TMath::Sqrt(3) / 2.0;
		xloop += 1.5 * a;
	}
}

int main(int argc, char** argv) {
	TH2Poly *hcomb = new TH2Poly();
	double width = 622.5;
	double side = 4.5936;
	unsigned ncellwide = width / (2. * side);
	unsigned ny = ncellwide + 1;
	unsigned nx = ncellwide + 4;
	double xstart = -((double) ncellwide) * side ;
	double ystart = -((double) ncellwide) * side * sqrt(3) / 2.0;
	myHoneycomb(hcomb, xstart, ystart, side, ny, nx);
	TFile *infile = TFile::Open(argv[1]);
	TTree *tree = (TTree*) infile->Get("HGCSSTree");
	freopen("log.txt", "w", stdout);

	HGCSSEvent* evt_ = 0;
	tree->SetBranchAddress("HGCSSEvent", &evt_);

	std::vector<HGCSSSimHit> * hitVec_ = 0;
	tree->SetBranchAddress("HGCSSSimHitVec", &hitVec_);

	std::vector<HGCSSGenParticle> * genVec = 0;
	tree->SetBranchAddress("HGCSSGenAction", &genVec);

	TFile hfile("analyzed_tuple.root", "RECREATE");
	TTree t1("sampling", "Hadronic Study");


	unsigned nHits = 0,cellID[50000],cellLayer[50000];
	Float_t cellEnergy[50000],cellParentID[50000],cellParentKE[50000],cellParentTrack[50000],cellRellIso[50000],engDep;
	unsigned initEng,nHadrons,nGammas,nElectrons,nProtons,nMuons;
	t1.Branch("nHits", &nHits, "nHits/I");
	t1.Branch("cellID", &cellID, "cellID[nHits]/I");
	t1.Branch("cellLayer", &cellLayer, "cellLayer[nHits]/I");
	t1.Branch("cellEnergy", &cellEnergy, "cellEnergy[nHits]/F");
	t1.Branch("cellParentID", &cellParentID, "cellParentID[nHits]/F");
	t1.Branch("cellParentKE", &cellParentKE, "cellParentKE[nHits]/F");
	t1.Branch("cellParentTrack", &cellParentTrack, "cellParentTrack[nHits]/F");
	t1.Branch("cellRellIso", &cellRellIso, "cellRellIso[nHits]/F");

	t1.Branch("initEng", &initEng, "initEng/I");
	t1.Branch("engDep", &engDep, "engDep/F");
	t1.Branch("nHadrons", &nHadrons, "nHadrons/I");
	t1.Branch("nGammas", &nGammas, "nGammas/I");
	t1.Branch("nElectrons", &nElectrons, "nElectrons/I");
	t1.Branch("nProtons", &nProtons, "nProtons/I");
	t1.Branch("nMuons", &nMuons, "nMuons/I");
	TH2PolyBin *centerCell = 0;
	TH2PolyBin *neighborCell = 0;
	unsigned nEvts = tree->GetEntries();
	for (unsigned ievt(0); ievt < nEvts; ++ievt) { //loop on entries
		tree->GetEntry(ievt);
		nHits = hitVec_->size();
		//std::cout << "The size of hitVec_ is " << nHits << std::endl;
		initEng =  genVec->at(0).vertexKE();
		//if (initEng < 10) continue;
		engDep = evt_->wgtDep();
		/*if (engDep > 15.) {
			t1.Fill();
			continue;
		}*/
		for (unsigned j = 0; j < hitVec_->size(); j++) {
			HGCSSSimHit& hit = (*hitVec_)[j];
			cellLayer[j] 		= hit.layer_;
			cellID[j]			= hit.cellid_;
			cellEnergy[j]		= hit.energy_;
			cellParentID[j]		= hit.pdgIDMainParent_;
			cellParentKE[j]		= hit.KEMainParent_;
			cellParentTrack[j]	= hit.trackIDMainParent_;
			nHadrons 			= hit.nHadrons_;
			nGammas 			= hit.nGammas_;
			nElectrons 			= hit.nElectrons_;
			nProtons 			= hit.nProtons_;
			nMuons				= hit.nMuons_;
			if (cellEnergy[j] > .075){
				double outerDep = 0;
				for (unsigned k = 0; k < hitVec_->size(); k++) {
					HGCSSSimHit& nbr = (*hitVec_)[k];
					if (nbr.layer_ != cellLayer[j]) continue;
					if (nbr.cellid_ != cellID[j] +1 && nbr.cellid_ != cellID[j] - 1 &&
							nbr.cellid_ != cellID[j] +  67 && nbr.cellid_ != cellID[j] - 67 &&
							nbr.cellid_ != cellID[j] +  68 && nbr.cellid_ != cellID[j] - 67) continue;
					std::cout << "Getting the center bin " << std::endl;
					centerCell = (TH2PolyBin*) hcomb->GetBins()->At(cellID[j]-1);
					std::cout << "Getting the nbr bin " << std::endl;
					neighborCell = (TH2PolyBin*) hcomb->GetBins()->At(nbr.cellid_ - 1);
					if (centerCell != nullptr and neighborCell != nullptr){
						std::cout << "Computing the radius and summing the energy" << std::endl;
						double x_1 = (centerCell->GetXMax() + centerCell->GetXMin()) / 2.;
						double x_2 = (neighborCell->GetXMax() + neighborCell->GetXMin()) / 2.;

						double y_1 = (centerCell->GetYMax() + centerCell->GetYMin()) / 2.;
						double y_2 = (neighborCell->GetYMax() + neighborCell->GetYMin()) / 2.;

						if (pow( pow((x_1 - x_2),2) +pow((y_1 - y_2),2),.5) < 8)
							outerDep += nbr.energy_;
					}
				}
				if (outerDep > 0)
					cellRellIso[j] = cellEnergy[j]/(cellEnergy[j]+outerDep);
			}
		}
		t1.Fill();
	}
	t1.Write();

	return 1;
}
