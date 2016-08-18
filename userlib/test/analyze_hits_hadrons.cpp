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

bool checkDuplicate(std::vector<double> oldEng,double newEng){
	for (unsigned i = 0; i < oldEng.size(); i++){
		if ( (oldEng.at(i) - newEng) < .01) return false;
	}
	return true;
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

	std::vector<HGCSSGenParticle> * incVec = 0;
	tree->SetBranchAddress("HGCSSIncAction", &incVec);

	std::vector<HGCSSGenParticle> * genVec = 0;
	tree->SetBranchAddress("HGCSSGenAction", &genVec);

	std::vector<HGCSSGenParticle> * hadVec = 0;
	tree->SetBranchAddress("HGCSSHadAction", &hadVec);

	std::vector<HGCSSGenParticle> * escapeVec = 0;
	tree->SetBranchAddress("HGCSSEscapeAction", &escapeVec);

	std::vector<HGCSSGenParticle> * novelVec = 0;
	tree->SetBranchAddress("HGCSSNovelAction", &novelVec);


	std::vector<HGCSSSimHit> * hitVec_ = 0;
	tree->SetBranchAddress("HGCSSSimHitVec", &hitVec_);
	
	TFile hfile("analyzed_tuple.root", "RECREATE");
	TTree t1("sampling", "Hadronic Study");

	Int_t nInteractions,nSecondaries[50000],nProtonSecondaries[50000],nNeutronSecondaries[50000],
	nOtherSecondaries[50000],nContainedSecondaries[50000],nUncontainedSecondaries[50000],seeds[4];

	t1.Branch("nInteractions", &nInteractions, "nInteractions/I");
	t1.Branch("nSecondaries", &nSecondaries, "nSecondaries[nInteractions]/I");
	t1.Branch("nProtonSecondaries", &nProtonSecondaries, "nProtonSecondaries[nInteractions]/I");
	t1.Branch("nNeutronSecondaries", &nNeutronSecondaries, "nNeutronSecondaries[nInteractions]/I");
	t1.Branch("nOtherSecondaries", &nOtherSecondaries, "nOtherSecondaries[nInteractions]/I");
	t1.Branch("nContainedSecondaries", &nContainedSecondaries, "nContainedSecondaries[nInteractions]/I");
	t1.Branch("nUncontainedSecondaries", &nUncontainedSecondaries, "nUncontainedSecondaries[nInteractions]/I");

	t1.Branch("seeds", &seeds, "seeds[4]/I");

	Float_t inc_KE[50000],inc_zpos[50000],inc_theta[50000];
	Int_t   inc_pdgid[50000];

	t1.Branch("inc_pdgid", &inc_pdgid, "inc_pdgid[nInteractions]/I");
	t1.Branch("inc_KE", &inc_KE, "inc_KE[nInteractions]/F");
	t1.Branch("inc_zpos", &inc_zpos, "inc_zpos[nInteractions]/F");
	t1.Branch("inc_theta", &inc_theta, "inc_theta[nInteractions]/F");

	Float_t out_KE[50000],out_NE[50000],out_PE[50000],out_Eff[50000],out_OE[50000];

	t1.Branch("out_KE", &out_KE, "out_KE[nInteractions]/F");
	t1.Branch("out_NE", &out_NE, "out_NE[nInteractions]/F");
	t1.Branch("out_PE", &out_PE, "out_PE[nInteractions]/F");
	t1.Branch("out_OE", &out_OE, "out_OE[nInteractions]/F");
	t1.Branch("out_Eff", &out_Eff, "out_Eff[nInteractions]/F");

	Float_t hadron_zpos[50000],
	hadron_theta[50000],hadron_px[50000]  , hadron_py[50000]  ,hadron_pz[50000],
	hadron_KE[50000];
	Int_t nHadrons,hadron_pdgid[50000];
    Int_t hadron_int[50000];

	t1.Branch("nHadrons", &nHadrons, "nHadrons/I");
	t1.Branch("hadron_pdgid", &hadron_pdgid, "hadron_pdgid[nHadrons]/I");
	t1.Branch("hadron_zpos", &hadron_zpos, "hadron_zpos[nHadrons]/F");
	t1.Branch("hadron_px", &hadron_px, "hadron_px[nHadrons]/F");
	t1.Branch("hadron_py", &hadron_py, "hadron_py[nHadrons]/F");
	t1.Branch("hadron_pz", &hadron_pz, "hadron_pz[nHadrons]/F");
	t1.Branch("hadron_theta", &hadron_theta, "hadron_theta[nHadrons]/F");
	t1.Branch("hadron_KE", &hadron_KE, "hadron_KE[nHadrons]/F");
    t1.Branch("hadron_int", &hadron_int, "hadron_int[nHadrons]/I");



	Float_t novel_zpos[50000],
	novel_theta[50000],novel_px[50000]  , novel_py[50000]  ,novel_pz[50000],
	novel_KE[50000],novel_parentKE[50000];
	Int_t nNovels,novel_pdgid[50000],novel_parentID[50000];

	t1.Branch("nNovels", &nNovels, "nNovels/I");
	t1.Branch("novel_pdgid", &novel_pdgid, "novel_pdgid[nNovels]/I");
	t1.Branch("novel_zpos", &novel_zpos, "novel_zpos[nNovels]/F");
	t1.Branch("novel_px", &novel_px, "novel_px[nNovels]/F");
	t1.Branch("novel_py", &novel_py, "novel_py[nNovels]/F");
	t1.Branch("novel_pz", &novel_pz, "novel_pz[nNovels]/F");
	t1.Branch("novel_theta", &novel_theta, "novel_theta[nNovels]/F");
	t1.Branch("novel_KE", &novel_KE, "novel_KE[nNovels]/F");
	t1.Branch("novel_parentID", &novel_parentID, "novel_parentID[nNovels]/I");
	t1.Branch("novel_parentKE", &novel_parentKE, "novel_parentKE[nNovels]/F");

	Float_t escape_zpos[50000],escape_xpos[50000],escape_ypos[50000],
	escape_theta[50000],escape_px[50000]  , escape_py[50000]  ,escape_pz[50000],
	escape_VKE[50000],escape_FKE[50000];
	Int_t nEscapes,escape_pdgid[50000];
	t1.Branch("nEscapes", &nEscapes, "nEscapes/I");
	t1.Branch("escape_pdgid", &escape_pdgid, "escape_pdgid[nEscapes]/I");
	t1.Branch("escape_theta", &escape_theta, "escape_theta[nEscapes]/F");
	t1.Branch("escape_px", &escape_px, "escape_px[nEscapes]/F");
	t1.Branch("escape_py", &escape_py, "escape_py[nEscapes]/F");
	t1.Branch("escape_pz", &escape_pz, "escape_pz[nEscapes]/F");
	t1.Branch("escape_xpos", &escape_xpos, "escape_xpos[nEscapes]/F");
	t1.Branch("escape_ypos", &escape_ypos, "escape_ypos[nEscapes]/F");
	t1.Branch("escape_zpos", &escape_zpos, "escape_zpos[nEscapes]/F");
	t1.Branch("escape_VKE", &escape_VKE, "escape_VKE[nEscapes]/F");
	t1.Branch("escape_FKE", &escape_FKE, "escape_FKE[nEscapes]/F");
	Float_t summedSen,summedSenWgt,convEng_1,accconvEng_1,lostEng_1,convEng_2,accconvEng_2,lostEng_2,initEng,initX,initY,initZ;

	t1.Branch("initEng", &initEng, "initEng/F");
	t1.Branch("initX", &initX, "initX/F");
	t1.Branch("initY", &initY, "initY/F");
	t1.Branch("initZ", &initZ, "initZ/F");

	t1.Branch("convEng_1", &convEng_1, "convEng_1/F");
	t1.Branch("accconvEng_1", &accconvEng_1, "accconvEng_1/F");
	t1.Branch("lostEng_1", &lostEng_1, "lostEng_1/F");

	t1.Branch("convEng_2", &convEng_2, "convEng_2/F");
	t1.Branch("accconvEng_2", &accconvEng_2, "accconvEng_2/F");
	t1.Branch("lostEng_2", &lostEng_2, "lostEng_2/F");


	t1.Branch("summedSen", &summedSen, "summedSen/F");
	t1.Branch("summedSenWgt", &summedSenWgt, "summedSenWgt/F");
	TH2PolyBin *centerCell = 0;
	TH2PolyBin *neighborCell = 0;
	unsigned nHits = 0,cellID[50000],cellLayer[50000];
	Float_t cellEnergy[50000],cellParentID[50000],cellParentKE[50000],cellParentTrack[50000],cellRellIso[50000],engDep;
	unsigned nHadronsHit,nGammas,nElectrons,nProtons,nMuons;
	t1.Branch("nHits", &nHits, "nHits/I");
	t1.Branch("cellID", &cellID, "cellID[nHits]/I");
	t1.Branch("cellLayer", &cellLayer, "cellLayer[nHits]/I");
	t1.Branch("cellEnergy", &cellEnergy, "cellEnergy[nHits]/F");
	t1.Branch("cellParentID", &cellParentID, "cellParentID[nHits]/F");
	t1.Branch("cellParentKE", &cellParentKE, "cellParentKE[nHits]/F");
	t1.Branch("cellParentTrack", &cellParentTrack, "cellParentTrack[nHits]/F");
	t1.Branch("cellRellIso", &cellRellIso, "cellRellIso[nHits]/F");

	t1.Branch("nHadronsHit", &nHadronsHit, "nHadronsHit/I");
	t1.Branch("nGammas", &nGammas, "nGammas/I");
	t1.Branch("nElectrons", &nElectrons, "nElectrons/I");
	t1.Branch("nProtons", &nProtons, "nProtons/I");
	t1.Branch("nMuons", &nMuons, "nMuons/I");
	


	std::vector<double> hadronKEs,novelKEs;

	unsigned nEvts = tree->GetEntries();
	for (unsigned ievt(0); ievt < nEvts; ++ievt) { //loop on entries
		//std::cout << "The event is ievt = " << ievt << std::endl;
		tree->GetEntry(ievt);
		summedSen = evt_->dep();
		summedSenWgt = evt_->wgtDep();
		if (summedSen == 0) continue;
		seeds[0] = evt_->status().x();
		seeds[1] = evt_->status().y();
		seeds[2] = evt_->seeds().x();
		seeds[3] = evt_->seeds().y();
		lostEng_1 = 0;
		convEng_1 = 0;
		accconvEng_1 = 0;

		lostEng_2 = 0;
		convEng_2 = 0;
		accconvEng_2 = 0;


		nInteractions = 0;
		nHadrons = 0;
		nEscapes = 0;
		nNovels = 0;

		nHits = hitVec_->size();

		for (Int_t j = 0; j < incVec->size(); j++) {

			HGCSSGenParticle& incPart = (*incVec)[j];
			TVector3 momVec = incPart.vertexMom();
			TVector3 posVec = incPart.vertexPos();
			unsigned iLoc		=	-incPart.layer() - 1;
			inc_KE[iLoc] = incPart.vertexKE();
			inc_zpos[iLoc] = posVec[2];
			inc_theta[iLoc] = acos(momVec[2])*180/3.14;
			inc_pdgid[iLoc] = incPart.pdgid();

			nInteractions = iLoc + 1;
			nSecondaries[iLoc] = 0;
			nNeutronSecondaries[iLoc] = 0;
			nProtonSecondaries[iLoc] = 0;
			nOtherSecondaries[iLoc] = 0;
			nContainedSecondaries[iLoc] = 0;
			nUncontainedSecondaries[iLoc] = 0;
			out_PE[iLoc] = 0;
			out_NE[iLoc] = 0;
			out_OE[iLoc] = 0;
			out_Eff[iLoc] = 0;
			out_KE[iLoc] = 0;

		}

		for (Int_t j = 0; j < novelVec->size(); j++) {
			HGCSSGenParticle& novel = (*novelVec)[j];

			if (checkDuplicate(novelKEs,novel.vertexKE()) == false) continue;
			novelKEs.push_back(novel.vertexKE());
			nNovels = nNovels + 1;
			TVector3 momVec = novel.vertexMom();
			TVector3 posVec = novel.vertexPos();

			novel_zpos[j]   	= posVec[2];
			novel_px[j]   		= momVec[0];
			novel_py[j]   		= momVec[1];
			novel_pz[j]   		= momVec[2];
			novel_theta[j]   	= acos(momVec[2]) * 180/3.14;
			novel_pdgid[j]   	= novel.pdgid();
			novel_KE[j]			= novel.vertexKE();
			novel_parentID[j]			= novel.parentPdgId();
			novel_parentKE[j]			= novel.parentKE();

		}

		for (Int_t j = 0; j < hadVec->size(); j++) {
			HGCSSGenParticle& hadron = (*hadVec)[j];

			if (checkDuplicate(hadronKEs,hadron.vertexKE()) == false) continue;
			hadronKEs.push_back(hadron.vertexKE());
			nHadrons = nHadrons + 1;
			TVector3 momVec = hadron.vertexMom();
			TVector3 posVec = hadron.vertexPos();

			hadron_zpos[j]   	= posVec[2];
			hadron_px[j]   		= momVec[0];
			hadron_py[j]   		= momVec[1];
			hadron_pz[j]   		= momVec[2];
			hadron_theta[j]   	= acos(momVec[2]) * 180/3.14;
			hadron_pdgid[j]   	= hadron.pdgid();
			hadron_KE[j]		= hadron.vertexKE();
			unsigned iLoc		=	hadron.layer() - 1;
            if (iLoc > nInteractions) nInteractions = iLoc +1;

			out_KE[iLoc] += hadron_KE[j];
            hadron_int[j] = iLoc;
            /* Don't count photons or anti-electrons*/
            if((abs(hadron_pdgid[j]) != 22) && (abs(hadron_pdgid[j]) != 11)){
                nSecondaries[iLoc] += 1;
            }
			bool acc = false;
			if (hadron_theta[j] < 30){
				nContainedSecondaries[iLoc] += 1;
				acc = true;
			}
			else
				nUncontainedSecondaries[iLoc] += 1;

			//Do hadron stuff
			if (abs(hadron_pdgid[j]) == 2112 || abs(hadron_pdgid[j])  == 2212){
				convEng_1 += hadron_KE[j];
				convEng_2 += hadron_KE[j];
				if (acc) {
					accconvEng_1 += hadron_KE[j];
					accconvEng_2 += hadron_KE[j];
				}
				out_Eff[iLoc] += hadron_KE[j];


				// Do neutron stuff
				if (abs(hadron_pdgid[j]) == 2112){
					nNeutronSecondaries[iLoc] = nNeutronSecondaries[iLoc] + 1;
					if (hadron_pdgid[j] > 0)
						out_NE[iLoc] += hadron_KE[j];
					else{
						out_NE[iLoc] += hadron_KE[j] +  2*hadron.mass() ;
						convEng_1 +=  2*hadron.mass();

						out_Eff[iLoc]+=  2*hadron.mass();

						if (acc) accconvEng_1 += 2*hadron.mass();

					}

				}
				//Do Proton stuff
				if (abs(hadron_pdgid[j]) == 2212){
					nProtonSecondaries[iLoc] += 1;
					if (hadron_pdgid[j] > 0)
						out_PE[iLoc] += hadron_KE[j];
					else{
						out_PE[iLoc] += hadron_KE[j] +  hadron.mass() ;
						convEng_1 +=  2*hadron.mass();

						out_Eff[iLoc] +=  2*hadron.mass();

						if (acc) accconvEng_1 += 2*hadron.mass();
					}

				}

			}
			else if (abs(hadron_pdgid[j]) <10000  && abs(hadron_pdgid[j]) != 22 && abs(hadron_pdgid[j]) != 11
					&& hadron_pdgid[j] != 0){
				convEng_1 += hadron_KE[j] + hadron.mass();
				convEng_2 += hadron_KE[j];

				if (acc) {
					accconvEng_1 += hadron_KE[j]+hadron.mass();
					accconvEng_2 += hadron_KE[j];
				}

				out_OE[iLoc] += hadron_KE[j]+hadron.mass();
				out_Eff[iLoc] +=  hadron.mass();
				nOtherSecondaries[iLoc] += 1;

			}
		}

		for (Int_t j = 0; j < escapeVec->size(); j++) {

			HGCSSGenParticle& escape = (*escapeVec)[j];
			nEscapes = nEscapes + 1;
			TVector3 momVec = escape.vertexMom();
			TVector3 posVec = escape.vertexPos();
			escape_xpos[j]   	= posVec[0];
			escape_ypos[j]   	= posVec[1];
			escape_zpos[j]   	= posVec[2];
			escape_px[j]   		= momVec[0];
			escape_py[j]   		= momVec[1];
			escape_pz[j]   		= momVec[2];
			escape_theta[j]   	= acos(momVec[2]) * 180/3.14;

			escape_pdgid[j]   	= escape.pdgid();
			escape_VKE[j]		= escape.vertexKE();
			escape_FKE[j]			= escape.finalKE();

			if (escape_pdgid[j] == -2112 || escape_pdgid[j]  == -2212){
				lostEng_1 += 2* escape.mass();
				lostEng_1 += escape_FKE[j];
				lostEng_2 += escape_FKE[j];

			}
			else if (escape_pdgid[j] == 2112 || escape_pdgid[j]  == 2212){
				lostEng_1 += escape_FKE[j];
				lostEng_2 += escape_FKE[j];

			}
			else if (abs(escape_pdgid[j]) < 10000 && escape_pdgid[j] != 0){
				lostEng_1 += escape_FKE[j]+escape.mass();
				lostEng_2 += escape_FKE[j];

			}
		}
		initEng = genVec->at(0).vertexKE();
		TVector3 initPos = genVec->at(0).vertexPos();

		initX = initPos[0];
		initY = initPos[1];
		initZ = initPos[2];
		for (unsigned j = 0; j < hitVec_->size(); j++) {
			HGCSSSimHit& hit = (*hitVec_)[j];
			cellLayer[j] 		= hit.layer_;
			cellID[j]			= hit.cellid_;
			cellEnergy[j]		= hit.energy_;
			cellParentID[j]		= hit.pdgIDMainParent_;
			cellParentKE[j]		= hit.KEMainParent_;
			cellParentTrack[j]	= hit.trackIDMainParent_;
			nHadronsHit 			= hit.nHadrons_;
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
		hadronKEs.clear();novelKEs.clear();
	}
	t1.Write();

	return 1;
}
