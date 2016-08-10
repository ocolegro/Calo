#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "HGCSSInfo.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>
#include "TVector3.h"
#include "TMath.h"

#include "EventActionMessenger.hh"
#include "StackingAction.hh"
#include "SamplingSection.hh"

//
EventAction::EventAction(G4bool doFast) {
	eventMessenger = new EventActionMessenger(this);
	printModulo = 100;
	doFast_ = doFast;
	outF_ = TFile::Open("PFcal.root", "RECREATE");

	outF_->cd();
	double xysize =
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCalorSizeXY();
	initLayer = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer();

	//save some info
	HGCSSInfo *info = new HGCSSInfo();
	info->calorSizeXY(xysize);
	info->cellSize(CELL_SIZE_X);
	info->model(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getModel());
	info->version(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getVersion());

	std::cout << " -- check Info: version = " << info->version() << " model = "
			<< info->model() << std::endl;
	outF_->WriteObjectAny(info, "HGCSSInfo", "Info");

	//honeycomb
	geomConv_ = new HGCSSGeometryConversion(info->model(), CELL_SIZE_X);
	geomConv_->initialiseHoneyComb(xysize, CELL_SIZE_X);
	//square map for BHCAL
	geomConv_->initialiseSquareMap(xysize, 10.);

	tree_ = new TTree("HGCSSTree", "HGC Standalone simulation tree");
	tree_->Branch("HGCSSEvent", "HGCSSEvent", &event_);
	tree_->Branch("HGCSSSimHitVec", "std::vector<HGCSSSimHit>", &hitvec_);
	tree_->Branch("HGCSSGenAction", "std::vector<HGCSSGenParticle>",
			&genvec_);
	tree_->Branch("HGCSSHadAction", "std::vector<HGCSSGenParticle>",
			&hadvec_);
	tree_->Branch("HGCSSIncAction", "std::vector<HGCSSGenParticle>",
			&incvec_);
	tree_->Branch("HGCSSEscapeAction", "std::vector<HGCSSGenParticle>",
			&escapevec_);
	tree_->Branch("HGCSSNovelAction", "std::vector<HGCSSGenParticle>",
			&novelVec_);
}

//
EventAction::~EventAction() {
	outF_->cd();
	tree_->Write();
	outF_->Close();
	//fout_.close();
	delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt) {
	evtNb_ = evt->GetEventID();
	if (evtNb_ % printModulo == 0) {
		G4cout << "\n---> Begin of event: " << evtNb_ << G4endl;
		CLHEP::HepRandom::showEngineStatus();
	}

}

//
void EventAction::Detect(G4double eDepRaw, G4VPhysicalVolume *volume,G4Track* lTrack,const G4ThreeVector & position) {
	std::pair<G4bool,G4bool> stopIter = std::make_pair(false,false);
	for (size_t i = initLayer; i < detector_->size(); i++){
		if (stopIter.first) break;
		stopIter = (*detector_)[i].add( eDepRaw, volume,lTrack, position);
	}

}


void EventAction::EndOfEventAction(const G4Event* g4evt) {
	//return;
	event_.eventNumber(evtNb_);
	double totalSens = 0;
	double wgtTotalSens = 0;

	//if (summedDep < depCut){
	G4String fileN = "currentEvent.rndm";
	CLHEP::HepRandom::saveEngineStatus(fileN);
	std::ifstream input(fileN);
	std::string currentLine;
	Double_t stat_x = 0,stat_y = 0,seed_x = 0,seed_y = 0;
	for (size_t i = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer()
			; i < detector_->size(); i++) {
		Double_t weight = (i < 8) ? .8 : 1. ;
		totalSens += (*detector_)[i].getTotalSensE();
		wgtTotalSens += weight*(*detector_)[i].getTotalSensE();

		} //loop on sensitive layers
	for(int count = 0; count < 5; count++ ){
		getline( input, currentLine );
		if (count == 1)
			stat_x = std::atoi(currentLine.c_str());
		if (count == 2)
			stat_y = std::atoi(currentLine.c_str());

		if (count == 3)
			seed_x = std::atoi(currentLine.c_str());
		if (count == 4)
			seed_y = std::atoi(currentLine.c_str());
		//}
		TVector3 status(stat_x,stat_y,0);
		TVector3 seeds(seed_x,seed_y,0);

		event_.seeds(seeds);
		event_.status(status);
		//Changing initLayer because initial layers contain tracking sections.
	}
	for (size_t i =  ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer()
			; i < detector_->size(); i++) {
		for (unsigned idx(0); idx < (*detector_)[i].n_sens_elements; ++idx) {
					std::map<unsigned, HGCSSSimHit> lHitMap;
					std::pair<std::map<unsigned, HGCSSSimHit>::iterator, bool> isInserted;

					for (unsigned iSiHit(0);iSiHit < (*detector_)[i].getSiHitVec(idx).size();++iSiHit) {
						G4SiHit lSiHit = (*detector_)[i].getSiHitVec(idx)[iSiHit];
						bool is_scint = (*detector_)[i].hasScintillator;
						HGCSSSimHit lHit(lSiHit,is_scint ?geomConv_->squareMap() : geomConv_->hexagonMap());

						isInserted = lHitMap.insert(std::pair<unsigned, HGCSSSimHit>(lHit.cellid(), lHit));
						if (!isInserted.second)
							isInserted.first->second.Add(lSiHit);
					}

				std::map<unsigned, HGCSSSimHit>::iterator lIter = lHitMap.begin();
				hitvec_.reserve(hitvec_.size() + lHitMap.size());
				for (; lIter != lHitMap.end(); ++lIter) {
					hitvec_.push_back(lIter->second);
				}
			}
		(*detector_)[i].resetCounters();
	}
	//std::cout << "The initial layer is i = " << ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer() << std::endl;
	//std::cout << "The total deposited energy is " << totalSens << std::endl;
	event_.dep(totalSens);
	event_.wgtDep(wgtTotalSens);
	SetWait(false);
	tree_->Fill();

	//reset vectors
	genvec_.clear();
	hadvec_.clear();
	incvec_.clear();
	escapevec_.clear();
	novelVec_.clear();
	hitvec_.clear();
	targetPartEngs.clear();
	novelPartEngs.clear();

	hadronicInts = 0;

}
