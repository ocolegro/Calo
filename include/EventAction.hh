#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TTree.h"
#include "SamplingSection.hh"
#include "G4SiHit.hh"
#include "G4Track.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSGeometryConversion.hh"

#include <vector>
#include <map>
#include "fstream"

class EventActionMessenger;
class StackingAction;
class SamplingSection;

class EventAction: public G4UserEventAction {
public:
	EventAction(G4bool doFast);
	virtual ~EventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void Detect(G4double eRawDep,G4VPhysicalVolume *volume,G4Track* lTrack,const G4ThreeVector & position);

	void SetPrintModulo(G4int val) {
		printModulo = val;
	}
	;
	void Add(std::vector<SamplingSection> *newDetector) {
		detector_ = newDetector;
	}
	inline void SetWait(G4bool wait) {
		wait_ = wait;
	}
	;
	inline G4bool GetWait() {
		return wait_;
	}
	;
	inline G4bool firstPass() {
		return doFast_;
	}
	;
	template <typename T>
	void FreeAll( T & t )
	;
	HGCSSGenParticleVec genvec_,hadvec_,incvec_,escapevec_,novelVec_;
	G4int evtNb_;
	G4int hadronicInts;
	std::vector<G4double> targetPartEngs;
	std::vector<G4int> novelPartEngs;
private:
	HGCSSSamplingSectionVec ssvec_;
	std::vector<SamplingSection> *detector_;
	G4int printModulo,initLayer;
	G4bool wait_;
	G4bool doFast_;
	HGCSSSimHitVec hitvec_;
	HGCSSGeometryConversion* geomConv_;
	TFile *outF_;
	TTree *tree_;
	HGCSSEvent event_;

	EventActionMessenger* eventMessenger;

	//std::ofstream fout_;

};

#endif

