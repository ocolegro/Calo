#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4EmSaturation.hh"

class EventAction;

class SteppingAction: public G4UserSteppingAction {
public:
	SteppingAction(std::string data = "");
	virtual ~SteppingAction();

	void UserSteppingAction(const G4Step*);

	inline bool checkDuplicate(std::vector<double> engVec,double currentEng){
		for (unsigned i = 0; i < engVec.size(); i++){
			if ( (engVec.at(i) - currentEng) < .01) return false;
		}
		return true;
	};
	G4double stepPDGID,stepKE;
private:
	void printParticle(G4Track* aTrack);
	EventAction *eventAction_;
	//to correct the energy in the scintillator
	G4EmSaturation* saturationEngine;
	G4double timeLimit_;
	G4bool secondPass;
	G4int version_;
	G4double zOff;


};

#endif
