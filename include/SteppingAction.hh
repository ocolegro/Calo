#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4EmSaturation.hh"

class EventAction;

class SteppingAction: public G4UserSteppingAction {
public:
	SteppingAction();
	virtual ~SteppingAction();

	void UserSteppingAction(const G4Step*);

	inline bool checkDuplicate(std::vector<double> engVec,double currentEng){
		for (unsigned i = 0; i < engVec.size(); i++){
			if ( (engVec.at(i) - currentEng) < .01) return false;
		}
		return true;
	};
private:
	void printParticle(G4Track* aTrack);
	EventAction *eventAction_;
	//to correct the energy in the scintillator
	G4EmSaturation* saturationEngine;
	G4double zOff;
};

#endif
