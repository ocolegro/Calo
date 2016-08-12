#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "HGCSSGenParticle.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"

//
SteppingAction::SteppingAction(std::string data) {
	eventAction_ =
			(EventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
	eventAction_->Add(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getStructure());
	saturationEngine = new G4EmSaturation();
	version_ = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getVersion();

	DetectorConstruction*  Detector =
			(DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
	zOff = -0.5 * (Detector->GetCalorSizeZ());
	secondPass = (data == "") ? false : true;


}

//
SteppingAction::~SteppingAction() {
}

//
void SteppingAction::UserSteppingAction(const G4Step* aStep) {

	// get PreStepPoint
	const G4StepPoint *thePreStepPoint = aStep->GetPreStepPoint();
	const G4StepPoint *thePostStepPoint = aStep->GetPostStepPoint();

    G4double bFieldPost[3] = {999,0,0};
    G4double passPostPos[4] = {0,0,0,0};
    G4double passPrePos[4] = {0,0,0,0};
	G4Track* lTrack = aStep->GetTrack();
	G4double kinEng = lTrack->GetKineticEnergy();
	G4int pdgID = lTrack->GetDefinition()->GetPDGEncoding();
	G4int trackID = lTrack->GetTrackID();

	G4VPhysicalVolume* volume = thePreStepPoint->GetPhysicalVolume();

	const G4ThreeVector & postPos = thePostStepPoint->GetPosition();
	const G4ThreeVector & prePos   = thePreStepPoint->GetPosition();



	const G4ThreeVector &preMom = lTrack->GetMomentum() + -1.*aStep->GetDeltaMomentum();
	const G4ThreeVector &postMom = lTrack->GetMomentum();
	/*
    for(int i = 0; i < 3; i++)
    {
    	passPrePos[i] = prePos[i]/cm;
        passPostPos[i] = postPos[i]/cm;
    }

	*/

	//G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(passPrePos, bFieldPre);
	//G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(passPostPos, bFieldPost);



	//G4cout << "The pdgid is" << pdgID << G4endl;
	//G4cout << "The trackid is" << trackID << G4endl;
	/*G4cout << "The volume is " << volume->GetName() << G4endl;
	G4cout << "The preMom is" << preMom[0] << ", " << preMom[1] << ", " << preMom[2] << G4endl;
	G4cout << "The passPrePos is" << passPrePos[0] << ", " << passPrePos[1] << ", " << passPrePos[2] << G4endl;
	G4cout << "The bFieldPre is" << bFieldPre[0] << ", " << bFieldPre[1] << ", " << bFieldPre[2] << G4endl;
*/
   /* double lenUnit = centimeter;

	if(trackID == 1){
		G4cout << "The passPrePos is " << passPrePos[0]/lenUnit << ", " << passPrePos[1]/lenUnit << ", " << passPrePos[2]/lenUnit << G4endl;
		G4cout << "The passPostPos is " << passPostPos[0]/lenUnit << ", " << passPostPos[1]/lenUnit << ", " << passPostPos[2]/lenUnit << G4endl;

		G4cout << "The bFieldPre is " << bFieldPre[0]/gauss  << ", " << bFieldPre[1] /gauss << ", " << bFieldPre[2] /gauss << G4endl;
		G4cout << "The bFieldPost is " << bFieldPost[0] /gauss << ", " << bFieldPost[1] /gauss << ", " << bFieldPost[2]/gauss  << G4endl;

		G4cout << "The preMom is " << preMom[0] << ", " << preMom[1] << ", " << preMom[2] << G4endl;
		G4cout << "The postMom is " << postMom[0] << ", " << postMom[1] << ", " << postMom[2] << G4endl;

	}*/
	for (int i = 0; i <150; i ++)
	{
		for (double j = 0; j < 104; j++){
		double bFieldPre[3] = {999,0,0};
			passPrePos[0] = j/4.*centimeter;
			passPrePos[1] = j/4.*centimeter;
			passPrePos[2] = i*1*centimeter;
			G4cout << "making a field call w/ x,y,z = " << passPrePos[2]/centimeter << G4endl;
			G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(passPrePos, bFieldPre);
		}
	}
}





void SteppingAction::printParticle(G4Track* aTrack)
{
  G4cout << aTrack->GetParticleDefinition()->GetParticleName() << "  "
	<< aTrack->GetDefinition()->GetPDGEncoding() << "  "
	<< aTrack->GetTrackID()<< "  "

	<< aTrack->GetTotalEnergy() << "  "
	<< aTrack->GetKineticEnergy() << "  "
	<< aTrack->GetMomentum().x() << "  "
	<< aTrack->GetMomentum().y() << "  "
	<< aTrack->GetMomentum().z() << "  "
	<< aTrack->GetParticleDefinition()->GetPDGMass() << G4endl;
  return;
}
