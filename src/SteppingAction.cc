#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "HGCSSGenParticle.hh"
#include "G4TransportationManager.hh"
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

	double bFieldPre[3] = {0,0,0};
    double bFieldPost[3] = {0,0,0};
    double passPostPos[4] = {0,0,0,0};
    double passPrePos[4] = {0,0,0,0};
	G4Track* lTrack = aStep->GetTrack();
	G4double kinEng = lTrack->GetKineticEnergy();
	G4int pdgID = lTrack->GetDefinition()->GetPDGEncoding();
	G4VPhysicalVolume* volume = thePreStepPoint->GetPhysicalVolume();

	const G4ThreeVector & postPos = thePostStepPoint->GetPosition();
	const G4ThreeVector & prePos   = thePreStepPoint->GetPosition();



	G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(passPrePos, bFieldPre);
	G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(passPostPos, bFieldPost);


	const G4ThreeVector &preMom = lTrack->GetMomentum() + -1.*aStep->GetDeltaMomentum();
	const G4ThreeVector &postMom = lTrack->GetMomentum();

	G4cout << "The pdgid is" << pdgID << G4endl;
	G4cout << "The volume is " << volume->GetName() << G4endl;
	G4cout << "The preMom is" << preMom[0] << ", " << preMom[1] << ", " << preMom[2] << G4endl;
	G4cout << "The passPrePos is" << passPrePos[0] << ", " << passPrePos[1] << ", " << passPrePos[2] << G4endl;
	G4cout << "The bFieldPre is" << bFieldPre[0] << ", " << bFieldPre[1] << ", " << bFieldPre[2] << G4endl;
	G4cout << "The postMom is" << postMom[0] << ", " << postMom[1] << ", " << postMom[2] << G4endl;
	G4cout << "The passPostPos is" << passPostPos[0] << ", " << passPostPos[1] << ", " << passPostPos[2] << G4endl;
	G4cout << "The bFieldPost is" << bFieldPost[0] << ", " << bFieldPost[1] << ", " << bFieldPost[2] << G4endl;

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
