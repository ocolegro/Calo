#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "HGCSSGenParticle.hh"

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



	G4double eRawDep = aStep->GetTotalEnergyDeposit()  * MeV;

	stepPDGID = pdgID;
	stepKE = kinEng - aStep->GetDeltaEnergy();

	const G4ThreeVector & postPos = thePostStepPoint->GetpostPos();
	const G4ThreeVector & prePos   = thePreStepPoint->GetpostPos();

	HGCSSGenParticle genPart;
	eventAction_->Detect(eRawDep, volume,lTrack,postPos);
	const G4TrackVector* secondaries= aStep->GetSecondary();

	bool trackEscapes = ((lTrack->GetTrackStatus()!=fAlive
			&& (lTrack->GetKineticEnergy() > 10)
			&& (secondPass || (!eventAction_->firstPass()))
			&& (volume->GetName() == "expHall")
			&& secondaries->size() == 0)
			 );
	if (trackEscapes){

		HGCSSGenParticle escapePart;
		escapePart.vertexKE(lTrack->GetVertexKineticEnergy()); //- aStep->GetDeltaEnergy());
		escapePart.finalKE(lTrack->GetKineticEnergy()); //- aStep->GetDeltaEnergy());

		const G4ThreeVector &p = lTrack->GetVertexMomentumDirection() ;//- aStep->GetDeltaMomentum();
		const G4ThreeVector &pos = lTrack->GetVertexpostPos();
		TVector3 momVec(p[0], p[1], p[2]);
		escapePart.vertexMom(momVec);
		TVector3 posVec(pos[0], pos[1], pos[2] - zOff);
		escapePart.vertexPos(posVec);
		escapePart.pdgid(pdgID);
		escapePart.mass(lTrack->GetDefinition()->GetPDGMass());

		eventAction_->escapevec_.push_back(escapePart);
	}

	if(secondaries->size() > 0){
		G4String theProcessName=secondaries->at(0)->GetCreatorProcess()->GetProcessName();
		if ( (theProcessName == "PhotonInelastic" || theProcessName == "ElectroNuclear" || theProcessName == "PositronNuclear")
				&& (abs(pdgID) == 22 || abs(pdgID) == 11)){

			//sloppy fix for strange geant4 stepping action
			bool checkSecEngs = true;
			bool someHadrons = false;
			for(G4TrackVector::const_iterator i=secondaries->begin(); i!=secondaries->end(); ++i){
				G4Track* iTrack = *i;
				if (iTrack->GetKineticEnergy() < 10) continue;
				if (!checkDuplicate(eventAction_->targetPartEngs,iTrack->GetKineticEnergy()))
					checkSecEngs = false;
				if (abs(iTrack->GetDefinition()->GetPDGEncoding()) != 11 ||
						abs(iTrack->GetDefinition()->GetPDGEncoding()) != 22)
					someHadrons = true;

			}
			if (checkSecEngs && someHadrons){
				eventAction_->hadronicInts = eventAction_->hadronicInts  + 1;
				HGCSSGenParticle targPart;

				for(G4TrackVector::const_iterator i=secondaries->begin(); i!=secondaries->end(); ++i){
					G4Track* iTrack = *i;
					HGCSSGenParticle genPart;
					if (iTrack->GetKineticEnergy() < 10) continue;
						if (targPart.vertexKE() == 0){
							targPart.vertexKE(lTrack->GetKineticEnergy() - aStep->GetDeltaEnergy());


							const G4ThreeVector &p = lTrack->GetMomentum() + -1.*aStep->GetDeltaMomentum();
							const G4ThreeVector &pos = lTrack->GetpostPos();
							if (p.mag() > 0){
								TVector3 momVec(p[0]/p.mag(), p[1]/p.mag(), p[2]/p.mag());
								targPart.vertexMom(momVec);
							}
							else{
								G4cout << "This vertex momentum was broken" << G4endl;
							}
							TVector3 posVec(pos[0], pos[1], pos[2] - zOff);
							targPart.vertexPos(posVec);
							targPart.pdgid(pdgID);
							targPart.layer(-eventAction_->hadronicInts);
							eventAction_->incvec_.push_back(targPart);

						}
						eventAction_->targetPartEngs.push_back(iTrack->GetKineticEnergy());

								genPart.vertexKE(iTrack->GetKineticEnergy());
								const G4ThreeVector &p = iTrack->GetMomentumDirection();
								const G4ThreeVector &pos = iTrack->GetpostPos();
								TVector3 momVec(p[0], p[1], p[2]);
								genPart.vertexMom(momVec);
								TVector3 posVec(pos[0], pos[1], pos[2] - zOff);
								genPart.vertexPos(posVec);
								genPart.mass(iTrack->GetDefinition()->GetPDGMass());
								genPart.pdgid(iTrack->GetDefinition()->GetPDGEncoding());
								genPart.layer(eventAction_->hadronicInts);
								eventAction_->hadvec_.push_back(genPart);


						//}
				}

            //Get Magnetic Field Vectors at pre and post postPoss
                for(int i = 0; i < 3; i++)
            {
                passPostPos[i] = postPos[i];
                passPrePos[i] = prePos[i];
            }
            genPart.setBField(bFieldPost[0],bFieldPost[1], bFieldPost[2]);

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
