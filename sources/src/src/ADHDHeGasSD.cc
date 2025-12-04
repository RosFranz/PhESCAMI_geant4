// $Id: ADHDHeGasSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDHeGasSD.cc
/// \brief Implementation of the ADHDHeGasSD class

#include "ADHDHeGasSD.hh"
#include "ADHDHeGasHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "MyDet.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


ADHDHeGasSD::ADHDHeGasSD(G4String name): G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
  G4String HCname = "HeGasColl";
  collectionName.insert(HCname);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


ADHDHeGasSD::~ADHDHeGasSD(){}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void ADHDHeGasSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new ADHDHeGasHitsCollection(SensitiveDetectorName,collectionName[0]);

  if (fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

  hce->AddHitsCollection(fHCID,fHitsCollection);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool ADHDHeGasSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  MyDet *MD = MyDet::gethead();

  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  G4StepPoint* preStepPoint = step->GetPreStepPoint();

  G4Track* trk = (G4Track*) step->GetTrack();
  G4TrackStatus trkstat = trk->GetTrackStatus();

  G4double charge = postStepPoint->GetCharge();
  G4double mass = postStepPoint->GetMass();

  if (trkstat==fStopButAlive) {
    
    // ---- CAPTURE OF A NEGATIVE PARTICLE WITH MASS > MYMASS (900 MeV/c^2) ---
    if (charge < 0 && mass > MD->MyMass){
      
      MD->MyCaptureT.push_back(postStepPoint->GetGlobalTime());

      MD->MyCaptureX.push_back((postStepPoint->GetPosition())[0]);
      MD->MyCaptureY.push_back((postStepPoint->GetPosition())[1]);
      MD->MyCaptureZ.push_back((postStepPoint->GetPosition())[2]);

      MD->MyCaptureE.push_back(preStepPoint->GetKineticEnergy());

      MD->MyCaptureF.push_back(mass);
      MD->MyCaptureM.push_back(preStepPoint->GetMass());

      MD->MyCaptureCopyNo.push_back(postStepPoint->GetTouchableHandle()->GetCopyNumber());

      // DELAY add in capture
      trk->SetGlobalTime(trk->GetGlobalTime()+MD->MyDelay);
    }

  }
  G4double edep = step->GetTotalEnergyDeposit();

  G4double hitTime = 0.5*(preStepPoint->GetGlobalTime() + postStepPoint->GetGlobalTime());
  G4ThreeVector worldPos = 0.5*(preStepPoint->GetPosition()+postStepPoint->GetPosition());


  if(edep>0){

    G4TouchableHistory* touchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4int copyNo = touchable->GetVolume()->GetCopyNo();

    // Check on the kinetic energy before and after one Helium tank
    //if(trk->GetParentID()==0&&copyNo==12 && (step->IsFirstStepInVolume()||step->IsLastStepInVolume())) G4cout << "Primary Kinetic Energy: " << postStepPoint->GetKineticEnergy() << G4endl;
    
    // --- Checking for possible hit merging --- 

    // check if a spacetime-nearby hit is existing
    G4int imerge = 0;

    for (G4int i=0; i < int(fHitsCollection->entries()) ; i++){

      imerge = 1; // default is merging

      if( fabs( hitTime - ((*fHitsCollection)[i]->GetTime()) ) > MD->MyHeTres ) imerge=0; //not same hit 
      if( imerge == 0 ) continue; // try next hit

      double dist2 = worldPos.diff2( (*fHitsCollection)[i]->GetPos() );

      if( dist2 > ( (MD->MyHeDist)*(MD->MyHeDist) ) ) imerge=0; //not same hit
      if(imerge == 0) continue; // try next hit

      // here the Hit is narrow in spacetime -> HIT MERGING
      double toten = edep + ((*fHitsCollection)[i]->GetEne());

      ((*fHitsCollection)[i]) -> SetPos( ( edep * worldPos + ((*fHitsCollection)[i]->GetEne()) * ((*fHitsCollection)[i]->GetPos()) ) / toten );
      ((*fHitsCollection)[i]) -> SetTime( ( edep * hitTime + ((*fHitsCollection)[i]->GetEne()) * ((*fHitsCollection)[i]->GetTime()) ) / toten );
      ((*fHitsCollection)[i]) -> SetEne( toten ); 
      ((*fHitsCollection)[i]) -> SetCopyNo( copyNo );
      break; // hit merged

    }
    
    // create a new hit and set it to the collection
    if (imerge==0){

      ADHDHeGasHit* hit = new ADHDHeGasHit(copyNo);

      G4VPhysicalVolume* physical = touchable->GetVolume();
      hit->SetLogV(physical->GetLogicalVolume());

      hit->SetEne(edep);
      hit->SetPos(worldPos);
      hit->SetTime(hitTime);
      hit->SetCopyNo(copyNo);

      // if(copyNo==12) hit->Print(); //event 16
      // if(copyNo==37) hit->Print(); // generic event

      fHitsCollection->insert(hit);

    }

  }

  return true;
}