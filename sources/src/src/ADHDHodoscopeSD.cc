// $Id: ADHDHodoscopeSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDHodoscopeSD.cc
/// \brief Implementation of the ADHDHodoscopeSD class

#include "ADHDHodoscopeSD.hh"
#include "ADHDHodoscopeHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "MyDet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDHodoscopeSD::ADHDHodoscopeSD(G4String name): G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
  G4String HCname = "hodoscopeColl";
  collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDHodoscopeSD::~ADHDHodoscopeSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDHodoscopeSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new ADHDHodoscopeHitsCollection (SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ADHDHodoscopeSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  MyDet *MD = MyDet::gethead();
  G4double edep = step->GetTotalEnergyDeposit();

  if (edep==0.) return true;
  
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  G4Track* trk = (G4Track*) step->GetTrack();
  G4TrackStatus trkstat = trk->GetTrackStatus();

  G4double charge = postStepPoint->GetCharge();
  G4double mass = postStepPoint->GetMass();

  if (trkstat==fStopButAlive) {

    if (charge < 0 && mass > MD->MyMass){

      MD->HodCaptureT.push_back(postStepPoint->GetGlobalTime());

      MD->HodCaptureX.push_back((postStepPoint->GetPosition())[0]);
      MD->HodCaptureY.push_back((postStepPoint->GetPosition())[1]);
      MD->HodCaptureZ.push_back((postStepPoint->GetPosition())[2]);

      MD->HodCaptureE.push_back(preStepPoint->GetKineticEnergy());
      MD->HodCaptureF.push_back(mass);
      MD->HodCaptureM.push_back(preStepPoint->GetMass());

      G4LogicalVolume *lv = postStepPoint->GetPhysicalVolume()->GetLogicalVolume();

      if( lv->GetName() == "hodoscope2Logical") MD->HodCaptureVolume.push_back(1); //inner scint == 1
      if( lv->GetName() == "hodoscope3Logical") MD->HodCaptureVolume.push_back(2); //outer scint == 2

    }

  }


  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4int copyNo = touchable->GetVolume()->GetCopyNo();

  G4double hitTime = 0.5*(preStepPoint->GetGlobalTime() + postStepPoint->GetGlobalTime());
  G4ThreeVector worldPos = 0.5*(preStepPoint->GetPosition()+postStepPoint->GetPosition());

  // --- ATTEMPT TO MERGE THE HITS ---
  // check if a spacetime-nearby hit is existing
  G4int iadd = 0;

  for (G4int i=0 ; i < int(fHitsCollection->entries()) ; i++){

    iadd = 1; // default is adding

    if(fabs(hitTime-((*fHitsCollection)[i]->GetTime())) > MD->MySCTres) iadd=0; //not same hit
    if(iadd == 0) continue; // try next hit
    
    double dist2 = worldPos.diff2((*fHitsCollection)[i]->GetPos());

    if(dist2>((MD->MySCDist)*(MD->MySCDist))) iadd=0; //not same hit
    if(iadd == 0) continue; // try next hit

    // here the Hit is narrow in spacetime -> HIT MERGING
    double toten = edep+(*fHitsCollection)[i]->GetEne();

    ((*fHitsCollection)[i])->SetPos( ( edep * worldPos + ((*fHitsCollection)[i]->GetEne()) * ((*fHitsCollection)[i]->GetPos()) ) /toten);
    ((*fHitsCollection)[i])->SetTime( ( edep * hitTime + ((*fHitsCollection)[i]->GetEne()) * ((*fHitsCollection)[i]->GetTime()) ) /toten);
    ((*fHitsCollection)[i])->SetEne(toten);
    ((*fHitsCollection)[i])->SetCopyNo(copyNo);

    break; // hit merged

  }
  
  // create a new hit and set it to the collection
  if (iadd==0){

    ADHDHodoscopeHit* hit = new ADHDHodoscopeHit(copyNo);
    G4VPhysicalVolume* physical = touchable->GetVolume();
    hit->SetLogV(physical->GetLogicalVolume());

    hit->SetEne(edep);
    hit->SetPos(worldPos);
    hit->SetTime(hitTime);
    hit->SetCopyNo(copyNo);
    
    fHitsCollection->insert(hit);

  }  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
