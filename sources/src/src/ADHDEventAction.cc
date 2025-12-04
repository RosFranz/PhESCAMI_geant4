// $Id: ADHDEventAction.cc 94486 2015-11-19 08:33:37Z gcosmo $
//
/// \file ADHDEventAction.cc
/// \brief Implementation of the ADHDEventAction class

#include "ADHDEventAction.hh"
#include "ADHDHodoscopeHit.hh"
#include "ADHDHeGasHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "MyDet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDEventAction::ADHDEventAction(): G4UserEventAction(),
  fHHC2ID(-1), fHHC3ID(-1), fHGHCID(-1),
  fHod2T(), fHod2X(), fHod2Y(), fHod2Z(), fHod2E(), fHod2CopyNo(),
  fHod3T(), fHod3X(), fHod3Y(), fHod3Z(), fHod3E(), fHod3CopyNo(),
  fHGasT(), fHGasX(), fHGasY(), fHGasZ(), fHGasE(), fHGasCopyNo(),
  fHod3_E_slabPrompt(), fHod3_E_slabAcap(), fHod3_E_slabDelayed(),
  fHod3_E_slabMyPrompt(), fHod3_E_slabLate(),
  fHod2_E_slabPrompt(), fHod2_E_slabAcap(), fHod2_E_slabDelayed(),
  fHod2_E_slabMyPrompt(), fHod2_E_slabLate(),
  fHeCal_E_tankPrompt(),  fHeCal_E_tankAcap(), fHeCal_E_tankDelayed(),
  fHeCal_E_tankMyPrompt(),  fHeCal_E_tankLate()
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(10000);

  MyDet *MD = MyDet::gethead();
  
  // initialize the vectors
  MD->MyCaptureT.resize(1, 0.);
  MD->MyCaptureX.resize(1, 0.);
  MD->MyCaptureY.resize(1, 0.);
  MD->MyCaptureZ.resize(1, 0.);
  MD->MyCaptureE.resize(1, 0.);
  MD->MyCaptureF.resize(1, 0.);
  MD->MyCaptureM.resize(1, 0.);
  MD->MyCaptureCopyNo.resize(1, 0.);

  MD->MyCaptureDaughterM.resize(1, 0.);
  MD->MyCaptureDaughterX.resize(1, 0.);
  MD->MyCaptureDaughterY.resize(1, 0.);
  MD->MyCaptureDaughterZ.resize(1, 0.);
  MD->MyCaptureDaughterPDGcode.resize(1, 0.);

  MD->HodCaptureT.resize(1, 0.);
  MD->HodCaptureX.resize(1, 0.);
  MD->HodCaptureY.resize(1, 0.);
  MD->HodCaptureZ.resize(1, 0.);
  MD->HodCaptureE.resize(1, 0.);
  MD->HodCaptureF.resize(1, 0.);
  MD->HodCaptureM.resize(1, 0.);
  MD->HodCaptureVolume.resize(1, 0.);



  fHod2T.resize(1, 0.);
  fHod2X.resize(1, 0.);
  fHod2Y.resize(1, 0.);
  fHod2Z.resize(1, 0.);
  fHod2E.resize(1, 0.);
  fHod2CopyNo.resize(1, 0.);

  fHod3T.resize(1, 0.);
  fHod3X.resize(1, 0.);
  fHod3Y.resize(1, 0.);
  fHod3Z.resize(1, 0.);
  fHod3E.resize(1, 0.);
  fHod3CopyNo.resize(1, 0.);

  fHGasT.resize(1, 0.);
  fHGasX.resize(1, 0.);
  fHGasY.resize(1, 0.);
  fHGasZ.resize(1, 0.);
  fHGasE.resize(1, 0.);
  fHGasCopyNo.resize(1, 0.);

  //number of slabs = 64*6
  fHod3_E_slabPrompt.resize(384, 0.);
  fHod3_E_slabAcap.resize(384, 0.);
  fHod3_E_slabDelayed.resize(384, 0.);
  fHod3_E_slabMyPrompt.resize(384, 0.);
  fHod3_E_slabLate.resize(384, 0.);

  fHod2_E_slabPrompt.resize(384, 0.);
  fHod2_E_slabAcap.resize(384, 0.);
  fHod2_E_slabDelayed.resize(384, 0.);
  fHod2_E_slabMyPrompt.resize(384, 0.);
  fHod2_E_slabLate.resize(384, 0.);

  //number of tanks = 75
  fHeCal_E_tankPrompt.resize(75, 0.);
  fHeCal_E_tankAcap.resize(75, 0.);
  fHeCal_E_tankDelayed.resize(75, 0.);
  fHeCal_E_tankMyPrompt.resize(75, 0.);
  fHeCal_E_tankLate.resize(75, 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDEventAction::~ADHDEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDEventAction::BeginOfEventAction(const G4Event*)
{
  if (fHHC2ID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fHHC2ID = sdManager->GetCollectionID("hodoscope2/hodoscopeColl");
    fHHC3ID = sdManager->GetCollectionID("hodoscope3/hodoscopeColl");
    fHGHCID = sdManager->GetCollectionID("HeGas/HeGasColl");
  }

  MyDet *MD = MyDet::gethead();

  MD->MyCaptureT.resize(1, 0.);
  MD->MyCaptureX.resize(1, 0.);
  MD->MyCaptureY.resize(1, 0.);
  MD->MyCaptureZ.resize(1, 0.);
  MD->MyCaptureE.resize(1, 0.);
  MD->MyCaptureF.resize(1, 0.);
  MD->MyCaptureM.resize(1, 0.);
  MD->MyCaptureCopyNo.resize(1, 0.);

  MD->MyCaptureDaughterM.resize(1, 0.);
  MD->MyCaptureDaughterX.resize(1, 0.);
  MD->MyCaptureDaughterY.resize(1, 0.);
  MD->MyCaptureDaughterZ.resize(1, 0.);
  MD->MyCaptureDaughterPDGcode.resize(1, 0.);

  MD->HodCaptureT.resize(1, 0.);
  MD->HodCaptureX.resize(1, 0.);
  MD->HodCaptureY.resize(1, 0.);
  MD->HodCaptureZ.resize(1, 0.);
  MD->HodCaptureE.resize(1, 0.);
  MD->HodCaptureF.resize(1, 0.);
  MD->HodCaptureM.resize(1, 0.);
  MD->HodCaptureVolume.resize(1, 0.);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ADHDEventAction::HeGasHitCluster(const G4Event* event, int fHHC)
{
  //hit collection of the event
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  ADHDHeGasHitsCollection* hHC = static_cast<ADHDHeGasHitsCollection*>(hce->GetHC(fHHC));

  //resizing the vectors on the hit collection entries
  fHGasT.resize(hHC->entries(), 0.);
  fHGasX.resize(hHC->entries(), 0.);
  fHGasY.resize(hHC->entries(), 0.);
  fHGasZ.resize(hHC->entries(), 0.);
  fHGasE.resize(hHC->entries(), 0.);
  fHGasCopyNo.resize(hHC->entries(), 0.);

  //intializing the vectors
  for(int  i=0 ; i < int(hHC->entries()) ; i++){
    fHGasT.at(i) =0;
    fHGasX.at(i) =0;
    fHGasY.at(i) =0;
    fHGasZ.at(i) =0;
    fHGasE.at(i) =0;
    fHGasCopyNo.at(i) =0;
  }


  //filling and temporally sorting the vectors
  for (int i = 0 ; i < int(hHC->entries()) ; i++){

    double ee = (*hHC)[i]->GetEne();
    double tt = (*hHC)[i]->GetTime();

    if (tt<=0. || ee<=0.) continue;
    
    double xx = (*hHC)[i]->GetPos()[0];
    double yy = (*hHC)[i]->GetPos()[1];
    double zz = (*hHC)[i]->GetPos()[2];

    int copyno = (*hHC)[i]->GetCopyNo();

    for (int ist = int(hHC->entries()) -1 ; int(ist)>=0 ; ist--){
      
      int iplace=0;
      
      if(fHGasT.at(ist)==0 && ist>0 && fHGasT[ist-1]==0) continue; // c'e' ancora spazio vai avanti
      if(fHGasT.at(ist)>tt && ist < int(hHC->entries())-1 && int(ist)>0 && fHGasT[ist-1]<=tt) iplace = 2; // fatti spazio nel mezzo
      if(fHGasT.at(ist)>tt && ist==0) iplace = 2; // fatti spazio in cima
      if(fHGasT.at(ist)>tt && ist == int(hHC->entries())-1 && fHGasT[ist-1]<=tt) iplace = 1; // spazio finito devi sostituire l'ultimo
      if(fHGasT.at(ist)==0 && ist==0) iplace = 1; // arrivato all'inizio vettore tutto vuoto piazzalo li 
      if(fHGasT.at(ist)==0 && ist>0 && fHGasT[ist-1]<=tt) iplace = 1; // PIAZZA

      if(iplace==2){

        // fai spazio e piazza in mezzo
        for (int ist2 = int(hHC->entries())-2;ist2>=ist;ist2--){

          fHGasT.at(ist2+1)=fHGasT.at(ist2);
          fHGasX.at(ist2+1)=fHGasX.at(ist2);
          fHGasY.at(ist2+1)=fHGasY.at(ist2);
          fHGasZ.at(ist2+1)=fHGasZ.at(ist2);
          fHGasE.at(ist2+1)=fHGasE.at(ist2);
          fHGasCopyNo.at(ist2+1)=fHGasCopyNo.at(ist2);

          }

        iplace=1;

      }

      if(iplace==1){
        fHGasT.at(ist)=tt;
        fHGasX.at(ist)=xx;
        fHGasY.at(ist)=yy;
        fHGasZ.at(ist)=zz;
        fHGasE.at(ist)=ee;
        fHGasCopyNo.at(ist)=copyno;
        break; // piazzato esci dal ciclo
      }

    }

  }  
  G4int NHits=0; 
  NHits=fHGasT.size();
  return NHits;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ADHDEventAction::HodHitCluster(const G4Event* event, int fHHC, int layer){
  //Hits collection of this event
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  ADHDHodoscopeHitsCollection* hHC = static_cast<ADHDHodoscopeHitsCollection*>(hce->GetHC(fHHC));
  

  std::vector<G4double>*fHodE;
  std::vector<G4double>*fHodT;
  std::vector<G4double>*fHodX;
  std::vector<G4double>*fHodY;
  std::vector<G4double>*fHodZ;
  std::vector<G4int>*fHodCopyNo;

  //INNER LAYER "Hodoscope 2" 
  if(layer==2) {
    
    //resizing the vectors on the hit collection entries
    fHod2T.resize(hHC->entries(), 0.);
    fHod2X.resize(hHC->entries(), 0.);
    fHod2Y.resize(hHC->entries(), 0.);
    fHod2Z.resize(hHC->entries(), 0.);
    fHod2E.resize(hHC->entries(), 0.);
    fHod2CopyNo.resize(hHC->entries(), 0.);

    //intializing the vectors
    for(int i=0 ; i < int(hHC->entries()) ; i++){
      fHod2T.at(i) =0;
      fHod2X.at(i) =0;
      fHod2Y.at(i) =0;
      fHod2Z.at(i) =0;
      fHod2E.at(i) =0;
      fHod2CopyNo.at(i) =0;
    }

    fHodE = &GetHod2E();
    fHodT = &GetHod2T();
    fHodX = &GetHod2X();
    fHodY = &GetHod2Y();
    fHodZ = &GetHod2Z();
    fHodCopyNo = &GetHod2CopyNo();
  }

  //OUTER LAYER "Hodoscope 3" 
  else{

    //reisizing the vectors on the hit collection entries
    fHod3T.resize(hHC->entries(), 0.);
    fHod3X.resize(hHC->entries(), 0.);
    fHod3Y.resize(hHC->entries(), 0.);
    fHod3Z.resize(hHC->entries(), 0.);
    fHod3E.resize(hHC->entries(), 0.);
    fHod3CopyNo.resize(hHC->entries(), 0.);

    //intializing the vectors
    for(int i = 0; i < int(hHC->entries()); i++){
      fHod3E.at(i) = 0.;
      fHod3T.at(i) = 0.;
      fHod3X.at(i) = 0.;
      fHod3Y.at(i) = 0.;
      fHod3Z.at(i) = 0.;
      fHod3CopyNo.at(i) = 0.;
    }

    fHodE = &GetHod3E();
    fHodT = &GetHod3T();
    fHodX = &GetHod3X();
    fHodY = &GetHod3Y();
    fHodZ = &GetHod3Z();
    fHodCopyNo = &GetHod3CopyNo();
  }


  //filling and temporally sorting the vectors
  for (int i=0;i<int(hHC->entries());i++){

    //getting the info from the hits collection
    double ee = (*hHC)[i]->GetEne();
    double tt = (*hHC)[i]->GetTime();

    if (tt<=0. || ee<=0.) continue;

    double xx = ((*hHC)[i]->GetPos())[0];
    double yy = ((*hHC)[i]->GetPos())[1];
    double zz = ((*hHC)[i]->GetPos())[2];

    G4int copyno = (*hHC)[i]->GetCopyNo();


    //second cycle for sorting temporally the vectors
    for (int ist=int(hHC->entries())-1;int(ist)>=0;ist--){

      int iplace=0;

      if((*fHodT).at(ist)==0 && ist>0 && (*fHodT).at(ist-1)==0) continue;                            // c'e' ancora spazio vai avanti
      if((*fHodT).at(ist)>tt && ist<int(hHC->entries())-1 && ist>0 && (*fHodT).at(ist-1)<=tt) iplace = 2; // fatti spazio nel mezzo
      if((*fHodT).at(ist)>tt && ist==0) iplace = 2;                                                  // fatti spazio in cima
      if((*fHodT).at(ist)>tt && ist==int(hHC->entries())-1 && (*fHodT).at(ist-1)<=tt) iplace = 1;         // spazio finito devi sostituire l'ultimo
      if((*fHodT).at(ist)==0 && ist==0) iplace = 1;                                                  // arrivato all'inizio vettore tutto vuoto piazzalo li 
      if((*fHodT).at(ist)==0 && ist>0 && (*fHodT).at(ist-1)<=tt) iplace = 1;                         // PIAZZA

      if(iplace==2){                                                                                 // fai spazio e piazza in mezzo

        for (int ist2=int(hHC->entries())-2;ist2>=ist;ist2--){

          (*fHodT).at(ist2+1)=(*fHodT).at(ist2);
          (*fHodE).at(ist2+1)=(*fHodE).at(ist2);
          (*fHodX).at(ist2+1)=(*fHodX).at(ist2);
          (*fHodY).at(ist2+1)=(*fHodY).at(ist2);
          (*fHodZ).at(ist2+1)=(*fHodZ).at(ist2);
          (*fHodCopyNo).at(ist2+1)=(*fHodCopyNo).at(ist2);
          iplace=1;

        }

      }

      if(iplace==1){
        
        (*fHodT).at(ist)=tt;
        (*fHodE).at(ist)=ee;
        (*fHodX).at(ist)=xx;
        (*fHodY).at(ist)=yy;
        (*fHodZ).at(ist)=zz;
        (*fHodCopyNo).at(ist)=copyno;
        break;                                                                                        // piazzato esci dal ciclo

      }

    }
    
  }

  G4int NHod=0;
  NHod=fHodT->size();
  return NHod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDEventAction::EndOfEventAction(const G4Event* event)
{
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl; 
    G4Exception("ADHDEventAction::EndOfEventAction()", "ADHDCode001", JustWarning, msg);
    return;
  }

  MyDet *MD = MyDet::gethead();
    
  // information on primary particle
  G4PrimaryVertex* priv = event->GetPrimaryVertex();
  G4PrimaryParticle* primary = priv->GetPrimary();
  G4double MCEnergy = primary->GetKineticEnergy();
  G4double MCMomentum = primary->GetTotalMomentum();
  G4double MCMass = primary->GetMass();
  G4double MCCharge = primary->GetCharge();

  //Check which particle is the primary
  if(event->GetEventID() == 0){
    NofNucleons = 1.;
    //check for antiparticles
    if(primary->GetPDGcode()<0.) MD->antiparticle = 1;
    // for nueclei (deutons) get the number of nucloens
    if(fabs(primary->GetPDGcode()) == 1000010020) NofNucleons = 2.;
    //He4
    if(fabs(primary->GetPDGcode()) == 1000020040) NofNucleons = 4.;
    //C12
    if(fabs(primary->GetPDGcode()) == 1000060120) NofNucleons = 12.;
  }
  
  


  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 


  //sorting temporarly HeCal vectors
  int NHG = HeGasHitCluster(event, fHGHCID);
    
  //trigger on gas for granular tree "Hits"
  if(NHG>1) { 

    //temporarly sorting the hodoscopes vectors
    int NHod3 = HodHitCluster(event, fHHC3ID, 3);
    int NHod2 = HodHitCluster(event, fHHC2ID, 2);

    // counting the triggers
    MD->TRGEvnum=MD->TRGEvnum+1;

    //filling ntuple's columns
    // vector columns are automatically filled
    analysisManager->FillNtupleDColumn(0, 0, MCEnergy/NofNucleons);
    analysisManager->FillNtupleDColumn(0, 1, MCMomentum/NofNucleons);
    analysisManager->FillNtupleDColumn(0, 9, float(MD->MCEvnum)/MD->GenSurface);
    analysisManager->FillNtupleIColumn(0, 10, NHod2);
    analysisManager->FillNtupleIColumn(0, 11, NHod3);
    analysisManager->FillNtupleIColumn(0, 12, NHG);
    analysisManager->FillNtupleDColumn(0, 13, (primary->GetMomentum())[0]/NofNucleons);
    analysisManager->FillNtupleDColumn(0, 14, (primary->GetMomentum())[1]/NofNucleons);
    analysisManager->FillNtupleDColumn(0, 15, (primary->GetMomentum())[2]/NofNucleons);
    analysisManager->FillNtupleDColumn(0, 16, (priv->GetPosition())[0]);
    analysisManager->FillNtupleDColumn(0, 17, (priv->GetPosition())[1]);
    analysisManager->FillNtupleDColumn(0, 18, (priv->GetPosition())[2]);
    analysisManager->FillNtupleDColumn(0, 35, MCMass);
    analysisManager->FillNtupleDColumn(0, 36, MCCharge);
    analysisManager->FillNtupleIColumn(0, 45, event->GetEventID());

    analysisManager->AddNtupleRow(0); 

    //number of slabs = 64*6
    G4int NScint = 64*6;
    fHod3_E_slabPrompt.resize(NScint, 0.);
    fHod3_E_slabAcap.resize(NScint, 0.);
    fHod3_E_slabDelayed.resize(NScint, 0.);

    fHod3_E_slabMyPrompt.resize(NScint, 0.);
    fHod3_E_slabLate.resize(NScint, 0.);

    fHod2_E_slabPrompt.resize(NScint, 0.);
    fHod2_E_slabAcap.resize(NScint, 0.);
    fHod2_E_slabDelayed.resize(NScint, 0.);

    fHod2_E_slabMyPrompt.resize(NScint, 0.);
    fHod2_E_slabLate.resize(NScint, 0.);

    for(int i=0; i<NScint; i++){
      fHod3_E_slabPrompt.at(i)=0;
      fHod3_E_slabAcap.at(i)=0;
      fHod3_E_slabDelayed.at(i)=0;

      fHod3_E_slabMyPrompt.at(i)=0;
      fHod3_E_slabLate.at(i)=0;

      fHod2_E_slabPrompt.at(i)=0;
      fHod2_E_slabAcap.at(i)=0;
      fHod2_E_slabDelayed.at(i)=0;

      fHod2_E_slabMyPrompt.at(i)=0;
      fHod2_E_slabLate.at(i)=0;
    }

    //number of tanks = 75
    G4int NTanks = 75;
    fHeCal_E_tankPrompt.resize(NTanks, 0.);
    fHeCal_E_tankAcap.resize(NTanks, 0.);
    fHeCal_E_tankDelayed.resize(NTanks, 0.);

    fHeCal_E_tankMyPrompt.resize(NTanks, 0.);
    fHeCal_E_tankLate.resize(NTanks, 0.);

    for(int i=0; i<NTanks; i++){
      fHeCal_E_tankPrompt.at(i)=0;
      fHeCal_E_tankAcap.at(i)=0;
      fHeCal_E_tankDelayed.at(i)=0;

      fHeCal_E_tankMyPrompt.at(i)=0;
      fHeCal_E_tankLate.at(i)=0;
    }


    //---- REONSTRUCTED TREE ----//

    // >>> ANTIPARTICLES <<< //

    // at least one capture in HeCal
    // if( MD->MyCaptureT.size() > 0 ) { // avoid any condition in order to determine the trigger logic
    
    if( MD->antiparticle == 1 ) { 

      // initialization
      G4double HeCal_prompt_E=0, HeCal_Acap_E=0, HeCal_Delayed_E=0;
      G4double HeCal_Myprompt_E=0, HeCal_Late_E=0;
      G4double Hod3_prompt_E=0, Hod3_Acap_E=0, Hod3_Delayed_E=0;
      G4double Hod3_Myprompt_E=0, Hod3_Late_E=0;
      G4double Hod2_prompt_E=0, Hod2_Acap_E=0, Hod2_Delayed_E=0;
      G4double Hod2_Myprompt_E=0, Hod2_Late_E=0;
      G4double beta=0;
      G4int Hod2Hit_prompt=0, Hod2Hit_Acap=0, Hod2Hit_Delayed=0;
      G4int Hod2Hit_Myprompt=0, Hod2Hit_Late=0;
      G4int Hod3Hit_prompt=0, Hod3Hit_Acap=0, Hod3Hit_Delayed=0;
      G4int Hod3Hit_Myprompt=0, Hod3Hit_Late=0;
      int CopyNo = 999;
      int MyCopyNo = 0;

      // Incidence on ToF shell (F.Nozzoli)
      // G4float Cosi[3];



      // ---- FILLING NTUPLE ---- //
      //                       (ntuple ID, column ID, value)  
      analysisManager->FillNtupleDColumn(1, 0, MCEnergy/NofNucleons);
      analysisManager->FillNtupleDColumn(1, 1, MCMomentum/NofNucleons);
      analysisManager->FillNtupleDColumn(1, 9, float(MD->MCEvnum)/MD->GenSurface);
      analysisManager->FillNtupleIColumn(1, 10, NHod2);
      analysisManager->FillNtupleIColumn(1, 11, NHod3);
      analysisManager->FillNtupleIColumn(1, 12, NHG);
      analysisManager->FillNtupleDColumn(1, 13, (primary->GetMomentum())[0]/NofNucleons);
      analysisManager->FillNtupleDColumn(1, 14, (primary->GetMomentum())[1]/NofNucleons);
      analysisManager->FillNtupleDColumn(1, 15, (primary->GetMomentum())[2]/NofNucleons);
      analysisManager->FillNtupleDColumn(1, 16, (priv->GetPosition())[0]);
      analysisManager->FillNtupleDColumn(1, 17, (priv->GetPosition())[1]);
      analysisManager->FillNtupleDColumn(1, 18, (priv->GetPosition())[2]);
      analysisManager->FillNtupleDColumn(1, 35, MCMass);
      analysisManager->FillNtupleDColumn(1, 36, MCCharge);
      analysisManager->FillNtupleIColumn(1, 45, event->GetEventID());
      

      // --- CYCLE ON THE CAPTURES --- //
      // entrance in the for cycle is granted by the CaptureT size that is always at least 1
      // In the absence of a capture the value inside CaptureT is 0
      // The reconstruction for anti particles is still made here, even if no capture happens
      
      for(int k=0; k < int(MD->MyCaptureT.size()) ; k++){
        
        double time = MD->MyCaptureT.at(k);  // = 0 if no capture
        if(time==0 && MD->MyCaptureT.size()>1) continue; // if one capture happens, skip the first entry in the vector
        double NextCaptureTime = 0;
        if( k+1 < int(MD->MyCaptureT.size()) && time!=0 ) NextCaptureTime = MD->MyCaptureT.at(k+1);


        // if no capture, prompt is < 50 [ns] 
        // printf("Check on the size of Hod3T %ld value of NHG: %d value of NHod3: %d\n", fHod3T.size(), NHG, NHod3);
        // G4RunManager::GetRunManager()->rndmSaveThisEvent();


        double time_th = 50 + fHod3T.at(0); // 50 [ns] after the first hit in the hodoscope 3


        // ---- HECAL ---- //

          //prompt energy
        for(int j=0; j<NHG; j++){
            int copy = fHGasCopyNo.at(j);
            if(fHGasE.at(j)>0. && fHGasT.at(j)<time) fHeCal_E_tankPrompt.at(copy) += fHGasE.at(j);
            if(fHGasE.at(j)>0. && fHGasT.at(j)<time_th) fHeCal_E_tankMyPrompt.at(copy) += fHGasE.at(j);
        }
        for(int tank=0; tank<NTanks; tank++) {
          if(fHeCal_E_tankPrompt.at(tank)>HeCal_prompt_E) {
            HeCal_prompt_E = fHeCal_E_tankPrompt.at(tank);
            CopyNo = tank;
          }
          if(fHeCal_E_tankMyPrompt.at(tank)>HeCal_Myprompt_E) {
            HeCal_Myprompt_E = fHeCal_E_tankMyPrompt.at(tank);
            MyCopyNo = tank;
          }
        }

      
        analysisManager->FillNtupleIColumn(1, 48, CopyNo); // HeCal CopyNo for the MAX prompt energy release 
        analysisManager->FillNtupleDColumn(1, 49, HeCal_prompt_E);
        analysisManager->FillNtupleIColumn(1, 100, MyCopyNo); // HeCal CopyNo for the MAX prompt energy release 
        analysisManager->FillNtupleDColumn(1, 101, HeCal_Myprompt_E);


          //after capture and before delayed events
        for(int j=0; j<NHG; j++){
          int copy = fHGasCopyNo.at(j);
          if(fHGasE.at(j)>0. && fHGasT.at(j)>time && fHGasT.at(j) < (time+MD->MyDelay)) fHeCal_E_tankAcap.at(copy) += fHGasE.at(j);
        }
        for(int tank=0; tank<NTanks; tank++) {
          if(fHeCal_E_tankAcap.at(tank)>HeCal_Acap_E){
            HeCal_Acap_E = fHeCal_E_tankAcap.at(tank);
            CopyNo = tank;
          }
        }
        analysisManager->FillNtupleIColumn(1, 50, CopyNo); // HeCal CopyNo for the MAX Acap energy release 
        analysisManager->FillNtupleDColumn(1, 51, HeCal_Acap_E);

      


          // capture delayed events

        if(NextCaptureTime==0){ // no double capture
          for(int j=0; j<NHG; j++){
            int copy = fHGasCopyNo.at(j);
            //getting the energy released in each tank
            if(fHGasE.at(j)>0. && fHGasT.at(j) > (time+MD->MyDelay) ) fHeCal_E_tankDelayed.at(copy) += fHGasE.at(j);
            if(fHGasE.at(j)>0. && fHGasT.at(j) > time_th ) fHeCal_E_tankLate.at(copy) += fHGasE.at(j);
          }
          for(int tank=0; tank<NTanks; tank++) {
            if(fHeCal_E_tankDelayed.at(tank)>HeCal_Delayed_E){
              HeCal_Delayed_E = fHeCal_E_tankDelayed.at(tank);
              CopyNo = tank;
            }
            if(fHeCal_E_tankLate.at(tank)>HeCal_Late_E){
              HeCal_Late_E = fHeCal_E_tankLate.at(tank);
              MyCopyNo = tank;
            }
          }
        }

        else{ // double capture
          for(int j=0; j<NHG; j++){
            int copy = fHGasCopyNo.at(j);
            if(fHGasE.at(j)>0. && (NextCaptureTime+MD->MyDelay) > fHGasT.at(j) && fHGasT.at(j) > (time+MD->MyDelay)) fHeCal_E_tankDelayed.at(copy) += fHGasE.at(j);
            if(fHGasE.at(j)>0. && fHGasT.at(j) > time_th ) fHeCal_E_tankLate.at(copy) += fHGasE.at(j);
          }
          for(int tank=0; tank<NTanks; tank++) {
            if(fHeCal_E_tankDelayed.at(tank)>HeCal_Delayed_E){
              HeCal_Delayed_E = fHeCal_E_tankDelayed.at(tank);
              CopyNo = tank;
            }
            if(fHeCal_E_tankLate.at(tank)>HeCal_Late_E){
              HeCal_Late_E = fHeCal_E_tankLate.at(tank);
              MyCopyNo = tank;
            }
          }
        }
        
        analysisManager->FillNtupleIColumn(1, 52, CopyNo); // HeCal CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 53, HeCal_Delayed_E);

        analysisManager->FillNtupleIColumn(1, 102, MyCopyNo); // HeCal CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 103, HeCal_Late_E);

      



        // ---- HODOSCOPES ---- //

          //prompt energy and hit Hod3
        for(int j=0; j<NHod3; j++){
          int copy3 = fHod3CopyNo.at(j);
          if(fHod3T.at(j)<time && fHod3E.at(j)!=0){
            fHod3_E_slabPrompt.at(copy3) += fHod3E.at(j);
            Hod3Hit_prompt++;
          }
          if(fHod3T.at(j)<time_th && fHod3E.at(j)!=0){
            fHod3_E_slabMyPrompt.at(copy3) += fHod3E.at(j);
            Hod3Hit_Myprompt++;
          }
        }


        for(int hod=0; hod<NScint; hod++){
          if(fHod3_E_slabPrompt.at(hod)>0 && fHod3_E_slabPrompt.at(hod)>Hod3_prompt_E){
            Hod3_prompt_E = fHod3_E_slabPrompt.at(hod);
            CopyNo = hod;
          }
          if(fHod3_E_slabMyPrompt.at(hod)>0 && fHod3_E_slabMyPrompt.at(hod)>Hod3_Myprompt_E){
            Hod3_Myprompt_E = fHod3_E_slabMyPrompt.at(hod);
            MyCopyNo = hod;
          }
        }
        analysisManager->FillNtupleIColumn(1, 63, CopyNo); // Hod3 CopyNo for the MAX Prompt energy release 
        analysisManager->FillNtupleDColumn(1, 64, Hod3_prompt_E);
        analysisManager->FillNtupleIColumn(1, 65, Hod3Hit_prompt);

        analysisManager->FillNtupleIColumn(1, 79, MyCopyNo); // Hod3 CopyNo for the MAX Prompt energy release 
        analysisManager->FillNtupleDColumn(1, 80, Hod3_Myprompt_E);
        analysisManager->FillNtupleIColumn(1, 81, Hod3Hit_Myprompt);

        // if(Hod3_Myprompt_E<0.6) printf("\n\nATTENTIZIONE MAX EDEP HOD3 < 0.6 \n Pari a: %f , copyNo = %d\n", Hod3_Myprompt_E, MyCopyNo);


          //prompt energy and hit Hod2
        for(int j=0; j<NHod2; j++){
          int copy2 = fHod2CopyNo.at(j);
          if(fHod2T.at(j)<time && fHod2E.at(j)!=0) {
            fHod2_E_slabPrompt.at(copy2) += fHod2E.at(j);
            Hod2Hit_prompt++;
          }
          if(fHod2T.at(j)<time_th && fHod2E.at(j)!=0) {
            fHod2_E_slabMyPrompt.at(copy2) += fHod2E.at(j);
            Hod2Hit_Myprompt++;
          }
        }
        for(int hod=0; hod<NScint; hod++){
          if(fHod2_E_slabPrompt.at(hod)>0 && fHod2_E_slabPrompt.at(hod)>Hod2_prompt_E){
            Hod2_prompt_E = fHod2_E_slabPrompt.at(hod);
            CopyNo = hod;
          }
          if(fHod2_E_slabMyPrompt.at(hod)>0 && fHod2_E_slabMyPrompt.at(hod)>Hod2_Myprompt_E){
            Hod2_Myprompt_E = fHod2_E_slabMyPrompt.at(hod);
            MyCopyNo = hod;
          }
        }
        analysisManager->FillNtupleIColumn(1, 54, CopyNo); // Hod2 CopyNo for the MAX Prompt energy release 
        analysisManager->FillNtupleDColumn(1, 55, Hod2_prompt_E);
        analysisManager->FillNtupleIColumn(1, 56, Hod2Hit_prompt);

        analysisManager->FillNtupleIColumn(1, 90, MyCopyNo); // Hod2 CopyNo for the MAX Prompt energy release 
        analysisManager->FillNtupleDColumn(1, 91, Hod2_Myprompt_E);
        analysisManager->FillNtupleIColumn(1, 92, Hod2Hit_Myprompt);

        // if(Hod2_Myprompt_E<0.6) printf("\n\nATTENTIZIONE MAX EDEP HOD2 < 0.6 \n Pari a: %f , copyNo = %d\n", Hod2_Myprompt_E, MyCopyNo);

        // if(Hod3_Myprompt_E<0.6 && Hod2_Myprompt_E<0.6) printf("\n\nATTENTIZIONE EVENTO INTERESSANTE!!! \n \t MAX EDEP HOD3 < 0.6 e MAX EDEP HOD2 < 0.6 \n Pari a: %f , copyNo = %d\n", Hod3_Myprompt_E, MyCopyNo);


      


          // after capture and before delayed events Hod3
        for(int j=0; j<NHod3; j++){
          int copy3 = fHod3CopyNo.at(j);
          if(fHod3T.at(j)>time && fHod3T.at(j)<(time+MD->MyDelay) && fHod3E.at(j)!=0){
            fHod3_E_slabAcap.at(copy3) += fHod3E.at(j);
            Hod3Hit_Acap++;
          }
        }
        for(int hod=0; hod<NScint; hod++){
          if(fHod3_E_slabAcap.at(hod)>0 && fHod3_E_slabAcap.at(hod)>Hod3_Acap_E) {
            Hod3_Acap_E = fHod3_E_slabAcap.at(hod);
            CopyNo = hod;
          }
        }
        analysisManager->FillNtupleIColumn(1, 66, CopyNo); // Hod3 CopyNo for the MAX Acap energy release 
        analysisManager->FillNtupleDColumn(1, 67, Hod3_Acap_E);
        analysisManager->FillNtupleIColumn(1, 68, Hod3Hit_Acap);

      


          // after capture and before delayed events Hod2
        for(int j=0; j<NHod2; j++){
          int copy2 = fHod2CopyNo.at(j);
          if(fHod2T.at(j)>time && fHod2T.at(j)<(time+MD->MyDelay) && fHod2E.at(j)!=0) {
            fHod2_E_slabAcap.at(copy2) += fHod2E.at(j);
            Hod2Hit_Acap++;
          }
        }
        for(int hod=0; hod<NScint; hod++){
            if(fHod2_E_slabAcap.at(hod)>0 && fHod2_E_slabAcap.at(hod)>Hod2_Acap_E) {
              Hod2_Acap_E = fHod2_E_slabAcap.at(hod);
              CopyNo = hod;
            }
        }
        analysisManager->FillNtupleIColumn(1, 57, CopyNo); // Hod2 CopyNo for the MAX Acap energy release 
        analysisManager->FillNtupleDColumn(1, 58, Hod2_Acap_E);
        analysisManager->FillNtupleIColumn(1, 59, Hod2Hit_Acap);



      

        

          // capture delayed events Hod3
        if(NextCaptureTime==0){ // no double capture
          for(int j=0; j<NHod3; j++){
            int copy3 = fHod3CopyNo.at(j);
            if(fHod3T.at(j)>(time+MD->MyDelay) && fHod3E.at(j)!=0){
              fHod3_E_slabDelayed.at(copy3) += fHod3E.at(j);
              Hod3Hit_Delayed++;
            }
            if(fHod3T.at(j)>time_th && fHod3E.at(j)!=0){
              fHod3_E_slabLate.at(copy3) += fHod3E.at(j);
              Hod3Hit_Late++;
            }
          }
          for(int hod=0; hod<NScint; hod++){
            if(fHod3_E_slabDelayed.at(hod)>0 && fHod3_E_slabDelayed.at(hod)>Hod3_Delayed_E) {
              Hod3_Delayed_E = fHod3_E_slabDelayed.at(hod);
              CopyNo = hod;
            }
            if(fHod3_E_slabLate.at(hod)>0 && fHod3_E_slabLate.at(hod)>Hod3_Late_E) {
              Hod3_Late_E = fHod3_E_slabLate.at(hod);
              MyCopyNo = hod;
            }
          }
        }

        else{// double capture
          for(int j=0; j<NHod3; j++){
            int copy3 = fHod3CopyNo.at(j);
            if( (NextCaptureTime+MD->MyDelay)>fHod3T.at(j) && fHod3T.at(j)>(time+MD->MyDelay) && fHod3E.at(j)!=0 ){
              fHod3_E_slabDelayed.at(copy3) += fHod3E.at(j);
              Hod3Hit_Delayed++;
            }
            if(fHod3T.at(j)>time_th && fHod3E.at(j)!=0){
              fHod3_E_slabLate.at(copy3) += fHod3E.at(j);
              Hod3Hit_Late++;
            }
          }
          for(int hod=0; hod<NScint; hod++){
            if(fHod3_E_slabDelayed.at(hod)>0 && fHod3_E_slabDelayed.at(hod)>Hod3_Delayed_E) {
              Hod3_Delayed_E = fHod3_E_slabDelayed.at(hod);
              CopyNo = hod;
            }
            if(fHod3_E_slabLate.at(hod)>0 && fHod3_E_slabLate.at(hod)>Hod3_Late_E) {
              Hod3_Late_E = fHod3_E_slabLate.at(hod);
              MyCopyNo = hod;
            }
          }
        }
        analysisManager->FillNtupleIColumn(1, 69, CopyNo); // Hod3 CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 70, Hod3_Delayed_E);
        analysisManager->FillNtupleIColumn(1, 71, Hod3Hit_Delayed);

        analysisManager->FillNtupleIColumn(1, 82, MyCopyNo); // Hod3 CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 83, Hod3_Late_E);
        analysisManager->FillNtupleIColumn(1, 84, Hod3Hit_Late);



      

      


          // capture delayed events Hod2
        if(NextCaptureTime==0){ // no double capture
          for(int j=0; j<NHod2; j++){
            int copy2 = fHod2CopyNo.at(j);
            if(fHod2T.at(j)>(time+MD->MyDelay) && fHod2E.at(j)!=0) {
              fHod2_E_slabDelayed.at(copy2) += fHod2E.at(j);
              Hod2Hit_Delayed++;
            }
            if(fHod2T.at(j)>time_th && fHod2E.at(j)!=0) {
              fHod2_E_slabLate.at(copy2) += fHod2E.at(j);
              Hod2Hit_Late++;
            }
          }

          for(int hod=0; hod<NScint; hod++){
            if(fHod2_E_slabDelayed.at(hod)>0 && fHod2_E_slabDelayed.at(hod)>Hod2_Delayed_E) {
              Hod2_Delayed_E = fHod2_E_slabDelayed.at(hod);
              CopyNo = hod;
            }
            if(fHod2_E_slabLate.at(hod)>0 && fHod2_E_slabLate.at(hod)>Hod2_Late_E) {
              Hod2_Late_E = fHod2_E_slabLate.at(hod);
              MyCopyNo = hod;
            }
          }
        }

        else{ // double capture
          for(int j=0; j<NHod2; j++){
            int copy2 = fHod2CopyNo.at(j);
            if( (NextCaptureTime+MD->MyDelay)>fHod2T.at(j) && fHod2T.at(j)>(time+MD->MyDelay) && fHod2E.at(j)!=0) {
              fHod2_E_slabDelayed.at(copy2) += fHod2E.at(j);
              Hod2Hit_Delayed++;
            }
            if(fHod2T.at(j)>time_th && fHod2E.at(j)!=0) {
              fHod2_E_slabLate.at(copy2) += fHod2E.at(j);
              Hod2Hit_Late++;
            }
          }

          for(int hod=0; hod<NScint; hod++){
            if(fHod2_E_slabDelayed.at(hod)>0 && fHod2_E_slabDelayed.at(hod)>Hod2_Delayed_E) {
              Hod2_Delayed_E = fHod2_E_slabDelayed.at(hod);
              CopyNo = hod;
            }
            if(fHod2_E_slabLate.at(hod)>0 && fHod2_E_slabLate.at(hod)>Hod2_Late_E) {
              Hod2_Late_E = fHod2_E_slabLate.at(hod);
              MyCopyNo = hod;
            }
          }
        }
        
        analysisManager->FillNtupleIColumn(1, 60, CopyNo); // Hod2 CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 61, Hod2_Delayed_E);
        analysisManager->FillNtupleIColumn(1, 62, Hod2Hit_Delayed);

        analysisManager->FillNtupleIColumn(1, 93, MyCopyNo); // Hod2 CopyNo for the MAX Delayed energy release 
        analysisManager->FillNtupleDColumn(1, 94, Hod2_Late_E);
        analysisManager->FillNtupleIColumn(1, 95, Hod2Hit_Late);



        


        // --- Beta evaluatuion --- //

        if(NHod3>0 && NHod2>0){
          G4double dist = sqrt(pow((fHod2Y.at(0)-fHod3Y.at(0)),2) + pow((fHod2X.at(0)-fHod3X.at(0)),2) + pow((fHod2Z.at(0)-fHod3Z.at(0)),2));
          beta = (1./299.792458)*dist/(fHod2T.at(0)-fHod3T.at(0));  
        }
      
        // evaluate incidence on ToF shell (F.Nozzoli)
        // Cosi[0]=1;
        // Cosi[1]=1;
        // Cosi[2]=1; // cosine WRT outher shell
        // double dir[3];
        // if(NHod3>0 && NHod2>0){
        //   dir[0] = fHod2X.at(0)-fHod3X.at(0);
        //   dir[1] = fHod2Y.at(0)-fHod3Y.at(0);
        //   dir[2] = fHod2Z.at(0)-fHod3Z.at(0); 
        // }
        // else {dir[0]=1;dir[1]=1;dir[2]=1;}
        // double len = sqrt(pow(dir[0],2) + pow(dir[1],2) + pow(dir[2],2)); 
        // for(int jj = 0;jj<3;jj++) dir[jj]=dir[jj]/len;
        // if(NHod3>0){
        //   if(fabs(fabs(fHod3Z.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[2]); // select SIDE
        //   if(fabs(fabs(fHod3Y.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[1]);
        //   if(fabs(fabs(fHod3X.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[0]);
        // }
        // if(NHod2>0){
        //   if(fabs(fabs(fHod2Z.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[2]); // select SIDE
        //   if(fabs(fabs(fHod2Y.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[1]);
        //   if(fabs(fabs(fHod2X.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[0]);
        // }

        analysisManager->FillNtupleDColumn(1, 72, beta);
        // analysisManager->FillNtupleDColumn(1, Cosi);

       

        analysisManager->AddNtupleRow(1);

        
        break; // ONLY THE FIRST CAPTURE

      } // end of the cycle on the first captures
    }


    
    else{
      // >>> NORMAL PARTICLES <<< //
      // initialization
      G4double HeCal_prompt_E=0, HeCal_Late_E=0;
      G4double Hod3_prompt_E=0, Hod3_Late_E=0;
      G4double Hod2_prompt_E=0, Hod2_Late_E=0;
      G4double beta=0;
      G4int Hod2Hit_prompt=0, Hod2Hit_Late=0;
      G4int Hod3Hit_prompt=0, Hod3Hit_Late=0;
      int CopyNo = 999;

      const G4double time = 50 + fHod3T.at(0); // ns

      // Incidence on ToF shell (F.Nozzoli)
      // G4float Cosi[3];



      // ---- HECAL ---- //

        //prompt energy
      for(int j=0; j<NHG; j++){
          int copy = fHGasCopyNo.at(j);
          //getting the energy released in each tank
          if(fHGasE.at(j)>0. && fHGasT.at(j)<time) fHeCal_E_tankPrompt.at(copy) += fHGasE.at(j);
      }

      for(int tank=0; tank<NTanks; tank++) {
        if(fHeCal_E_tankPrompt.at(tank)>HeCal_prompt_E) {
          HeCal_prompt_E = fHeCal_E_tankPrompt.at(tank);
          CopyNo = tank;
        }
      }
      analysisManager->FillNtupleIColumn(2, 48, CopyNo); // HeCal CopyNo for the MAX prompt energy release 


        //late events
      for(int j=0; j<NHG; j++){
        int copy = fHGasCopyNo.at(j);
        if(fHGasE.at(j)>0. && fHGasT.at(j)>time) fHeCal_E_tankDelayed.at(copy) += fHGasE.at(j);
      }

      for(int tank=0; tank<NTanks; tank++) {
        if(fHeCal_E_tankDelayed.at(tank)>HeCal_Late_E){
          HeCal_Late_E = fHeCal_E_tankDelayed.at(tank) ;
          CopyNo = tank;
        }
      }
      analysisManager->FillNtupleIColumn(2, 50, CopyNo); // HeCal CopyNo for the MAX late energy release 



      // ---- HODOSCOPES ---- //

        //prompt energy and hit Hod3
      for(int j=0; j<NHod3; j++){
        int copy3 = fHod3CopyNo.at(j);
        if(fHod3T.at(j)<time && fHod3E.at(j)!=0){
          fHod3_E_slabPrompt.at(copy3) += fHod3E.at(j);
          Hod3Hit_prompt++;
        }
      }

      for(int hod=0; hod<NScint; hod++){
        if(fHod3_E_slabPrompt.at(hod)>0 && fHod3_E_slabPrompt.at(hod)>Hod3_prompt_E){
          Hod3_prompt_E = fHod3_E_slabPrompt.at(hod);
          CopyNo = hod;
        }
      }
      analysisManager->FillNtupleIColumn(2, 58, CopyNo); // Hod3 CopyNo for the MAX Prompt energy release 


        //prompt energy and hit Hod2
      for(int j=0; j<NHod2; j++){
        int copy2 = fHod2CopyNo.at(j);
        if(fHod2T.at(j)<time && fHod2E.at(j)!=0) {
          fHod2_E_slabPrompt.at(copy2) += fHod2E.at(j);
          Hod2Hit_prompt++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(fHod2_E_slabPrompt.at(hod)>0 && fHod2_E_slabPrompt.at(hod)>Hod2_prompt_E){
          Hod2_prompt_E = fHod2_E_slabPrompt.at(hod);
          CopyNo = hod;
        }
      }
      analysisManager->FillNtupleIColumn(2, 52, CopyNo); // Hod2 CopyNo for the MAX Prompt energy release 



        // late events Hod3
      for(int j=0; j<NHod3; j++){
        int copy3 = fHod3CopyNo.at(j);
        if(fHod3T.at(j)>time && fHod3E.at(j)!=0){
          fHod3_E_slabDelayed.at(copy3) += fHod3E.at(j);
          Hod3Hit_Late++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(fHod3_E_slabDelayed.at(hod)>0 && fHod3_E_slabDelayed.at(hod)>Hod3_Late_E) {
          Hod3_Late_E = fHod3_E_slabDelayed.at(hod);
          CopyNo = hod;
        }
      }
      analysisManager->FillNtupleIColumn(2, 61, CopyNo); // Hod3 CopyNo for the MAX Late energy release 


        // late events Hod2
      for(int j=0; j<NHod2; j++){
        int copy2 = fHod2CopyNo.at(j);
        if(fHod2T.at(j)>time && fHod2E.at(j)!=0) {
          fHod2_E_slabDelayed.at(copy2) += fHod2E.at(j);
          Hod2Hit_Late++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(fHod2_E_slabDelayed.at(hod)>0 && fHod2_E_slabDelayed.at(hod)>Hod2_Late_E) {
          Hod2_Late_E = fHod2_E_slabDelayed.at(hod);
          CopyNo = hod;
        }
      }
      analysisManager->FillNtupleIColumn(2, 55, CopyNo); // Hod2 CopyNo for the MAX Late energy release 


      // --- Beta evaluatuion --- //

      if(NHod3>0 && NHod2>0){
        G4double dist = sqrt(pow((fHod2Y.at(0)-fHod3Y.at(0)),2) + pow((fHod2X.at(0)-fHod3X.at(0)),2) + pow((fHod2Z.at(0)-fHod3Z.at(0)),2));
        beta = (1./299.792458)*dist/(fHod2T.at(0)-fHod3T.at(0));  
      }
    
      // evaluate incidence on ToF shell (F.Nozzoli)
      // Cosi[0]=1;
      // Cosi[1]=1;
      // Cosi[2]=1; // cosine WRT outher shell
      // double dir[3];
      // if(NHod3>0 && NHod2>0){
      //   dir[0] = fHod2X.at(0)-fHod3X.at(0);
      //   dir[1] = fHod2Y.at(0)-fHod3Y.at(0);
      //   dir[2] = fHod2Z.at(0)-fHod3Z.at(0); 
      // }
      // else {dir[0]=1;dir[1]=1;dir[2]=1;}
      // double len = sqrt(pow(dir[0],2) + pow(dir[1],2) + pow(dir[2],2)); 
      // for(int jj = 0;jj<3;jj++) dir[jj]=dir[jj]/len;
      // if(NHod3>0){
      //   if(fabs(fabs(fHod3Z.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[2]); // select SIDE
      //   if(fabs(fabs(fHod3Y.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[1]);
      //   if(fabs(fabs(fHod3X.at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[0]);
      // }
      // if(NHod2>0){
      //   if(fabs(fabs(fHod2Z.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[2]); // select SIDE
      //   if(fabs(fabs(fHod2Y.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[1]);
      //   if(fabs(fabs(fHod2X.at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[0]);
      // }


      // ---- FILLING NTUPLE ---- //

      analysisManager->FillNtupleDColumn(2, 0, MCEnergy/NofNucleons);
      analysisManager->FillNtupleDColumn(2, 1, MCMomentum/NofNucleons);
      analysisManager->FillNtupleDColumn(2, 9, float(MD->MCEvnum)/MD->GenSurface);
      analysisManager->FillNtupleIColumn(2, 10,NHod2);
      analysisManager->FillNtupleIColumn(2, 11,NHod3);
      analysisManager->FillNtupleIColumn(2, 12,NHG);
      analysisManager->FillNtupleDColumn(2, 13, (primary->GetMomentum())[0]/NofNucleons);
      analysisManager->FillNtupleDColumn(2, 14, (primary->GetMomentum())[1]/NofNucleons);
      analysisManager->FillNtupleDColumn(2, 15, (primary->GetMomentum())[2]/NofNucleons);
      analysisManager->FillNtupleDColumn(2, 16, (priv->GetPosition())[0]);
      analysisManager->FillNtupleDColumn(2, 17, (priv->GetPosition())[1]);
      analysisManager->FillNtupleDColumn(2, 18, (priv->GetPosition())[2]);
      analysisManager->FillNtupleDColumn(2, 35, MCMass);
      analysisManager->FillNtupleDColumn(2, 36, MCCharge);
      analysisManager->FillNtupleIColumn(2, 45, event->GetEventID());

      analysisManager->FillNtupleDColumn(2, 49, HeCal_prompt_E);
      analysisManager->FillNtupleDColumn(2, 51, HeCal_Late_E);

      analysisManager->FillNtupleDColumn(2, 53, Hod2_prompt_E);
      analysisManager->FillNtupleIColumn(2, 54, Hod2Hit_prompt);
      analysisManager->FillNtupleDColumn(2, 56, Hod2_Late_E);
      analysisManager->FillNtupleIColumn(2, 57, Hod2Hit_Late);

      analysisManager->FillNtupleDColumn(2, 59, Hod3_prompt_E);
      analysisManager->FillNtupleIColumn(2, 60, Hod3Hit_prompt);
      analysisManager->FillNtupleDColumn(2, 62, Hod3_Late_E);
      analysisManager->FillNtupleIColumn(2, 63, Hod3Hit_Late);

      analysisManager->FillNtupleDColumn(2, 64, beta);
      // analysisManager->FillNtupleDColumn(2, Cosi);

      analysisManager->AddNtupleRow(2); 

    }




    // reclean needed
    fHod2T.clear();
    fHod2X.clear();
    fHod2Y.clear();
    fHod2Z.clear();
    fHod2E.clear();
    fHod2CopyNo.clear();

    fHod3T.clear();
    fHod3X.clear();
    fHod3Y.clear();
    fHod3Z.clear();
    fHod3E.clear();
    fHod3CopyNo.clear();

    fHGasT.clear();
    fHGasX.clear();
    fHGasY.clear();
    fHGasZ.clear();
    fHGasE.clear();
    fHGasCopyNo.clear();

    fHod3_E_slabPrompt.clear();
    fHod3_E_slabAcap.clear();
    fHod3_E_slabDelayed.clear();
    fHod3_E_slabMyPrompt.clear();
    fHod3_E_slabLate.clear();

    fHod2_E_slabPrompt.clear();
    fHod2_E_slabAcap.clear();
    fHod2_E_slabDelayed.clear();
    fHod2_E_slabMyPrompt.clear();
    fHod2_E_slabLate.clear();

    fHeCal_E_tankPrompt.clear();
    fHeCal_E_tankAcap.clear();
    fHeCal_E_tankDelayed.clear();
    fHeCal_E_tankMyPrompt.clear();
    fHeCal_E_tankLate.clear();
    

    MD->MyCaptureT.clear();
    MD->MyCaptureX.clear();
    MD->MyCaptureY.clear();
    MD->MyCaptureZ.clear();
    MD->MyCaptureE.clear();
    MD->MyCaptureF.clear();
    MD->MyCaptureM.clear();
    MD->MyCaptureCopyNo.clear();

    MD->HodCaptureT.clear();
    MD->HodCaptureX.clear();
    MD->HodCaptureY.clear();
    MD->HodCaptureZ.clear();
    MD->HodCaptureE.clear();
    MD->HodCaptureF.clear();
    MD->HodCaptureM.clear();
    MD->HodCaptureVolume.clear();

  }
  

  // Print diagnostics
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();

  if (printModulo==0 || event->GetEventID() % printModulo != 0) return;


  
  G4cout << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " PDG code " << primary->GetPDGcode()
    << " momentum vector " << primary->GetMomentum()
	  << " position " << event->GetPrimaryVertex(0)->GetPosition()
	<< G4endl;
    
  G4cout << MD->MCEvnum << " Generated " << MD->TRGEvnum << " triggered SurF=" <<MD->GenSurface << " Nucleons number = " << NofNucleons << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
