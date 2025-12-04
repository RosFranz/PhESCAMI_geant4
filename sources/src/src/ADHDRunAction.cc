// $Id: ADHDRunAction.cc 74204 2013-10-01 07:04:43Z ihrivnac $
//
/// \file ADHDRunAction.cc
/// \brief Implementation of the ADHDRunAction class

#include "ADHDRunAction.hh"
#include "ADHDEventAction.hh"


#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "MyDet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDRunAction::ADHDRunAction(ADHDEventAction* eventAction): fEventAction(eventAction)
{ 
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetActivation(true); // allows to de-/activate histograms and ntuples
  analysisManager->SetDefaultFileType("root");

  MyDet *MD = MyDet::gethead();

  // Default settings
  analysisManager->SetVerboseLevel(2);
  analysisManager->SetFileName("Events");
  G4int nbins = 400;
  G4double Ek_min = 10, Ek_max = 100000; // MeV
  // G4cout << "index = " << index << G4endl;

  // if(MD->antiparticle == 1) analysisManager->SetFileName("Events_bar");
  analysisManager->SetNtupleFileName(0,"ID0");

  // if(MD->antiparticle == 1) analysisManager->SetNtupleFileName(1,"ID1");


  // --------------------------------------------------
  
  // Creating ntuple
  if ( fEventAction ) {

    G4int index = analysisManager->CreateH1("Gen_spectra", "Generated spectra", nbins, Ek_min, Ek_max, "none", "none", "log");
    if(index == -1) G4cout << "Problem creating histogram" << G4endl;


    // >>> CREATING NTUPLES FOR ANTIPARTICLES SIMULATION <<< //

    // --- FIRST Ntuple  ID = 0  with all hits (raw data) ---
    analysisManager->CreateNtuple("Hits_bar", "Hits_bar");

    analysisManager->CreateNtupleDColumn(0, "MCEnergy"); // column Id = 0
    analysisManager->CreateNtupleDColumn(0, "MCMomentum"); // column Id = 1
    analysisManager->CreateNtupleDColumn(0, "CaptureT", MD->MyCaptureT); // column Id = 2
    analysisManager->CreateNtupleDColumn(0, "CaptureX", MD->MyCaptureX); // column Id = 3
    analysisManager->CreateNtupleDColumn(0, "CaptureY", MD->MyCaptureY); // column Id = 4
    analysisManager->CreateNtupleDColumn(0, "CaptureZ", MD->MyCaptureZ); // column Id = 5
    analysisManager->CreateNtupleDColumn(0, "CaptureE", MD->MyCaptureE); // column Id = 6
    analysisManager->CreateNtupleDColumn(0, "CaptureF", MD->MyCaptureF); // column Id = 7
    analysisManager->CreateNtupleDColumn(0, "CaptureM", MD->MyCaptureM); // column Id = 8
    analysisManager->CreateNtupleDColumn(0, "MCEvnum_over_GenSurface");  // column Id = 9

    analysisManager->CreateNtupleIColumn(0, "Hod2N");  // column Id = 10
    analysisManager->CreateNtupleIColumn(0, "Hod3N");  // column Id = 11
    analysisManager->CreateNtupleIColumn(0, "HGasN");  // column Id = 12

    analysisManager->CreateNtupleDColumn(0, "MCMomX"); // column Id = 13
    analysisManager->CreateNtupleDColumn(0, "MCMomY"); // column Id = 14
    analysisManager->CreateNtupleDColumn(0, "MCMomZ"); // column Id = 15
    analysisManager->CreateNtupleDColumn(0, "MCPosX"); // column Id = 16
    analysisManager->CreateNtupleDColumn(0, "MCPosY"); // column Id = 17
    analysisManager->CreateNtupleDColumn(0, "MCPosZ"); // column Id = 18

    // vector columns are automatically filled
    analysisManager->CreateNtupleDColumn(0, "Hod2E", fEventAction->GetHod2E()); // column Id = 19
    analysisManager->CreateNtupleDColumn(0, "Hod2T", fEventAction->GetHod2T()); // column Id = 20
    analysisManager->CreateNtupleDColumn(0, "Hod2X", fEventAction->GetHod2X()); // column Id = 21
    analysisManager->CreateNtupleDColumn(0, "Hod2Y", fEventAction->GetHod2Y()); // column Id = 22
    analysisManager->CreateNtupleDColumn(0, "Hod2Z", fEventAction->GetHod2Z()); // column Id = 23

    analysisManager->CreateNtupleDColumn(0, "Hod3E", fEventAction->GetHod3E()); // column Id = 24
    analysisManager->CreateNtupleDColumn(0, "Hod3T", fEventAction->GetHod3T()); // column Id = 25
    analysisManager->CreateNtupleDColumn(0, "Hod3X", fEventAction->GetHod3X()); // column Id = 26
    analysisManager->CreateNtupleDColumn(0, "Hod3Y", fEventAction->GetHod3Y()); // column Id = 27
    analysisManager->CreateNtupleDColumn(0, "Hod3Z", fEventAction->GetHod3Z()); // column Id = 28

    analysisManager->CreateNtupleDColumn(0, "HGasT", fEventAction->GetHGasT()); // column Id = 29
    analysisManager->CreateNtupleDColumn(0, "HGasX", fEventAction->GetHGasX()); // column Id = 30
    analysisManager->CreateNtupleDColumn(0, "HGasY", fEventAction->GetHGasY()); // column Id = 31
    analysisManager->CreateNtupleDColumn(0, "HGasZ", fEventAction->GetHGasZ()); // column Id = 32
    analysisManager->CreateNtupleDColumn(0, "HGasE", fEventAction->GetHGasE()); // column Id = 33
    analysisManager->CreateNtupleIColumn(0, "HGasCopyNo", fEventAction->GetHGasCopyNo()); // column Id = 34

    analysisManager->CreateNtupleDColumn(0, "MCMass"); // column Id = 35
    analysisManager->CreateNtupleDColumn(0, "MCCharge"); // column Id = 36

    analysisManager->CreateNtupleDColumn(0, "HodCaptureT", MD->HodCaptureT); // column Id = 37
    analysisManager->CreateNtupleDColumn(0, "HodCaptureX", MD->HodCaptureX); // column Id = 38
    analysisManager->CreateNtupleDColumn(0, "HodCaptureY", MD->HodCaptureY); // column Id = 39
    analysisManager->CreateNtupleDColumn(0, "HodCaptureZ", MD->HodCaptureZ); // column Id = 40
    analysisManager->CreateNtupleDColumn(0, "HodCaptureE", MD->HodCaptureE); // column Id = 41
    analysisManager->CreateNtupleDColumn(0, "HodCaptureF", MD->HodCaptureF); // column Id = 42
    analysisManager->CreateNtupleDColumn(0, "HodCaptureM", MD->HodCaptureM); // column Id = 43
    analysisManager->CreateNtupleIColumn(0, "HodCaptureVolume", MD->HodCaptureVolume); // column Id = 44
    analysisManager->CreateNtupleIColumn(0, "EventID"); // column Id = 45
    analysisManager->CreateNtupleIColumn(0, "Hod2CopyNo", fEventAction->GetHod2CopyNo()); // column Id = 46
    analysisManager->CreateNtupleIColumn(0, "Hod3CopyNo", fEventAction->GetHod3CopyNo()); // column Id = 47
    analysisManager->CreateNtupleIColumn(0, "CaptureCopyNo", MD->MyCaptureCopyNo); // column Id = 48"
    
    analysisManager->FinishNtuple(0);

    

    // --- SECOND Ntuple ID = 1 reconstructed data for antiparticles (only events with capture in HeCal) ---
    analysisManager->CreateNtuple("RHits_bar", "RHits_bar");

    analysisManager->CreateNtupleDColumn(1, "MCEnergy"); // column Id = 0
    analysisManager->CreateNtupleDColumn(1, "MCMomentum"); // column Id = 1
    analysisManager->CreateNtupleDColumn(1, "CaptureT", MD->MyCaptureT); // column Id = 2
    analysisManager->CreateNtupleDColumn(1, "CaptureX", MD->MyCaptureX); // column Id = 3
    analysisManager->CreateNtupleDColumn(1, "CaptureY", MD->MyCaptureY); // column Id = 4
    analysisManager->CreateNtupleDColumn(1, "CaptureZ", MD->MyCaptureZ); // column Id = 5
    analysisManager->CreateNtupleDColumn(1, "CaptureE", MD->MyCaptureE); // column Id = 6
    analysisManager->CreateNtupleDColumn(1, "CaptureF", MD->MyCaptureF); // column Id = 7
    analysisManager->CreateNtupleDColumn(1, "CaptureM", MD->MyCaptureM); // column Id = 8
    analysisManager->CreateNtupleDColumn(1, "MCEvnum_over_GenSurface");  // column Id = 9

    analysisManager->CreateNtupleIColumn(1, "Hod2N");  // column Id = 10
    analysisManager->CreateNtupleIColumn(1, "Hod3N");  // column Id = 11
    analysisManager->CreateNtupleIColumn(1, "HGasN");  // column Id = 12

    analysisManager->CreateNtupleDColumn(1, "MCMomX"); // column Id = 13
    analysisManager->CreateNtupleDColumn(1, "MCMomY"); // column Id = 14
    analysisManager->CreateNtupleDColumn(1, "MCMomZ"); // column Id = 15
    analysisManager->CreateNtupleDColumn(1, "MCPosX"); // column Id = 16
    analysisManager->CreateNtupleDColumn(1, "MCPosY"); // column Id = 17
    analysisManager->CreateNtupleDColumn(1, "MCPosZ"); // column Id = 18

    // vector columns are automatically filled
    analysisManager->CreateNtupleDColumn(1, "Hod2E", fEventAction->GetHod2E()); // column Id = 19
    analysisManager->CreateNtupleDColumn(1, "Hod2T", fEventAction->GetHod2T()); // column Id = 20
    analysisManager->CreateNtupleDColumn(1, "Hod2X", fEventAction->GetHod2X()); // column Id = 21
    analysisManager->CreateNtupleDColumn(1, "Hod2Y", fEventAction->GetHod2Y()); // column Id = 22
    analysisManager->CreateNtupleDColumn(1, "Hod2Z", fEventAction->GetHod2Z()); // column Id = 23

    analysisManager->CreateNtupleDColumn(1, "Hod3E", fEventAction->GetHod3E()); // column Id = 24
    analysisManager->CreateNtupleDColumn(1, "Hod3T", fEventAction->GetHod3T()); // column Id = 25
    analysisManager->CreateNtupleDColumn(1, "Hod3X", fEventAction->GetHod3X()); // column Id = 26
    analysisManager->CreateNtupleDColumn(1, "Hod3Y", fEventAction->GetHod3Y()); // column Id = 27
    analysisManager->CreateNtupleDColumn(1, "Hod3Z", fEventAction->GetHod3Z()); // column Id = 28

    analysisManager->CreateNtupleDColumn(1, "HGasT", fEventAction->GetHGasT()); // column Id = 29
    analysisManager->CreateNtupleDColumn(1, "HGasX", fEventAction->GetHGasX()); // column Id = 30
    analysisManager->CreateNtupleDColumn(1, "HGasY", fEventAction->GetHGasY()); // column Id = 31
    analysisManager->CreateNtupleDColumn(1, "HGasZ", fEventAction->GetHGasZ()); // column Id = 32
    analysisManager->CreateNtupleDColumn(1, "HGasE", fEventAction->GetHGasE()); // column Id = 33
    analysisManager->CreateNtupleIColumn(1, "HGasCopyNo", fEventAction->GetHGasCopyNo()); // column Id = 34

    analysisManager->CreateNtupleDColumn(1, "MCMass"); // column Id = 35
    analysisManager->CreateNtupleDColumn(1, "MCCharge"); // column Id = 36

    analysisManager->CreateNtupleDColumn(1, "HodCaptureT", MD->HodCaptureT); // column Id = 37
    analysisManager->CreateNtupleDColumn(1, "HodCaptureX", MD->HodCaptureX); // column Id = 38
    analysisManager->CreateNtupleDColumn(1, "HodCaptureY", MD->HodCaptureY); // column Id = 39
    analysisManager->CreateNtupleDColumn(1, "HodCaptureZ", MD->HodCaptureZ); // column Id = 40
    analysisManager->CreateNtupleDColumn(1, "HodCaptureE", MD->HodCaptureE); // column Id = 41
    analysisManager->CreateNtupleDColumn(1, "HodCaptureF", MD->HodCaptureF); // column Id = 42
    analysisManager->CreateNtupleDColumn(1, "HodCaptureM", MD->HodCaptureM); // column Id = 43
    analysisManager->CreateNtupleIColumn(1, "HodCaptureVolume", MD->HodCaptureVolume); // column Id = 44
    analysisManager->CreateNtupleIColumn(1, "EventID"); // column Id = 45
    analysisManager->CreateNtupleIColumn(1, "Hod2CopyNo", fEventAction->GetHod2CopyNo()); // column Id = 46
    analysisManager->CreateNtupleIColumn(1, "Hod3CopyNo", fEventAction->GetHod3CopyNo()); // column Id = 47


    //ADHD detector reconstruction

    //HeCal reconstruction for captures
    analysisManager->CreateNtupleIColumn(1, "HeCal_MaxEprompt_CopyNo"); // column Id = 48
    analysisManager->CreateNtupleDColumn(1, "HeCal_MaxEprompt"); // column Id = 49

    analysisManager->CreateNtupleIColumn(1, "HeCal_MaxEAcap_CopyNo"); // column Id = 50
    analysisManager->CreateNtupleDColumn(1, "HeCal_MaxEAcap"); // column Id = 51

    analysisManager->CreateNtupleIColumn(1, "HeCal_MaxEDelayed_CopyNo"); // column Id = 52
    analysisManager->CreateNtupleDColumn(1, "HeCal_MaxEDelayed"); // column Id = 53

    //Hod2 reconstruction for captures
    analysisManager->CreateNtupleIColumn(1, "Hod2_MaxEprompt_CopyNo"); // column Id = 54
    analysisManager->CreateNtupleDColumn(1, "Hod2_MaxEprompt"); // column Id = 55
    analysisManager->CreateNtupleIColumn(1, "Hod2Hit_prompt"); // column Id = 56

    analysisManager->CreateNtupleIColumn(1, "Hod2_MaxEAcap_CopyNo"); // column Id = 57
    analysisManager->CreateNtupleDColumn(1, "Hod2_MaxEAcap"); // column Id = 58
    analysisManager->CreateNtupleIColumn(1, "Hod2Hit_Acap"); // column Id = 59

    analysisManager->CreateNtupleIColumn(1, "Hod2_MaxEDelayed_CopyNo"); // column Id = 60
    analysisManager->CreateNtupleDColumn(1, "Hod2_MaxEDelayed"); // column Id = 61
    analysisManager->CreateNtupleIColumn(1, "Hod2Hit_Delayed"); // column Id = 62

    //Hod3 reconstruction for captures
    analysisManager->CreateNtupleIColumn(1, "Hod3_MaxEprompt_CopyNo"); // column Id = 63
    analysisManager->CreateNtupleDColumn(1, "Hod3_MaxEprompt"); // column Id = 64
    analysisManager->CreateNtupleIColumn(1, "Hod3Hit_prompt"); // column Id = 65

    analysisManager->CreateNtupleIColumn(1, "Hod3_MaxEAcap_CopyNo"); // column Id = 66
    analysisManager->CreateNtupleDColumn(1, "Hod3_MaxEAcap"); // column Id = 67
    analysisManager->CreateNtupleIColumn(1, "Hod3Hit_Acap"); // column Id = 68

    analysisManager->CreateNtupleIColumn(1, "Hod3_MaxEDelayed_CopyNo"); // column Id = 69
    analysisManager->CreateNtupleDColumn(1, "Hod3_MaxEDelayed"); // column Id = 70
    analysisManager->CreateNtupleIColumn(1, "Hod3Hit_Delayed"); // column Id = 71

    //beta
    analysisManager->CreateNtupleDColumn(1, "beta"); // column Id = 72

    //HeCal capture copy noumber
    analysisManager->CreateNtupleIColumn(1, "CaptureCopyNo", MD->MyCaptureCopyNo); // column Id = 73



    //Energy released for each Hod3 slab for captures
    analysisManager->CreateNtupleDColumn(1, "Hod3_E_slabPrompt", fEventAction->GetHod3_E_slabPrompt() ); // column Id 74
    analysisManager->CreateNtupleDColumn(1, "Hod3_E_slabAcap", fEventAction->GetHod3_E_slabAcap() ); // column Id 75
    analysisManager->CreateNtupleDColumn(1, "Hod3_E_slabDelayed", fEventAction->GetHod3_E_slabDelayed() ); // column Id 76

    //Energy released for each Hod3 slab WITHOUT captures
    analysisManager->CreateNtupleDColumn(1, "Hod3_E_slabMyPrompt", fEventAction->GetHod3_E_slabMyPrompt() ); // column Id 77
    analysisManager->CreateNtupleDColumn(1, "Hod3_E_slabLate", fEventAction->GetHod3_E_slabLate() ); // column Id 78

    //Hod3 reconstruction WITHOUT captures
    analysisManager->CreateNtupleIColumn(1, "Hod3_MaxEMyprompt_CopyNo"); // column Id = 79
    analysisManager->CreateNtupleDColumn(1, "Hod3_MaxEMyprompt"); // column Id = 80
    analysisManager->CreateNtupleIColumn(1, "Hod3Hit_Myprompt"); // column Id = 81

    analysisManager->CreateNtupleIColumn(1, "Hod3_MaxELate_CopyNo"); // column Id = 82
    analysisManager->CreateNtupleDColumn(1, "Hod3_MaxELate"); // column Id = 83
    analysisManager->CreateNtupleIColumn(1, "Hod3Hit_Late"); // column Id = 84



    // //Energy released for each Hod2 slab for captures
    analysisManager->CreateNtupleDColumn(1, "Hod2_E_slabPrompt", fEventAction->GetHod2_E_slabPrompt() ); // column Id 85
    analysisManager->CreateNtupleDColumn(1, "Hod2_E_slabAcap", fEventAction->GetHod2_E_slabAcap() ); // column Id 86
    analysisManager->CreateNtupleDColumn(1, "Hod2_E_slabDelayed", fEventAction->GetHod2_E_slabDelayed() ); // column Id 87

    //Energy released for each Hod2 slab WITHOUT captures
    analysisManager->CreateNtupleDColumn(1, "Hod2_E_slabMyPrompt", fEventAction->GetHod2_E_slabMyPrompt() ); // column Id 88
    analysisManager->CreateNtupleDColumn(1, "Hod2_E_slabLate", fEventAction->GetHod2_E_slabLate() ); // column Id 89

    //Hod2 reconstruction without captures
    analysisManager->CreateNtupleIColumn(1, "Hod2_MaxEMyprompt_CopyNo"); // column Id = 90
    analysisManager->CreateNtupleDColumn(1, "Hod2_MaxEMyprompt"); // column Id = 91
    analysisManager->CreateNtupleIColumn(1, "Hod2Hit_Myprompt"); // column Id = 92

    analysisManager->CreateNtupleIColumn(1, "Hod2_MaxELate_CopyNo"); // column Id = 93
    analysisManager->CreateNtupleDColumn(1, "Hod2_MaxELate"); // column Id = 94
    analysisManager->CreateNtupleIColumn(1, "Hod2Hit_Late"); // column Id = 95



    // //Energy released for each HeCal tank
    analysisManager->CreateNtupleDColumn(1, "HeCal_E_tankPrompt", fEventAction->GetHeCal_E_tankPrompt() ); // column Id 96
    analysisManager->CreateNtupleDColumn(1, "HeCal_E_tankDelayed", fEventAction->GetHeCal_E_tankDelayed() ); // column Id 97
    
    //Energy released for each HeCal tank WITHOUT captures
    analysisManager->CreateNtupleDColumn(1, "HeCal_E_tankMyPrompt", fEventAction->GetHeCal_E_tankMyPrompt() ); // column Id 98
    analysisManager->CreateNtupleDColumn(1, "HeCal_E_tankLate", fEventAction->GetHeCal_E_tankLate() ); // column Id 99

    //HeCal reconstruction without captures
    analysisManager->CreateNtupleIColumn(1, "HeCal_MaxEMyprompt_CopyNo"); // column Id = 100
    analysisManager->CreateNtupleDColumn(1, "HeCal_MaxEMyprompt"); // column Id = 101

    analysisManager->CreateNtupleIColumn(1, "HeCal_MaxELate_CopyNo"); // column Id = 102
    analysisManager->CreateNtupleDColumn(1, "HeCal_MaxELate"); // column Id = 103


    analysisManager->FinishNtuple(1);





    // >>> CREATING NTUPLES FOR PARTICLES SIMULATION <<< //

    // --- THIRS Ntuple ID = 2 reconstructed data for particles ---
    analysisManager->CreateNtuple("RHits", "RHits");

    analysisManager->CreateNtupleDColumn(2, "MCEnergy"); // column Id = 0
    analysisManager->CreateNtupleDColumn(2, "MCMomentum"); // column Id = 1
    analysisManager->CreateNtupleDColumn(2, "CaptureT", MD->MyCaptureT); // column Id = 2
    analysisManager->CreateNtupleDColumn(2, "CaptureX", MD->MyCaptureX); // column Id = 3
    analysisManager->CreateNtupleDColumn(2, "CaptureY", MD->MyCaptureY); // column Id = 4
    analysisManager->CreateNtupleDColumn(2, "CaptureZ", MD->MyCaptureZ); // column Id = 5
    analysisManager->CreateNtupleDColumn(2, "CaptureE", MD->MyCaptureE); // column Id = 6
    analysisManager->CreateNtupleDColumn(2, "CaptureF", MD->MyCaptureF); // column Id = 7
    analysisManager->CreateNtupleDColumn(2, "CaptureM", MD->MyCaptureM); // column Id = 8
    analysisManager->CreateNtupleDColumn(2, "MCEvnum_over_GenSurface");  // column Id = 9

    analysisManager->CreateNtupleIColumn(2, "Hod2N");  // column Id = 10
    analysisManager->CreateNtupleIColumn(2, "Hod3N");  // column Id = 11
    analysisManager->CreateNtupleIColumn(2, "HGasN");  // column Id = 12

    analysisManager->CreateNtupleDColumn(2, "MCMomX"); // column Id = 13
    analysisManager->CreateNtupleDColumn(2, "MCMomY"); // column Id = 14
    analysisManager->CreateNtupleDColumn(2, "MCMomZ"); // column Id = 15
    analysisManager->CreateNtupleDColumn(2, "MCPosX"); // column Id = 16
    analysisManager->CreateNtupleDColumn(2, "MCPosY"); // column Id = 17
    analysisManager->CreateNtupleDColumn(2, "MCPosZ"); // column Id = 18

    // vector columns are automatically filled
    analysisManager->CreateNtupleDColumn(2, "Hod2E", fEventAction->GetHod2E()); // column Id = 19
    analysisManager->CreateNtupleDColumn(2, "Hod2T", fEventAction->GetHod2T()); // column Id = 20
    analysisManager->CreateNtupleDColumn(2, "Hod2X", fEventAction->GetHod2X()); // column Id = 21
    analysisManager->CreateNtupleDColumn(2, "Hod2Y", fEventAction->GetHod2Y()); // column Id = 22
    analysisManager->CreateNtupleDColumn(2, "Hod2Z", fEventAction->GetHod2Z()); // column Id = 23

    analysisManager->CreateNtupleDColumn(2, "Hod3E", fEventAction->GetHod3E()); // column Id = 24
    analysisManager->CreateNtupleDColumn(2, "Hod3T", fEventAction->GetHod3T()); // column Id = 25
    analysisManager->CreateNtupleDColumn(2, "Hod3X", fEventAction->GetHod3X()); // column Id = 26
    analysisManager->CreateNtupleDColumn(2, "Hod3Y", fEventAction->GetHod3Y()); // column Id = 27
    analysisManager->CreateNtupleDColumn(2, "Hod3Z", fEventAction->GetHod3Z()); // column Id = 28

    analysisManager->CreateNtupleDColumn(2, "HGasT", fEventAction->GetHGasT()); // column Id = 29
    analysisManager->CreateNtupleDColumn(2, "HGasX", fEventAction->GetHGasX()); // column Id = 30
    analysisManager->CreateNtupleDColumn(2, "HGasY", fEventAction->GetHGasY()); // column Id = 31
    analysisManager->CreateNtupleDColumn(2, "HGasZ", fEventAction->GetHGasZ()); // column Id = 32
    analysisManager->CreateNtupleDColumn(2, "HGasE", fEventAction->GetHGasE()); // column Id = 33
    analysisManager->CreateNtupleIColumn(2, "HGasCopyNo", fEventAction->GetHGasCopyNo()); // column Id = 34

    analysisManager->CreateNtupleDColumn(2, "MCMass"); // column Id = 35
    analysisManager->CreateNtupleDColumn(2, "MCCharge"); // column Id = 36

    analysisManager->CreateNtupleDColumn(2, "HodCaptureT", MD->HodCaptureT); // column Id = 37
    analysisManager->CreateNtupleDColumn(2, "HodCaptureX", MD->HodCaptureX); // column Id = 38
    analysisManager->CreateNtupleDColumn(2, "HodCaptureY", MD->HodCaptureY); // column Id = 39
    analysisManager->CreateNtupleDColumn(2, "HodCaptureZ", MD->HodCaptureZ); // column Id = 40
    analysisManager->CreateNtupleDColumn(2, "HodCaptureE", MD->HodCaptureE); // column Id = 41
    analysisManager->CreateNtupleDColumn(2, "HodCaptureF", MD->HodCaptureF); // column Id = 42
    analysisManager->CreateNtupleDColumn(2, "HodCaptureM", MD->HodCaptureM); // column Id = 43
    analysisManager->CreateNtupleIColumn(2, "HodCaptureVolume", MD->HodCaptureVolume); // column Id = 44
    analysisManager->CreateNtupleIColumn(2, "EventID"); // column Id = 45
    analysisManager->CreateNtupleIColumn(2, "Hod2CopyNo", fEventAction->GetHod2CopyNo()); // column Id = 46
    analysisManager->CreateNtupleIColumn(2, "Hod3CopyNo", fEventAction->GetHod3CopyNo()); // column Id = 47

    //HeCal reconstruction
    analysisManager->CreateNtupleIColumn(2, "HeCal_MaxEprompt_CopyNo"); // column Id = 48
    analysisManager->CreateNtupleDColumn(2, "HeCal_MaxEprompt"); // column Id = 49

    analysisManager->CreateNtupleIColumn(2, "HeCal_MaxELate_CopyNo"); // column Id = 50
    analysisManager->CreateNtupleDColumn(2, "HeCal_MaxELate"); // column Id = 51

    //Hod2 reconstruction
    analysisManager->CreateNtupleIColumn(2, "Hod2_MaxEprompt_CopyNo"); // column Id = 52
    analysisManager->CreateNtupleDColumn(2, "Hod2_MaxEprompt"); // column Id = 53
    analysisManager->CreateNtupleIColumn(2, "Hod2Hit_prompt"); // column Id = 54

    analysisManager->CreateNtupleIColumn(2, "Hod2_MaxELate_CopyNo"); // column Id = 55
    analysisManager->CreateNtupleDColumn(2, "Hod2_MaxELate"); // column Id = 56
    analysisManager->CreateNtupleIColumn(2, "Hod2Hit_Late"); // column Id = 57

    //Hod3 reconstruction
    analysisManager->CreateNtupleIColumn(2, "Hod3_MaxEprompt_CopyNo"); // column Id = 58
    analysisManager->CreateNtupleDColumn(2, "Hod3_MaxEprompt"); // column Id = 59
    analysisManager->CreateNtupleIColumn(2, "Hod3Hit_prompt"); // column Id = 60

    analysisManager->CreateNtupleIColumn(2, "Hod3_MaxELate_CopyNo"); // column Id = 61
    analysisManager->CreateNtupleDColumn(2, "Hod3_MaxELate"); // column Id = 62
    analysisManager->CreateNtupleIColumn(2, "Hod3Hit_Late"); // column Id = 63

    //beta
    analysisManager->CreateNtupleDColumn(2, "beta"); // column Id = 64

    //HeCal capture copy noumber
    analysisManager->CreateNtupleIColumn(2, "CaptureCopyNo", MD->MyCaptureCopyNo); // column Id = 65"

    //Energy released for each Hod3 slab
    analysisManager->CreateNtupleDColumn(2, "Hod3_E_slabPrompt", fEventAction->GetHod3_E_slabPrompt() ); // column Id 66
    analysisManager->CreateNtupleDColumn(2, "Hod3_E_slabLate", fEventAction->GetHod3_E_slabDelayed() ); // column Id 67

    // //Energy released for each Hod2 slab
    analysisManager->CreateNtupleDColumn(2, "Hod2_E_slabPrompt", fEventAction->GetHod2_E_slabPrompt() ); // column Id 68
    analysisManager->CreateNtupleDColumn(2, "Hod2_E_slabLate", fEventAction->GetHod2_E_slabDelayed() ); // column Id 69
    // //Energy released for each HeCal tank
    analysisManager->CreateNtupleDColumn(2, "HeCal_E_tankPrompt", fEventAction->GetHeCal_E_tankPrompt() ); // column Id 70
    analysisManager->CreateNtupleDColumn(2, "HeCal_E_tankLate", fEventAction->GetHeCal_E_tankDelayed() ); // column Id 71

    
    analysisManager->FinishNtuple(2);    


  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDRunAction::~ADHDRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDRunAction::BeginOfRunAction(const G4Run* /*run*/)
{  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file The default file name is set in ADHDRunAction::ADHDRunAction(), it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  // G4cout << "RunAction ended" << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......