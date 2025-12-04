// $Id: exampleADHD.cc 70284 2013-05-28 17:26:43Z perl $
//
/// \file exampleADHD.cc
/// \brief Main program of the analysis/ADHD example

#include "ADHDDetectorConstruction.hh"
#include "ADHDActionInitialization.hh"

#include "G4RunManager.hh"

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "MyDet.hh"
MyDet* MyDet::_head=0;
#include <chrono>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  MyDet *MD = MyDet::gethead();

  MD->MyCaptureT.resize(1,0); // Default
  MD->MyCaptureX.resize(1,0); // Default
  MD->MyCaptureY.resize(1,0); // Default
  MD->MyCaptureZ.resize(1,0); // Default
  MD->MyCaptureE.resize(1,0); // Default
  MD->MyCaptureF.resize(1,0); // Default
  MD->MyCaptureM.resize(1,0); // Default
  MD->MyCaptureCopyNo.resize(1,0); // Default


  MD->HodCaptureT.resize(1,0);
  MD->HodCaptureX.resize(1,0);
  MD->HodCaptureY.resize(1,0);
  MD->HodCaptureZ.resize(1,0);
  MD->HodCaptureE.resize(1,0);
  MD->HodCaptureF.resize(1,0);
  MD->HodCaptureM.resize(1,0);
  MD->HodCaptureVolume.resize(1,0);

  //Setting the engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine());

  // Get the current time in milliseconds
  auto current_time_millis = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

  // Convert it to a long
  long seed = static_cast<long>(current_time_millis);

  // Set the seed
  CLHEP::HepRandom::setTheSeed(seed);
  G4Random::setTheSeed(seed);

  // ONLY IN SINGE THREAD MODE
  G4RunManager* runManager = new G4RunManager;

  // Mandatory user initialization classes
  runManager->SetUserInitialization(new ADHDDetectorConstruction);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;

  //in order to force step limit in Hod3 and Hod2
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());

  //standard EM physics
  physicsList->ReplacePhysics(new G4EmStandardPhysics());
  //more suitable for low and middle enrgies process (CPU consuming)
  // physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  // physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization(new ADHDActionInitialization());

  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) ui = new G4UIExecutive(argc, argv);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
    
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( !ui ) {
    // execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    //UImanager->ApplyCommand("/control/execute gui.mac");
    G4cout << "UI session starts ..." << G4endl;
    ui->SessionStart();
    // delete ui;
    
  }

  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
  // delete ui;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

