// $Id: ADHDActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file ADHDActionInitialization.cc
/// \brief Implementation of the ADHDActionInitialization class

#include "ADHDActionInitialization.hh"
#include "ADHDPrimaryGeneratorAction.hh"
#include "ADHDRunAction.hh"
#include "ADHDEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDActionInitialization::ADHDActionInitialization() : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDActionInitialization::~ADHDActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDActionInitialization::Build() const
{
  SetUserAction(new ADHDPrimaryGeneratorAction);

  ADHDEventAction* eventAction = new ADHDEventAction;
  SetUserAction(eventAction); // EventAction

  SetUserAction(new ADHDRunAction(eventAction)); // RunAction
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
