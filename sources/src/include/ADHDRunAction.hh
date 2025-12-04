// $Id: ADHDRunAction.hh 74204 2013-10-01 07:04:43Z ihrivnac $
// 
/// \file ADHDRunAction.hh
/// \brief Definition of the ADHDRunAction class

#ifndef ADHDRunAction_h
#define ADHDRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class ADHDEventAction;

class G4Run;

/// Run action class

class ADHDRunAction : public G4UserRunAction
{
  public:
    ADHDRunAction(ADHDEventAction* eventAction);
    virtual ~ADHDRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

  private:
    ADHDEventAction* fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
