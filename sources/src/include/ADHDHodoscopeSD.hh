// $Id: ADHDHodoscopeSD.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDHodoscopeSD.hh
/// \brief Definition of the ADHDHodoscopeSD class

#ifndef ADHDHodoscopeSD_h
#define ADHDHodoscopeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ADHDHodoscopeHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// Hodoscope sensitive detector

class ADHDHodoscopeSD : public G4VSensitiveDetector
{
public:
    ADHDHodoscopeSD(G4String name);
    virtual ~ADHDHodoscopeSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
private:
    ADHDHodoscopeHitsCollection* fHitsCollection;
    G4int fHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
