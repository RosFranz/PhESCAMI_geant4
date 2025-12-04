// $Id: ADHDHeGasSD.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \filG4TrackStatue ADHDHeGasSD.hh
/// \brief Definition of the ADHDHeGasSD class

#ifndef ADHDHeGasSD_h
#define ADHDHeGasSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ADHDHeGasHit.hh"

class G4Step;
class G4Track;
class G4HCofThisEvent;
class G4TouchableHistory;

/// HeGas sensitive detector

class ADHDHeGasSD : public G4VSensitiveDetector
{
public:
    ADHDHeGasSD(G4String name);
    virtual ~ADHDHeGasSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
private:
    ADHDHeGasHitsCollection* fHitsCollection;
    G4int fHCID;
    G4double PrimKin;

    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
