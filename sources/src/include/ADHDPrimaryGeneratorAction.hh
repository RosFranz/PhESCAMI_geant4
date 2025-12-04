//$Id: ADHDPrimaryGeneratorAction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDPrimaryGeneratorAction.hh
/// \brief Definition of the ADHDPrimaryGeneratorAction class

#ifndef ADHDPrimaryGeneratorAction_h
#define ADHDPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"

class G4GeneralParticleSource;
class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

/// Primary generator
///
/// A single particle is generated.
/// User can select 
/// - the initial momentum and angle
/// - the momentum and angle spreads
/// - random selection of a particle type from proton, kaon+, pi+, muon+, e+ 


class ADHDPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    ADHDPrimaryGeneratorAction();
    virtual ~ADHDPrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event*);
    
    void SetRandomizeGun(G4bool val) { fRandPrimaryGun = val; }
    void SetRandomizeGPS(G4bool val) { fRandPrimaryGPS = val; }
    G4bool GetRandomizeGun() const { return fRandPrimaryGun; }
    G4bool GetRandomizeGPS() const { return fRandPrimaryGPS; }
    
private:
    void DefineCommands();

    G4GeneralParticleSource *fParticleSource;

    G4ParticleGun* fParticleGun;
    G4ParticleDefinition* fPositron;
    G4ParticleDefinition* fMuon;
    G4ParticleDefinition* fProton;
    G4ParticleDefinition* fpbar;
    G4ParticleDefinition* fdbar;
    G4ParticleDefinition* fDeuton;
    G4ParticleDefinition* fHe4;
    G4double NNucleons;
    G4bool fRandPrimaryGun, fRandPrimaryGPS;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
