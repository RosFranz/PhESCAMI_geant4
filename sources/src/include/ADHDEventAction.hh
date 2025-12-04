// $Id: ADHDEventAction.hh 94486 2015-11-19 08:33:37Z gcosmo $
//
/// \file ADHDEventAction.hh
/// \brief Definition of the ADHDEventAction class

#ifndef ADHDEventAction_h
#define ADHDEventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>


/// Event action

class ADHDEventAction : public G4UserEventAction
{
public:
  ADHDEventAction();
  virtual ~ADHDEventAction();
  
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  int HodHitCluster(const G4Event*, int fHHC, int layer); 
  int HeGasHitCluster(const G4Event*, int fHHC); 

  std::vector<G4double>& GetHod2T() { return fHod2T; }
  std::vector<G4double>& GetHod2X() { return fHod2X; }
  std::vector<G4double>& GetHod2Y() { return fHod2Y; }
  std::vector<G4double>& GetHod2Z() { return fHod2Z; }
  std::vector<G4double>& GetHod2E() { return fHod2E; }
  std::vector<G4int>& GetHod2CopyNo() { return fHod2CopyNo; }
  
  std::vector<G4double>& GetHod3T() { return fHod3T; }
  std::vector<G4double>& GetHod3X() { return fHod3X; }
  std::vector<G4double>& GetHod3Y() { return fHod3Y; }
  std::vector<G4double>& GetHod3Z() { return fHod3Z; }
  std::vector<G4double>& GetHod3E() { return fHod3E; }
  std::vector<G4int>& GetHod3CopyNo() { return fHod3CopyNo; }

  std::vector<G4double>& GetHGasT() { return fHGasT; }
  std::vector<G4double>& GetHGasX() { return fHGasX; }
  std::vector<G4double>& GetHGasY() { return fHGasY; }
  std::vector<G4double>& GetHGasZ() { return fHGasZ; }
  std::vector<G4double>& GetHGasE() { return fHGasE; }
  std::vector<G4int>& GetHGasCopyNo() { return fHGasCopyNo; }

  std::vector<G4double>& GetHod3_E_slabPrompt() { return fHod3_E_slabPrompt; }
  std::vector<G4double>& GetHod3_E_slabAcap() { return fHod3_E_slabAcap; }
  std::vector<G4double>& GetHod3_E_slabDelayed() { return fHod3_E_slabDelayed; }

  std::vector<G4double>& GetHod3_E_slabMyPrompt() { return fHod3_E_slabMyPrompt; }
  std::vector<G4double>& GetHod3_E_slabLate() { return fHod3_E_slabLate; }

  std::vector<G4double>& GetHod2_E_slabPrompt() { return fHod2_E_slabPrompt; }
  std::vector<G4double>& GetHod2_E_slabAcap() { return fHod2_E_slabAcap; }
  std::vector<G4double>& GetHod2_E_slabDelayed() { return fHod2_E_slabDelayed; }

  std::vector<G4double>& GetHod2_E_slabMyPrompt() { return fHod2_E_slabMyPrompt; }
  std::vector<G4double>& GetHod2_E_slabLate() { return fHod2_E_slabLate; }

  std::vector<G4double>& GetHeCal_E_tankPrompt() { return fHeCal_E_tankPrompt; }
  std::vector<G4double>& GetHeCal_E_tankAcap() { return fHeCal_E_tankAcap; }
  std::vector<G4double>& GetHeCal_E_tankDelayed() { return fHeCal_E_tankDelayed; }

  std::vector<G4double>& GetHeCal_E_tankMyPrompt() { return fHeCal_E_tankMyPrompt; }
  std::vector<G4double>& GetHeCal_E_tankLate() { return fHeCal_E_tankLate; }

  
    
private:
  G4int fHHC2ID;
  G4int fHHC3ID;
  G4int fHGHCID;

  std::vector<G4double> fHod2T;
  std::vector<G4double> fHod2X;
  std::vector<G4double> fHod2Y;
  std::vector<G4double> fHod2Z;
  std::vector<G4double> fHod2E;
  std::vector<G4int> fHod2CopyNo;

  std::vector<G4double> fHod3T;
  std::vector<G4double> fHod3X;
  std::vector<G4double> fHod3Y;
  std::vector<G4double> fHod3Z;
  std::vector<G4double> fHod3E;
  std::vector<G4int> fHod3CopyNo;

  std::vector<G4double> fHGasT;
  std::vector<G4double> fHGasX;
  std::vector<G4double> fHGasY;
  std::vector<G4double> fHGasZ;
  std::vector<G4double> fHGasE;
  std::vector<G4int> fHGasCopyNo;

  std::vector<G4double> fHod3_E_slabPrompt;
  std::vector<G4double> fHod3_E_slabAcap;
  std::vector<G4double> fHod3_E_slabDelayed;

  std::vector<G4double> fHod3_E_slabMyPrompt;
  std::vector<G4double> fHod3_E_slabLate;

  std::vector<G4double> fHod2_E_slabPrompt;
  std::vector<G4double> fHod2_E_slabAcap;
  std::vector<G4double> fHod2_E_slabDelayed;

  std::vector<G4double> fHod2_E_slabMyPrompt;
  std::vector<G4double> fHod2_E_slabLate;

  std::vector<G4double> fHeCal_E_tankPrompt;
  std::vector<G4double> fHeCal_E_tankAcap;
  std::vector<G4double> fHeCal_E_tankDelayed;

  std::vector<G4double> fHeCal_E_tankMyPrompt;
  std::vector<G4double> fHeCal_E_tankLate;

  G4double NofNucleons;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
