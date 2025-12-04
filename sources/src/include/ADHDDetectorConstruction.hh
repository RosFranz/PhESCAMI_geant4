// $Id: ADHDDetectorConstruction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDDetectorConstruction.hh
/// \brief Definition of the ADHDDetectorConstruction class

#ifndef ADHDDetectorConstruction_h
#define ADHDDetectorConstruction_h
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"

#include <vector>

class ADHDMagneticField;

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

/// Detector construction

class ADHDDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    ADHDDetectorConstruction();
    virtual ~ADHDDetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    //void SetArmAngle(G4double val);
    //G4double GetArmAngle() { return fArmAngle; }
    
    void ConstructMaterials();
    
private:
    //void DefineCommands();
    
    static G4ThreadLocal ADHDMagneticField* fMagneticField;
    static G4ThreadLocal G4FieldManager* fFieldMgr;
    
    G4LogicalVolume* fHodoscope2Logical;
    G4LogicalVolume* fHodoscope3Logical;
    G4LogicalVolume* fHeGasLV;
    
    std::vector<G4VisAttributes*> fVisAttributes;

    //materials:
    G4Material *G4_Galactic, *HeGasPress, *Zylon, *scintillator, *Aluminium;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
