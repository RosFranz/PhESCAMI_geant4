// $Id: ADHDDetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file ADHDDetectorConstruction.cc
/// \brief Implementation of the ADHDDetectorConstruction class

#include "ADHDDetectorConstruction.hh"
#include "ADHDHodoscopeSD.hh"
#include "ADHDHeGasSD.hh"
#include "MyDet.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDDetectorConstruction::ADHDDetectorConstruction() : G4VUserDetectorConstruction(),
fHodoscope2Logical(0), fHodoscope3Logical(0), fHeGasLV(0), fVisAttributes(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDDetectorConstruction::~ADHDDetectorConstruction()
{
  for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) delete fVisAttributes[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ADHDDetectorConstruction::Construct()
{

  // Construct materials
  ConstructMaterials();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  G4VSolid* worldSolid = new G4Box("worldBox",10.*m,10.*m,10.*m);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,G4_Galactic,"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0, false,0,checkOverlaps);


  // VESSELS AND PLUGS SPECIFICATIONS 

  int NTankX = 3;
  int NTankYZ = 5;

  MyDet *MD = MyDet::gethead();
  G4double PlugRadius = MD->MyVesselPlugRadius;
  G4double PlugLength = MD->MyVesselPlugLength;
  G4double VesselThick = MD->MyVesselThick;
  G4double tubsRmax = MD->MyVesselRadius;
  G4double tubsRmin = tubsRmax - VesselThick;
  G4double tubsLength = MD->MyVesselTubsLength;
  G4cout << "tubs-length" << tubsLength << G4endl;

  //vessel
  G4VSolid* V_tubs = new G4Tubs("SV_tubs",tubsRmin, tubsRmax, tubsLength, 0., 360*deg);
  G4VSolid* V_Sph1 = new G4Sphere("SV_Sph1", tubsRmin, tubsRmax, 0., 360*deg, 0., 90.*deg);
  G4VSolid* V_Sph2 = new G4Sphere("SV_Sph2", tubsRmin, tubsRmax, 0., 360*deg, 90.*deg, 180.*deg);
  G4VSolid* Ve = new G4UnionSolid("SVe", V_tubs, V_Sph1, 0, G4ThreeVector(0,0,tubsLength));
  G4VSolid* vessel = new G4UnionSolid("SVvessel", Ve, V_Sph2, 0, G4ThreeVector(0,0,-tubsLength));
  G4LogicalVolume* vesselLV = new G4LogicalVolume(vessel, Zylon, "vesselLV");    

  // vessel plug
  G4VSolid* V_plug = new G4Tubs("SV_plug", 0., PlugRadius, PlugLength, 0., 360*deg);
  G4LogicalVolume* plug_LV = new G4LogicalVolume(V_plug, Aluminium, "plug_LV");
  
  // He gas scintillator
  G4VSolid* He_Tubs = new G4Tubs("SHe_tubs", 0., tubsRmin, tubsLength, 0., 360);
  G4VSolid* He_Sph1 = new G4Orb("SHe_Sph1", tubsRmin);
  G4VSolid* He1 = new G4UnionSolid("SHe1", He_Tubs, He_Sph1, 0, G4ThreeVector(0,0,tubsLength));
  G4VSolid* HeGas = new G4UnionSolid("HeGas", He1, He_Sph1, 0, G4ThreeVector(0,0,-tubsLength));
  
  fHeGasLV = new G4LogicalVolume( HeGas,HeGasPress,"HeGasLV");  

  //positioning Vessels, plugs, Helium volumes
  G4double xLength = tubsLength + tubsRmax + PlugLength;

  //old code
  /*G4double CoordinateX = RS2 - (xLength+PlugLength);
  G4double OffsetX = RS2 - NTankX*(xLength+PlugLength);
  G4double CoordinateYZ = RS2 - tubsRmax;
  G4double OffsetYZ = RS2 - NTankYZ*(tubsRmax);*/

  

  G4double CoordinateX = (NTankX-1)*(xLength+PlugLength)+40.; // guardare le due righe sopra per capire come la ottengo
  G4double DeltaX = 2*(xLength+PlugLength)+40.; // lo shift di 40 mm serve per separare maggiormente i tanks

  G4double CoordinateYZ = tubsRmax*(NTankYZ-1);// guardare le due righe sopra per capire come la ottengo
  G4double DeltaYZ = 2*tubsRmax;
  
  G4RotationMatrix* RotationX = new G4RotationMatrix();
  RotationX->rotateX(90.*deg);

  G4RotationMatrix* RotationY = new G4RotationMatrix();
  RotationY->rotateY(90.*deg);

  G4RotationMatrix* negRotationY = new G4RotationMatrix();
  negRotationY->rotateY(0.*deg);

  // vessel, plug and He copy numbers 
  int CopyNumberV = 0;
  int CopyNumberHe = 0;
  int CopyNumberPlug = 0;

  
  if (MD->OnlyHod != 1){ // invisibility
    for(int ix = 0; ix < NTankX; ++ix){
      for(int iy = 0; iy < NTankYZ; ++iy){
        for(int iz = 0; iz < NTankYZ; ++iz){
            // new G4PVPlacement(negRotationY, G4ThreeVector(-CoordinateYZ+(iz*DeltaYZ), -CoordinateYZ+(iy*DeltaYZ), CoordinateX-(ix*DeltaX)+xLength), plug_LV, "plug1", worldLogical, false, CopyNumberPlug++, checkOverlaps);
            // new G4PVPlacement(negRotationY, G4ThreeVector(-CoordinateYZ+(iz*DeltaYZ), -CoordinateYZ+(iy*DeltaYZ), CoordinateX-(ix*DeltaX)), fHeGasLV, "HeGas", worldLogical, false,CopyNumberHe++,checkOverlaps);
            // new G4PVPlacement(negRotationY, G4ThreeVector(-CoordinateYZ+(iz*DeltaYZ), -CoordinateYZ+(iy*DeltaYZ), CoordinateX-(ix*DeltaX)), vesselLV, "vessel", worldLogical, false,CopyNumberV++,checkOverlaps);
            // new G4PVPlacement(negRotationY, G4ThreeVector(-CoordinateYZ+(iz*DeltaYZ), -CoordinateYZ+(iy*DeltaYZ), CoordinateX-(ix*DeltaX)-xLength), plug_LV, "plug2", worldLogical, false, CopyNumberPlug++, checkOverlaps);

            new G4PVPlacement(RotationY, G4ThreeVector(CoordinateX-(ix*DeltaX)+xLength, -CoordinateYZ+(iy*DeltaYZ), -CoordinateYZ+(iz*DeltaYZ)), plug_LV, "plug1", worldLogical, false, CopyNumberPlug++, checkOverlaps);
            new G4PVPlacement(RotationY, G4ThreeVector(CoordinateX-(ix*DeltaX),-CoordinateYZ+(iy*DeltaYZ), -CoordinateYZ+(iz*DeltaYZ)), fHeGasLV, "HeGas", worldLogical, false,CopyNumberHe++,checkOverlaps);
            new G4PVPlacement(RotationY, G4ThreeVector(CoordinateX-(ix*DeltaX),-CoordinateYZ+(iy*DeltaYZ), -CoordinateYZ+(iz*DeltaYZ)), vesselLV, "vessel", worldLogical, false,CopyNumberV++,checkOverlaps);
            new G4PVPlacement(RotationY, G4ThreeVector(CoordinateX-(ix*DeltaX)-xLength, -CoordinateYZ+(iy*DeltaYZ),-CoordinateYZ+(iz*DeltaYZ)), plug_LV, "plug2", worldLogical, false, CopyNumberPlug++, checkOverlaps);
        }
      }
    }
  }



  // HODOSCOPE INNER LAYER

  G4double RS2 = MD->MyRS2; // hodoscope inner layer half length
  G4double Sthick = MD->MySthick;  
  int NBars2 = MD->MyNBars2;
  G4double Slat2 = RS2/float(NBars2);
  G4VSolid* hodoscope2Solid = new G4Box("hodoscope2Box",Slat2,RS2,Sthick); 
  fHodoscope2Logical = new G4LogicalVolume(hodoscope2Solid,scintillator,"hodoscope2Logical");

  G4double maxStep = 10*um; // Maximum step allowed inside hodoscopes
  G4UserLimits *fStepLimit = new G4UserLimits(maxStep);
  fHodoscope2Logical->SetUserLimits(fStepLimit);
  
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++){
      G4double x1 = -1.*RS2 + Slat2 + i*2.*Slat2;
      new G4PVPlacement(0, G4ThreeVector(x1,0.,-1.*(RS2+3*Sthick)), fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i,checkOverlaps);
    }
  // }

  // second side
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++){
      G4double x1 = -1.*RS2 + Slat2 + 2*Sthick + i*2.*Slat2;
      new G4PVPlacement(0, G4ThreeVector(x1,0.,(RS2-Sthick)), fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i + NBars2,checkOverlaps);
    }
  // }

  // third side layer
  if (MD->OnlyHod != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++){
      G4double x1 = -1.*RS2 + Slat2 + i*2.*Slat2;
      new G4PVPlacement(RotationX, G4ThreeVector(x1,-1.*(RS2+Sthick),0.), fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i + 2*NBars2,checkOverlaps);
    }
  }

  // fourth side layer
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++){
      G4double x1 = -1.*RS2 + Slat2 + i*2.*Slat2;
      new G4PVPlacement(RotationX, G4ThreeVector(x1,RS2+Sthick,0.), fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i + 3*NBars2,checkOverlaps);
    }
  // }

  // 5 side layer
  if (MD->OnlyHod != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++){
      G4double z1 = -1.*RS2 + Slat2 + i*2.*Slat2;
      new G4PVPlacement(RotationY, G4ThreeVector(-1.*(RS2-Sthick),0.,z1), fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i + 4*NBars2,checkOverlaps);
    }
  }

  // 6 side layer
  if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars2;i++)  {
      G4double z1 = -1.*(RS2 + 2*Sthick) + Slat2 + i*2.*Slat2;
      new G4PVPlacement(RotationY, G4ThreeVector(RS2-Sthick,0.,z1),fHodoscope2Logical,"hodoscope2Physical",worldLogical,false,i + 5*NBars2,checkOverlaps);
    }
  }

  // HODOSCIOPE EXTERNAL LAYER 

  G4double RS3 = MD->MyRS3;
  int NBars3 = MD->MyNBars3;
  G4double Slat3 = RS3/float(NBars3);
  G4VSolid* hodoscope3Solid = new G4Box("hodoscope3Box",Slat3,RS3,Sthick); 
  fHodoscope3Logical = new G4LogicalVolume(hodoscope3Solid,scintillator,"hodoscope3Logical");

  fHodoscope3Logical->SetUserLimits(fStepLimit);
  
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double x1 = -1.*RS3 + Slat3 + i*2.*Slat3;
      new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.*(RS3+3*Sthick)),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i,checkOverlaps);
    }
  // }

  // second side
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double x1 = -1.*RS3 + Slat3 + 2*Sthick + i*2.*Slat3;
      new G4PVPlacement(0,G4ThreeVector(x1,0.,(RS3-Sthick)),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i + NBars3,checkOverlaps);
    }
  // }

  // third side layer
  if (MD->OnlyHod != 1){ // invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double x1 = -1.*RS3 + Slat3 + i*2.*Slat3;
      new G4PVPlacement(RotationX,G4ThreeVector(x1,-1.*(RS3+Sthick),0.),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i + 2*NBars3,checkOverlaps);
    }
  }

  // fourth side layer
  // if (MD->Visible != 1){ // invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double x1 = -1.*RS3 + Slat3 + i*2.*Slat3;
      new G4PVPlacement(RotationX,G4ThreeVector(x1,(RS3+Sthick),0.),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i + 3*NBars3,checkOverlaps);
    }
  // }

  // 5 side layer
  if (MD->OnlyHod != 1){ // invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double z1 = -1.*RS3 + Slat3 + i*2.*Slat3;
      new G4PVPlacement(RotationY,G4ThreeVector(-1.*(RS3-Sthick),0.,z1),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i + 4*NBars3,checkOverlaps);
    }
  }
  

  // 6 side layer
  if (MD->Visible != 1){// invisibility
    for (G4int i=0;i<NBars3;i++){
      G4double z1 = -1.*(RS3 + 2*Sthick) + Slat3 + i*2.*Slat3;
      new G4PVPlacement(RotationY,G4ThreeVector((RS3-Sthick),0.,z1),fHodoscope3Logical,"hodoscope3Physical",worldLogical,false,i + 5*NBars3,checkOverlaps);
    }
  }


  // VISUALIZATION ATTRIBUTES ------------------------------------------------
  
  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.1,0.1, 1.));
  vesselLV->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(1.0,0.,0.));
  fHodoscope2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(1.,1.0,0.));
  fHodoscope3Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  // return the world physical volume 
  return worldPhysical;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void ADHDDetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  G4VSensitiveDetector* HeGas = new ADHDHeGasSD(SDname="/HeGas");
  SDman->AddNewDetector(HeGas);
  fHeGasLV->SetSensitiveDetector(HeGas);

  G4VSensitiveDetector* hodoscope2 = new ADHDHodoscopeSD(SDname="/hodoscope2");
  SDman->AddNewDetector(hodoscope2);
  fHodoscope2Logical->SetSensitiveDetector(hodoscope2);

  G4VSensitiveDetector* hodoscope3 = new ADHDHodoscopeSD(SDname="/hodoscope3");
  SDman->AddNewDetector(hodoscope3);
  fHodoscope3Logical->SetSensitiveDetector(hodoscope3);
}    


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void ADHDDetectorConstruction::ConstructMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();

  // Vacuum "Air with low density"
  //G4Material* air = G4Material::GetMaterial("G4_AIR");
  //G4double density = 1.0e-5*air->GetDensity();
  //nistManager->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

  G4_Galactic = nistManager->FindOrBuildMaterial("G4_Galactic");

  // He gas with a density different from the one defined in NIST
  //nistManager->FindOrBuildMaterial("G4_He");
  // 
  // G4double density = 0.0267*g/cm3;//(MD->MyHeDensity)*g/cm3; 
  // G4double temparature = 253.15*kelvin; //283.15*kelvin; //10 C
  // G4double pressure = 310.0*bar; // MEOP from tank datasheet
  // HeGasPress = nistManager->BuildMaterialWithNewDensity("ADHD_He","G4_He",density,temparature,pressure);

  G4double a = 4.002*g/mole;
  G4double z = 2;
  MyDet *MD = MyDet::gethead();
  G4double density = (MD->MyHeDensity)*g/cm3;
  // do we really need also kStateGas, tempature and pressure?
  HeGasPress = new G4Material("ADHD_He", z, a, density);// kStateGas, 283.15, 310*bar);//, 288.15*kelvin, 310.*bar);

  // He Vessel
  Zylon = nistManager->FindOrBuildMaterial("G4_NYLON-11_RILSAN"); //density 1.425 similar to zylon 1.5
  // He Vessel plug
  Aluminium = nistManager->FindOrBuildMaterial("G4_Al");
  
  // Scintillator (PolyVinylToluene, C_9H_10)
  scintillator = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
