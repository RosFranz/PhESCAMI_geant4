// $Id: ADHDHodoscopeHit.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDHodoscopeHit.cc
/// \brief Implementation of the ADHDHodoscopeHit class

#include "ADHDHodoscopeHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<ADHDHodoscopeHit>* ADHDHodoscopeHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDHodoscopeHit::ADHDHodoscopeHit(G4int i): G4VHit(), fId(i), fTime(0), fPos(0), fPLogV(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDHodoscopeHit::~ADHDHodoscopeHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDHodoscopeHit::ADHDHodoscopeHit(const ADHDHodoscopeHit &right): G4VHit()
{
  fId = right.fId;
  fTime = right.fTime;
  fEne = right.fEne;
  fPos = right.fPos;
  fPLogV = right.fPLogV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const ADHDHodoscopeHit& ADHDHodoscopeHit::operator=(const ADHDHodoscopeHit &right)
{
  fId = right.fId;
  fTime = right.fTime;
  fEne = right.fEne;
  fPos = right.fPos;
  fPLogV = right.fPLogV;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int ADHDHodoscopeHit::operator==(const ADHDHodoscopeHit &/*right*/) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDHodoscopeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(2);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* ADHDHodoscopeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store = G4AttDefStore::GetInstance("ADHDHodoscopeHit",isNew);

  if (isNew) {
    (*store)["HitType"] = G4AttDef("HitType","Hit Type","Physics","","G4String");
      
    (*store)["ID"] = G4AttDef("ID","ID","Physics","","G4int");
      
    (*store)["Time"] = G4AttDef("Time","Time","Physics","G4BestUnit","G4double");

    (*store)["Ene"] = G4AttDef("Ene","Energy","Physics","G4BestUnit","G4double");
      
    (*store)["Pos"] = G4AttDef("Pos","Position","Physics","G4BestUnit","G4ThreeVector");
      
    (*store)["LVol"] = G4AttDef("LVol","Logical Volume","Physics","","G4String");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* ADHDHodoscopeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
  
  values->push_back(G4AttValue("HitType","HodoscopeHit",""));
  values->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
  values->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  values->push_back(G4AttValue("Ene",G4BestUnit(fEne,"Energy"),""));
  values->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
  
  if (fPLogV) values->push_back(G4AttValue("LVol",fPLogV->GetName(),""));
  else values->push_back(G4AttValue("LVol"," ",""));
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDHodoscopeHit::Print()
{
  G4cout << "  Hodoscope[" << fId << "] " << fTime/ns << " (nsec) "  << fEne/MeV << " (MeV)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
