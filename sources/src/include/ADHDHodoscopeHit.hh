// $Id: ADHDHodoscopeHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file ADHDHodoscopeHit.hh
/// \brief Definition of the ADHDHodoscopeHit class

#ifndef ADHDHodoscopeHit_h
#define ADHDHodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hodoscope hit
///
/// It records:
/// - the strip ID
/// - the particle time
/// - the strip logical volume, its position
/// - the dep energy

class ADHDHodoscopeHit : public G4VHit
{
public:
    ADHDHodoscopeHit(G4int i);
    ADHDHodoscopeHit(const ADHDHodoscopeHit &right);
    virtual ~ADHDHodoscopeHit();

    const ADHDHodoscopeHit& operator=(const ADHDHodoscopeHit &right);
    int operator==(const ADHDHodoscopeHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    
    G4int GetID() const { return fId; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetEne(G4double val) { fEne = val; }
    G4double GetEne() const { return fEne; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }

    void SetCopyNo(G4int val) { fCopyNo = val; }
    G4int GetCopyNo() const { return fCopyNo; }
    
private:
    G4int fId;
    G4int fCopyNo;
    G4double fTime;
    G4double fEne;
    G4ThreeVector fPos;
    const G4LogicalVolume* fPLogV;
};

typedef G4THitsCollection<ADHDHodoscopeHit> ADHDHodoscopeHitsCollection;

extern G4ThreadLocal G4Allocator<ADHDHodoscopeHit>* ADHDHodoscopeHitAllocator;

inline void* ADHDHodoscopeHit::operator new(size_t)
{
    if (!ADHDHodoscopeHitAllocator) ADHDHodoscopeHitAllocator = new G4Allocator<ADHDHodoscopeHit>;
    return (void*)ADHDHodoscopeHitAllocator->MallocSingle();
}

inline void ADHDHodoscopeHit::operator delete(void*aHit)
{
    ADHDHodoscopeHitAllocator->FreeSingle((ADHDHodoscopeHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
