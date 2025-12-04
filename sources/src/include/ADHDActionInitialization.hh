// $Id: ADHDActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file ADHDActionInitialization.hh
/// \brief Definition of the ADHDActionInitialization class

#ifndef ADHDActionInitialization_h
#define ADHDActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class ADHDActionInitialization : public G4VUserActionInitialization
{
  public:
    ADHDActionInitialization();
    virtual ~ADHDActionInitialization();

    virtual void Build() const;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
