// $Id: ADHDPrimaryGeneratorAction.cc 77781 2013-11-28 07:54:07Z gcosmo $
//
/// \file ADHDPrimaryGeneratorAction.cc
/// \brief Implementation of the ADHDPrimaryGeneratorAction class

#include "ADHDPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "MyDet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDPrimaryGeneratorAction::ADHDPrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction(), fParticleGun(0), fParticleSource(0),
  fPositron(0), fMuon(0),
  fProton(0), fpbar(0),
  fdbar(0), fDeuton(0),
  fRandPrimaryGun(false),
  fRandPrimaryGPS(true)
{
  G4int n_particle = 1;

  fParticleSource = new G4GeneralParticleSource();
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fPositron = particleTable->FindParticle(particleName="e+");
  fMuon = particleTable->FindParticle(particleName="mu-");

  fProton = particleTable->FindParticle(particleName="proton");
  fpbar = particleTable->FindParticle(-2212);

  fDeuton = particleTable->FindParticle(1000010020);
  fdbar = particleTable->FindParticle(-1000010020);

  fHe4 = particleTable->FindParticle(particleName="alpha");
  


    // proton
  fParticleGun->SetParticleDefinition(fProton);
  NNucleons = 1;

  // Deuteron
  // fParticleGun->SetParticleDefinition(fDeuton);
  // NNucleons = 2; 

  //He4
  // fParticleGun->SetParticleDefinition(fHe4);
  // NNucleons = 4;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ADHDPrimaryGeneratorAction::~ADHDPrimaryGeneratorAction()
{
  delete fParticleSource;
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ADHDPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  MyDet *MD = MyDet::gethead();

  if(fRandPrimaryGPS){
    fParticleSource->GeneratePrimaryVertex(event);
    MD->MCEvnum=MD->MCEvnum+1;

    if(event->GetEventID() == 0){
      G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();
      if(primary->GetPDGcode()<0.) MD->antiparticle = 1;
      NNucleons=1;
      //(deutons) 
      if(fabs(primary->GetPDGcode()) == 1000010020) NNucleons = 2.;
      //He4
      if(fabs(primary->GetPDGcode()) == 1000020040) NNucleons = 4.;
      //C12
      if(fabs(primary->GetPDGcode()) == 1000060120) NNucleons = 12.;
    }

    //save the generated spectra
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4double EkinPerN = fParticleSource->GetParticleEnergy()/NNucleons;
	  man->FillH1(0, EkinPerN);
    return;
  }



  // RANDOMIZATIONN USING PARTICLE GUN

  //FLAT IN EKIN
  /*
  pp = fMomentum + 0.5*fSigmaMomentum; //fMomentum and fSigmaMomentum are not defined
  G4double Ekin2 = std::sqrt(pp*pp+mass*mass)-mass;
  pp = fMomentum - 0.5*fSigmaMomentum;
  G4double Ekin1 = std::sqrt(pp*pp+mass*mass)-mass;
  Ekin = 0.5*(Ekin1+Ekin2)+ (G4UniformRand()-0.5)*(Ekin2-Ekin1);
  */
  

  // ENERGY SPECTRA AS 1/Ekin per nucleon between 10-10'000 MeV ("MIP")
  // specific condition for ordinary particle spectra
  if(fRandPrimaryGun){
    G4double EkinMinPerN = 10. * MeV; // per Nucleon
    G4double EkinPerN = EkinMinPerN*pow(10.,3*G4UniformRand());
    G4double Ekin = EkinPerN*NNucleons;
    fParticleGun->SetParticleEnergy(Ekin);

    //save the generated spectra
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    // G4int id = man->GetH1Id("Gen_spectra");
	  man->FillH1(0, EkinPerN);
  }
  
  
  //ENERGY SPECTRA AS 1/Ekin per nucleon between 2'000-10'000 MeV ("MIP" p)
  // if(fRandPrimaryGun){
    // G4double EkinMinPerN = 2000. * MeV; // per Nucleon
    // G4double EkinPerN = EkinMinPerN * pow(10, log10(5.)*G4UniformRand());
    // G4double Ekin = EkinPerN*NNucleons;
    // fParticleGun->SetParticleEnergy(Ekin);
  // }
  
  
  

  if (fRandPrimaryGun){

    G4double phiangle = (G4UniformRand()-0.5)*360.*deg;
    double Lato = MD->MyRS3+5.*(MD->MySthick); // half length

    // randomize on a squared top plane L=(RS3+5.*Sthick) // just 3Sthick over top plane
    G4double cttq = G4UniformRand(); //estraggo d(cos^2) ovvero 2costheta*dcosteta
    G4double ctt = sqrt(cttq);
    G4double stt = sqrt(1.-cttq);
    G4ThreeVector mdir = G4ThreeVector(stt*std::cos(phiangle), stt*std::sin(phiangle),-1.*ctt);
    fParticleGun->SetParticleMomentumDirection(mdir);

    G4double LatoX = (G4UniformRand()-0.5)*Lato*2.*5.; // dimensioni piano = 5 dimensioni scint. esterno
    G4double LatoY = (G4UniformRand()-0.5)*Lato*2.*5.; // dimensioni piano = 5 dimensioni scint. esterno

    MD->GenSurface=pow(Lato*2.*5,2.); //area totale di generazione

    G4ThreeVector ppos = G4ThreeVector(LatoX, LatoY,Lato);
    fParticleGun->SetParticlePosition(ppos);

    fParticleGun->GeneratePrimaryVertex(event);
    MD->MCEvnum=MD->MCEvnum+1;
    return;
  }


  if(!fRandPrimaryGun && !fRandPrimaryGPS){ // NOT RANDOM position and not isotropic
    double height = MD->MyRS3+5.*(MD->MySthick);
    //position set in order to avoid plug n. 24
    // enters in hod3 n. 96 and hod2 n.96
    // meets the tank n. 12 
    fParticleGun->SetParticlePosition(G4ThreeVector(35.,35.,height));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
    fParticleGun->GeneratePrimaryVertex(event);
    MD->MCEvnum=MD->MCEvnum+1;
    return;
  }

  
}