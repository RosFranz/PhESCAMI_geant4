#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <map>
#include <TH2.h>
#include <TLine.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <iostream>
using namespace std;

/*
Programma con percorso al file da "ricostuire" come parametro di ingresso.
Legge il tree, seleziona solo gli eventi con catture in He e aggiunge allo
stesso tree delle branch in piu' come il max E release in HeCal per ogni evento.
i.e. max E in HeCal, max E (s1), max E(s2), #hit (S1+S2).
La trattazione diventa specifica per particelle prive di catture!
DA USARE per file CON DIVERSE SPECIE DI PARTICELLE.
*/

void ADHD_rec(TString fname){

  double MyDelay = 900000.;
  //vessel (He) dimensions
  double MyVesselThick = 11.; //  vessel thickness
  double MyVesselRadius = 216.;//432./2.; // 216mm vessel radius
  double MyVesselLength = 660.; // 660mm vessel length
  double MyVesselTubsLength = (MyVesselLength - (MyVesselRadius*2.))/2.; // 148mm vessel tubs half length
  double MyVesselPlugLength = 31.; // 62/2 mm vessel plug half legth
  double MyVesselPlugRadius = 31.; // 31 mm vessel plug radius

  double MyDistance = 200.; // 200 mm distance between the scintillator layers
  double MyRS3 = 1510.; // mm External ScG4int Half Lenght   
  double MyRS2 = MyRS3 - MyDistance; // mm internal ScG4int Half Lenght 

  const double time = 600*pow(10,3.); // prompt events < 600 micro-seconds    

  TFile *File = TFile::Open(fname.Data(),"update");
  TTree* ADHD = (TTree*) File->Get("Hits");
  
  double MCEnergy;
  int HGasN, Hod2N, Hod3N, EventID;

  // ---- SETTING THE BRANCHES ---- //
  ADHD->SetBranchAddress("MCEnergy",&MCEnergy);
  ADHD->SetBranchAddress("HGasN",&HGasN);
  ADHD->SetBranchAddress("Hod2N",&Hod2N);
  ADHD->SetBranchAddress("Hod3N",&Hod3N);
  ADHD->SetBranchAddress("EventID",&EventID);

    // Inner scintillator
  vector<double> *Hod2E=0, *Hod2T=0, *Hod2X=0, *Hod2Y=0, *Hod2Z=0;
  vector<int> *Hod2CopyNo=0;
  TBranch *b_Hod2E=0, *b_Hod2T=0, *b_Hod2CopyNo=0, *b_Hod2X=0, *b_Hod2Y=0, *b_Hod2Z=0;
  ADHD->SetBranchAddress("Hod2E",&Hod2E,&b_Hod2E);
  ADHD->SetBranchAddress("Hod2T",&Hod2T,&b_Hod2T);
  ADHD->SetBranchAddress("Hod2X",&Hod2X,&b_Hod2X);
  ADHD->SetBranchAddress("Hod2Y",&Hod2Y,&b_Hod2Y);
  ADHD->SetBranchAddress("Hod2Z",&Hod2Z,&b_Hod2Z);  
  ADHD->SetBranchAddress("Hod2CopyNo",&Hod2CopyNo,&b_Hod2CopyNo); 

    // Outer scintillator
  vector<double> *Hod3E=0, *Hod3T=0, *Hod3X=0, *Hod3Y=0, *Hod3Z=0;
  vector<int> *Hod3CopyNo=0;
  TBranch *b_Hod3E=0, *b_Hod3T=0, *b_Hod3CopyNo=0, *b_Hod3X=0, *b_Hod3Y=0, *b_Hod3Z=0;
  ADHD->SetBranchAddress("Hod3E",&Hod3E,&b_Hod3E);
  ADHD->SetBranchAddress("Hod3T",&Hod3T,&b_Hod3T);
  ADHD->SetBranchAddress("Hod3X",&Hod3X,&b_Hod3X);
  ADHD->SetBranchAddress("Hod3Y",&Hod3Y,&b_Hod3Y);
  ADHD->SetBranchAddress("Hod3Z",&Hod3Z,&b_Hod3Z);
  ADHD->SetBranchAddress("Hod3CopyNo",&Hod3CopyNo,&b_Hod3CopyNo);

      // array of the He scintillator
  vector<double> *HGasE=0, *HGasT=0;
  TBranch *b_HGasE=0, *b_HGasT=0;
  vector<int> *HGasCopyNo=0;
  TBranch *b_HGasCopyNo=0;
  ADHD->SetBranchAddress("HGasE",&HGasE,&b_HGasE);
  ADHD->SetBranchAddress("HGasT",&HGasT,&b_HGasT);
  ADHD->SetBranchAddress("HGasCopyNo",&HGasCopyNo,&b_HGasCopyNo);

    // array of the capture
  vector<double> *CaptureT=0, *CaptureM=0;
  TBranch *b_CaptureT=0, *b_CaptureM=0;
  ADHD->SetBranchAddress("CaptureT",&CaptureT,&b_CaptureT);
  ADHD->SetBranchAddress("CaptureM",&CaptureM,&b_CaptureM); 

  const int NScint = 64*6;
  double hitScint[NScint];

  const int NTanks = 75;
  double hitE[NTanks];


  // ---- CREATING NEW BRANCHES ---- //

    // hod3
  double Hod3_prompt_E=0, Hod3_Acap_E=0, Hod3_Delayed_E=0;
  TString branchName = "Hod3_prompt_E";
  TBranch *Br_Hod3_prompt_E = ADHD->Branch(branchName.Data(),&Hod3_prompt_E,branchName+"/D");
  branchName = "Hod3_Acap_E";
  TBranch *Br_Hod3_Acap_E = ADHD->Branch(branchName.Data(),&Hod3_Acap_E,branchName+"/D");
  branchName = "Hod3_Delayed_E";
  TBranch *Br_Hod3_Delayed_E = ADHD->Branch(branchName.Data(),&Hod3_Delayed_E,branchName+"/D");

    // hod2
  double Hod2_prompt_E=0, Hod2_Acap_E=0, Hod2_Delayed_E=0;
  branchName = "Hod2_prompt_E";
  TBranch *Br_Hod2_prompt_E = ADHD->Branch(branchName.Data(),&Hod2_prompt_E,branchName+"/D");
  branchName = "Hod2_Acap_E";
  TBranch *Br_Hod2_Acap_E = ADHD->Branch(branchName.Data(),&Hod2_Acap_E,branchName+"/D");
  branchName = "Hod2_Delayed_E";
  TBranch *Br_Hod2_Delayed_E = ADHD->Branch(branchName.Data(),&Hod2_Delayed_E,branchName+"/D");

    //hit number
  int HodHit_prompt=0, HodHit_Acap_new=0, HodHit_Delayed_new=0;
  branchName = "HodHit_prompt";
  TBranch *Br_HodHit_prompt = ADHD->Branch(branchName.Data(),&HodHit_prompt,branchName+"/I");
  branchName = "HodHit_Acap_new";
  TBranch *Br_HodHit_Acap_new = ADHD->Branch(branchName.Data(),&HodHit_Acap_new,branchName+"/I");
  branchName = "HodHit_Delayed_new";
  TBranch *Br_HodHit_Delayed_new = ADHD->Branch(branchName.Data(),&HodHit_Delayed_new,branchName+"/I");

    // beta measurements
  double beta=0;
  branchName = "beta";
  TBranch *Br_beta = ADHD->Branch(branchName.Data(),&beta,branchName+"/D");
  double dist;


    // HeCal
  double HeCal_prompt_E=0, HeCal_Acap_E=0, HeCal_Delayed_E=0;
  branchName= "HeCal_prompt_E";
  TBranch *Br_HeCal_prompt_E = ADHD->Branch(branchName.Data(),&HeCal_prompt_E,branchName+"/D");
  branchName = "HeCal_Acap_E";
  TBranch *Br_HeCal_Acap_E = ADHD->Branch(branchName.Data(),&HeCal_Acap_E,branchName+"/D");
  branchName = "HeCal_Delayed_E";
  TBranch *Br_HeCal_Delayed_E = ADHD->Branch(branchName.Data(),&HeCal_Delayed_E,branchName+"/D");


    //Nozzoli
  Float_t Cosi[3];
  branchName=Form("Cosi[%d]/F",3);
  TBranch *BrCosi = ADHD->Branch("BrCosi",Cosi,branchName.Data());


  // ---- CYCLE ON THE TREE ENTRIES ---- //

  Long64_t nev = ADHD->GetEntries();
  for (Int_t i = 0; i < nev; i++) {

    // initialization
    HeCal_prompt_E = 0, HeCal_Acap_E = 0, HeCal_Delayed_E=0;
    Hod3_prompt_E = 0, Hod3_Acap_E = 0, Hod3_Delayed_E = 0;
    Hod2_prompt_E = 0, Hod2_Acap_E = 0, Hod2_Delayed_E = 0;
    HodHit_prompt=0, HodHit_Acap_new=0, HodHit_Delayed_new=0;
    beta=0;
    for(int tank=0; tank<NTanks; tank++) hitE[tank] = 0;
    for(int hod=0; hod<NScint; hod++) hitScint[hod] = 0;
    
    // events selection and hit sum
    ADHD->GetEntry(i);

    // ---- cycle on the FIRST capture ---- //
    for(UInt_t k=0; k<CaptureT->size(); k++){
      double time = CaptureT->at(k);

      // ---- HECAL ---- //
        //prompt energy
      for(UInt_t j=0; j<HGasN; j++){
          int copy = HGasCopyNo->at(j);
          if(HGasE->at(j)>0. && HGasT->at(j)<time) hitE[copy] = hitE[copy] + HGasE->at(j);
      }
      for(int tank=0; tank<NTanks; tank++) {
        if(hitE[tank]>HeCal_prompt_E) HeCal_prompt_E = hitE[tank];
        hitE[tank] = 0;
      }
      Br_HeCal_prompt_E->Fill();

        //after capture and before delayed events
      for(UInt_t j=0; j<HGasN; j++){
        int copy = HGasCopyNo->at(j);
        if(HGasE->at(j)>0. && HGasT->at(j)>time && HGasT->at(j) < (time+MyDelay)) hitE[copy] = hitE[copy] + HGasE->at(j);
      }
      for(int tank=0; tank<NTanks; tank++) {
        if(hitE[tank]>HeCal_Acap_E) HeCal_Acap_E = hitE[tank];
        hitE[tank] = 0;
      }
      Br_HeCal_Acap_E->Fill();

        // capture delayed events
      for(UInt_t j=0; j<HGasN; j++){
        int copy = HGasCopyNo->at(j);
        if(HGasE->at(j)>0. && HGasT->at(j) > (time+MyDelay)) hitE[copy] = hitE[copy] + HGasE->at(j);
      }
      for(int tank=0; tank<NTanks; tank++) {
        if(hitE[tank]>HeCal_Delayed_E) HeCal_Delayed_E = hitE[tank];
        hitE[tank] = 0;
      }
      Br_HeCal_Delayed_E->Fill();


      // ---- HODOSCOPES ---- //

        //prompt energy and hit Hod3
      for(UInt_t j=0; j<Hod3N; j++){
        int copy3 = Hod3CopyNo->at(j);
        if(Hod3T->at(j)<time && Hod3E->at(j)!=0){
          hitScint[copy3] = hitScint[copy3] + Hod3E->at(j);
          HodHit_prompt++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(hitScint[hod]>0 && hitScint[hod]>Hod3_prompt_E) Hod3_prompt_E = hitScint[hod];
        hitScint[hod] = 0;
      }
        //prompt energy and hit Hod2
      for(UInt_t j=0; j<Hod2N; j++){
        int copy2 = Hod2CopyNo->at(j);
        if(Hod2T->at(j)<time && Hod2E->at(j)!=0) {
          hitScint[copy2] = hitScint[copy2] + Hod2E->at(j);
          HodHit_prompt++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
          if(hitScint[hod]>0 && hitScint[hod]>Hod2_prompt_E) Hod2_prompt_E = hitScint[hod];
          hitScint[hod] = 0;
      }

      Br_Hod2_prompt_E->Fill();
      Br_Hod3_prompt_E->Fill();
      Br_HodHit_prompt->Fill();


        // after capture and before delayed events Hod3
      for(UInt_t j=0; j<Hod3N; j++){
        int copy3 = Hod3CopyNo->at(j);
        if(Hod3T->at(j)>time && Hod3T->at(j)<(time+MyDelay) && Hod3E->at(j)!=0){
          hitScint[copy3] = hitScint[copy3] + Hod3E->at(j);
          HodHit_Acap_new++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(hitScint[hod]>0 && hitScint[hod]>Hod3_Acap_E) Hod3_Acap_E = hitScint[hod];
        hitScint[hod] = 0;
      }
        // after capture and before delayed events Hod2
      for(UInt_t j=0; j<Hod2N; j++){
        int copy2 = Hod2CopyNo->at(j);
        if(Hod2T->at(j)>time && Hod2T->at(j)<(time+MyDelay) && Hod2E->at(j)!=0) {
          hitScint[copy2] = hitScint[copy2] + Hod2E->at(j);
          HodHit_Acap_new++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
          if(hitScint[hod]>0 && hitScint[hod]>Hod2_Acap_E) Hod2_Acap_E = hitScint[hod];
          hitScint[hod] = 0;
      }

      Br_Hod2_Acap_E->Fill();
      Br_Hod3_Acap_E->Fill();
      Br_HodHit_Acap_new->Fill();

      

        // capture delayed events Hod3
      for(UInt_t j=0; j<Hod3N; j++){
        int copy3 = Hod3CopyNo->at(j);
        if(Hod3T->at(j)>(time+MyDelay) && Hod3E->at(j)!=0){
          hitScint[copy3] = hitScint[copy3] + Hod3E->at(j);
          HodHit_Delayed_new++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
        if(hitScint[hod]>0 && hitScint[hod]>Hod3_Delayed_E) Hod3_Delayed_E = hitScint[hod];
        hitScint[hod] = 0;
      }
        // capture delayed events Hod2
      for(UInt_t j=0; j<Hod2N; j++){
        int copy2 = Hod2CopyNo->at(j);
        if(Hod2T->at(j)>(time+MyDelay) && Hod2E->at(j)!=0) {
          hitScint[copy2] = hitScint[copy2] + Hod2E->at(j);
          HodHit_Delayed_new++;
        }
      }
      for(int hod=0; hod<NScint; hod++){
          if(hitScint[hod]>0 && hitScint[hod]>Hod2_Delayed_E) Hod2_Delayed_E = hitScint[hod];
          hitScint[hod] = 0;
      }

      Br_Hod2_Delayed_E->Fill();
      Br_Hod3_Delayed_E->Fill();
      Br_HodHit_Delayed_new->Fill();
      
      
      break; // ONLY THE FIRST CAPTURE

    } // end of the cycle on the first captures


    // ---- BETA EVALUATION ---- //
    if(Hod3N>0 && Hod2N>0){
      dist = sqrt(pow((Hod2Y->at(0)-Hod3Y->at(0)),2) + pow((Hod2X->at(0)-Hod3X->at(0)),2) + pow((Hod2Z->at(0)-Hod3Z->at(0)),2));
      beta = (1./299.792458)*dist/(Hod2T->at(0)-Hod3T->at(0));  
    }
    Br_beta->Fill();
  
    // evaluate incidence on ToF shell
    Cosi[0]=1;
    Cosi[1]=1;
    Cosi[2]=1; // cosine WRT outher shell
    double dir[3];
    if(Hod3N>0 && Hod2N>0){
      dir[0] = Hod2X->at(0)-Hod3X->at(0);
      dir[1] = Hod2Y->at(0)-Hod3Y->at(0);
      dir[2] = Hod2Z->at(0)-Hod3Z->at(0); 
    }
    else {dir[0]=1;dir[1]=1;dir[2]=1;}
    double len = sqrt(pow(dir[0],2) + pow(dir[1],2) + pow(dir[2],2)); 
    for(int jj = 0;jj<3;jj++) dir[jj]=dir[jj]/len;
    if(Hod3N>0){
      if(fabs(fabs(Hod3Z->at(0))-MyRS3)<0.5) Cosi[2]=fabs(dir[2]); // select SIDE
      if(fabs(fabs(Hod3Y->at(0))-MyRS3)<0.5) Cosi[2]=fabs(dir[1]);
      if(fabs(fabs(Hod3X->at(0))-MyRS3)<0.5) Cosi[2]=fabs(dir[0]);
    }
    if(Hod2N>0){
      if(fabs(fabs(Hod2Z->at(0))-MyRS2)<0.5) Cosi[1]=fabs(dir[2]); // select SIDE
      if(fabs(fabs(Hod2Y->at(0))-MyRS2)<0.5) Cosi[1]=fabs(dir[1]);
      if(fabs(fabs(Hod2X->at(0))-MyRS2)<0.5) Cosi[1]=fabs(dir[0]);
    }

    BrCosi->Fill();

    printf("\rProgress percentage %d / 100",int((i*1.0)/nev*100.));
    if(i==nev-1) printf("\n");
  }
  ADHD->Write("",TObject::kOverwrite);
  File->Write();
  delete File;
  return;
}
