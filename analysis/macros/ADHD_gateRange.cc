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
#include <TRandom3.h>
using namespace std;


/*
Programma senza argomenti in ingresso. Usa dei percorsi predefiniti ai file root "raw" per antiparticelle.
Recupera il tree dai file, e permette di customizzare la ricostruzione usando un valore variabile del gate
che distingue gli eventi prompt da quelli delayed.
Si aggiunge uno smearing gaussiano sulle energie: plastici = 5% HeCal = 10%
Si aggiunge uno smearing gaussiano in beta del 5%
Il valore di default e' 50 ns.
Le informazioni raccolte vengono salvate nel file di output ../sources/ROOT_files/particella_Neventi_rec_valoreGate.root
*/





//!!!!!!!!!!!!!!
// MAIN PROGRAM
//!!!!!!!!!!!!!!

void ADHD_gateRange()
{
    gStyle->SetOptStat(0); // no statistics

    // ----------------- INPUT FILE ----------------- //
    const int nFiles = 18;
    // const int nFiles = 2;
    // const int nFiles = 7;

    TString fInNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_STAT_gps_1raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_1raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_2raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_2raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_3raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_3raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_4raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_4raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_5raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_5raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_6raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_6raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_7raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_7raw.root",
                                "../../sources/ROOT_files/Dbar_10e6_STAT_gps_8raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_8raw.root", 
                                "../../sources/ROOT_files/Dbar_10e7_gps_raw.root", "../../sources/ROOT_files/pbar_10e7_gps_raw.root"};

    // TString fInNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_STAT_gps_1raw.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_1raw.root"};
    
    // TString fInNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_raw.root", "../../sources/ROOT_files/pbar_10e7_gps_raw.root",
    //                           "../../sources/ROOT_files/proton_10e7_gps_rec.root", "../../sources/ROOT_files/Deuteron_10e7_gps_rec.root",
    //                           "../../sources/ROOT_files/He4_10e7_gps_rec.root", "../../sources/ROOT_files/C12_10e7_gps_rec.root",
    //                           "../../sources/ROOT_files/electron_10e7_gps_rec.root"};



    // Gate intervals
    const double Cgate = 50; // ns
    double gate = 0;
    TString Sgate = to_string(int(Cgate));


    // ----------------- OUTPUT FILE ----------------- //
    // 

    TString fOutNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_gps_SM1rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM1rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM2rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM2rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM3rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM3rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM4rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM4rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM5rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM5rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM6rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM6rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM7rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM7rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e6_gps_SM8rec_TESTdir_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM8rec_TESTdir_.root",
                                 "../../sources/ROOT_files/Dbar_10e7_gps_rec_.root", "../../sources/ROOT_files/pbar_10e7_gps_rec_."};

    // TString fOutNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_gps_SM1rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM1rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM2rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM2rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM3rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM3rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM4rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM4rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM5rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM5rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM6rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM6rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM7rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM7rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM8rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM8rec_.root",
    //                              "../../sources/ROOT_files/Deuteron_10e7_gps_SMrec_.root", "../../sources/ROOT_files/proton_10e7_gps_SMrec_.root"};
    

    // TString fOutNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_rec_.root", "../../sources/ROOT_files/pbar_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/proton_10e7_gps_rec_.root", "../../sources/ROOT_files/Deuteron_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/He4_10e7_gps_rec_.root", "../../sources/ROOT_files/C12_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/electron_10e7_gps_rec_.root"};

    for(int i=0; i<nFiles; i++) fOutNames[i].ReplaceAll(".root", "_"+Sgate+"ns.root");

    TString branchName;


    // Constants
    const double MyDelay = 900000; // ns
    //number of slabs = 64*6
    const int NScint = 64*6;
    //number of tanks = 75
    const int NTanks = 75;

    //Random generator
    TRandom3 *r3 = new TRandom3(0);



    // CYCLE ON THE FILES
    for (int i = 0; i < nFiles; i++){

        TFile *f = new TFile(fInNames[i], "OPEN");
        if (!f){
            printf("File %s not found\n", fInNames[i].Data());
            return;
        }

         // --- getting the tree and setting specific branches --- //
        TTree *tr;
        if(fInNames[i].Contains("bar")) tr = (TTree *)f->Get("Hits_bar");
        //particle
        else tr = (TTree *)f->Get("RHits");

        std::printf("\nCopying tree from file %s \n", fInNames[i].Data());

        gROOT->cd();
        TFile* outputFile = TFile::Open(fOutNames[i], "RECREATE");
        outputFile->cd();
        TTree *t = (TTree*) tr->CloneTree();
        
        // ----------------- setting the branches from the tree ----------------- //

        double MCEnergy;
        int HGasN, Hod2N, Hod3N, EventID;

        t->SetBranchAddress("MCEnergy",&MCEnergy);
        t->SetBranchAddress("HGasN",&HGasN);
        t->SetBranchAddress("Hod2N",&Hod2N);
        t->SetBranchAddress("Hod3N",&Hod3N);
        t->SetBranchAddress("EventID",&EventID);

          // Inner scintillator
        vector<double> *Hod2E=0, *Hod2T=0, *Hod2X=0, *Hod2Y=0, *Hod2Z=0;
        vector<int> *Hod2CopyNo=0;
        TBranch *b_Hod2E=0, *b_Hod2T=0, *b_Hod2CopyNo=0, *b_Hod2X=0, *b_Hod2Y=0, *b_Hod2Z=0;
        t->SetBranchAddress("Hod2E",&Hod2E,&b_Hod2E);
        t->SetBranchAddress("Hod2T",&Hod2T,&b_Hod2T);
        t->SetBranchAddress("Hod2X",&Hod2X,&b_Hod2X);
        t->SetBranchAddress("Hod2Y",&Hod2Y,&b_Hod2Y);
        t->SetBranchAddress("Hod2Z",&Hod2Z,&b_Hod2Z);  
        t->SetBranchAddress("Hod2CopyNo",&Hod2CopyNo,&b_Hod2CopyNo); 

        // Outer scintillator
        vector<double> *Hod3E=0, *Hod3T=0, *Hod3X=0, *Hod3Y=0, *Hod3Z=0;
        vector<int> *Hod3CopyNo=0;
        TBranch *b_Hod3E=0, *b_Hod3T=0, *b_Hod3CopyNo=0, *b_Hod3X=0, *b_Hod3Y=0, *b_Hod3Z=0;
        t->SetBranchAddress("Hod3E",&Hod3E,&b_Hod3E);
        t->SetBranchAddress("Hod3T",&Hod3T,&b_Hod3T);
        t->SetBranchAddress("Hod3X",&Hod3X,&b_Hod3X);
        t->SetBranchAddress("Hod3Y",&Hod3Y,&b_Hod3Y);
        t->SetBranchAddress("Hod3Z",&Hod3Z,&b_Hod3Z);
        t->SetBranchAddress("Hod3CopyNo",&Hod3CopyNo,&b_Hod3CopyNo);

        // array of the He scintillator
        vector<double> *HGasE=0, *HGasT=0;
        TBranch *b_HGasE=0, *b_HGasT=0;
        vector<int> *HGasCopyNo=0;
        TBranch *b_HGasCopyNo=0;
        t->SetBranchAddress("HGasE",&HGasE,&b_HGasE);
        t->SetBranchAddress("HGasT",&HGasT,&b_HGasT);
        t->SetBranchAddress("HGasCopyNo",&HGasCopyNo,&b_HGasCopyNo);

            // array of the capture
        vector<double> *CaptureT=0, *CaptureM=0;
        TBranch *b_CaptureT=0, *b_CaptureM=0;
        t->SetBranchAddress("CaptureT",&CaptureT,&b_CaptureT);
        t->SetBranchAddress("CaptureM",&CaptureM,&b_CaptureM); 

        const int NScint = 64*6;
        double hitScint[NScint];

        const int NTanks = 75;
        double hitE[NTanks];



        //variable for the new branches
        //ANTIPARTICLES

        if(fInNames[i].Contains("bar")){
                        
            //HOD3 reconstruction WITHOUT captures
                //Energy released for each Hod3 slab
                vector<double> *Hod3_E_slabMyPrompt=0, *Hod3_E_slabLate=0;

                branchName = "Hod3_E_slabMyPrompt";
                TBranch *Br_Hod3_E_slabMyPrompt = t->Branch(branchName.Data(),&Hod3_E_slabMyPrompt);
                branchName = "Hod3_E_slabLate";
                TBranch *Br_Hod3_E_slabLate = t->Branch(branchName.Data(),&Hod3_E_slabLate);


                int Hod3_MaxEMyprompt_CopyNo=0, Hod3Hit_Myprompt=0,
                    Hod3_MaxELate_CopyNo=0, Hod3Hit_Late=0;
                double Hod3_MaxEMyprompt=0, Hod3_MaxELate=0;
                double TOTHod3Prompt=0;

                branchName = "Hod3_MaxEMyprompt_CopyNo";
                TBranch *Br_Hod3_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&Hod3_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "Hod3_MaxEMyprompt";
                TBranch *Br_Hod3_MaxEMyprompt = t->Branch(branchName.Data(),&Hod3_MaxEMyprompt,branchName+"/D");
                branchName = "Hod3Hit_Myprompt";
                TBranch *Br_Hod3Hit_Myprompt = t->Branch(branchName.Data(),&Hod3Hit_Myprompt,branchName+"/I");
                branchName = "TOTHod3Prompt";
                TBranch *Br_TOTHod3Prompt = t->Branch(branchName.Data(),&TOTHod3Prompt,branchName+"/D");
                
                branchName = "Hod3_MaxELate_CopyNo";
                TBranch *Br_Hod3_MaxELate_CopyNo = t->Branch(branchName.Data(),&Hod3_MaxELate_CopyNo,branchName+"/I");
                branchName = "Hod3_MaxELate";
                TBranch *Br_Hod3_MaxELate = t->Branch(branchName.Data(),&Hod3_MaxELate,branchName+"/D");
                branchName = "Hod3Hit_Late";
                TBranch *Br_Hod3Hit_Late = t->Branch(branchName.Data(),&Hod3Hit_Late,branchName+"/I");
            
            //HOD2 reconstruction WITHOUT captures
                //Energy released for each Hod2 slab WITHOUT captures
                vector<double> *Hod2_E_slabMyPrompt=0, *Hod2_E_slabLate=0;

                branchName = "Hod2_E_slabMyPrompt";
                TBranch *Br_Hod2_E_slabMyPrompt = t->Branch(branchName.Data(),&Hod2_E_slabMyPrompt);
                branchName = "Hod2_E_slabLate";
                TBranch *Br_Hod2_E_slabLate = t->Branch(branchName.Data(),&Hod2_E_slabLate);
                

                int Hod2_MaxEMyprompt_CopyNo=0, Hod2Hit_Myprompt=0,
                    Hod2_MaxELate_CopyNo=0, Hod2Hit_Late=0;
                double Hod2_MaxEMyprompt=0, Hod2_MaxELate=0;
                double TOTHod2Prompt=0;

                branchName = "Hod2_MaxEMyprompt_CopyNo";
                TBranch *Br_Hod2_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&Hod2_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "Hod2_MaxEMyprompt";
                TBranch *Br_Hod2_MaxEMyprompt = t->Branch(branchName.Data(),&Hod2_MaxEMyprompt,branchName+"/D");
                branchName = "Hod2Hit_Myprompt";
                TBranch *Br_Hod2Hit_Myprompt = t->Branch(branchName.Data(),&Hod2Hit_Myprompt,branchName+"/I");
                branchName = "TOTHod2Prompt";
                TBranch *Br_TOTHod2Prompt = t->Branch(branchName.Data(),&TOTHod2Prompt,branchName+"/D");
                
                branchName = "Hod2_MaxELate_CopyNo";
                TBranch *Br_Hod2_MaxELate_CopyNo = t->Branch(branchName.Data(),&Hod2_MaxELate_CopyNo,branchName+"/I");
                branchName = "Hod2_MaxELate";
                TBranch *Br_Hod2_MaxELate = t->Branch(branchName.Data(),&Hod2_MaxELate,branchName+"/D");
                branchName = "Hod2Hit_Late";
                TBranch *Br_Hod2Hit_Late = t->Branch(branchName.Data(),&Hod2Hit_Late,branchName+"/I");

            //HECAL reconstruction without captures
                //Energy released for each HeCal tank WITHOUT captures
                vector<double> *HeCal_E_tankMyPrompt=0, *HeCal_E_tankLate=0;

                branchName = "HeCal_E_tankMyPrompt";
                TBranch *Br_HeCal_E_tankMyPrompt = t->Branch(branchName.Data(),&HeCal_E_tankMyPrompt);
                branchName = "HeCal_E_tankLate";
                TBranch *Br_HeCal_E_tankLate = t->Branch(branchName.Data(),&HeCal_E_tankLate);


                int HeCal_MaxEMyprompt_CopyNo=0, HeCal_MaxELate_CopyNo=0;
                double HeCal_MaxEMyprompt=0, HeCal_MaxELate=0;
                double TOTHeCalPrompt=0, TOTHeCalLate=0; 

                branchName = "HeCal_MaxEMyprompt_CopyNo";
                TBranch *Br_HeCal_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&HeCal_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "HeCal_MaxEMyprompt";
                TBranch *Br_HeCal_MaxEMyprompt = t->Branch(branchName.Data(),&HeCal_MaxEMyprompt,branchName+"/D");
                branchName = "HeCal_MaxELate_CopyNo";
                TBranch *Br_HeCal_MaxELate_CopyNo = t->Branch(branchName.Data(),&HeCal_MaxELate_CopyNo,branchName+"/I");
                branchName = "HeCal_MaxELate";
                TBranch *Br_HeCal_MaxELate = t->Branch(branchName.Data(),&HeCal_MaxELate,branchName+"/D");
                branchName = "TOTHeCalPrompt";
                TBranch *Br_TOTHeCalPrompt = t->Branch(branchName.Data(),&TOTHeCalPrompt,branchName+"/D");
                branchName = "TOTHeCalLate";
                TBranch *Br_TOTHeCalLate = t->Branch(branchName.Data(),&TOTHeCalLate,branchName+"/D");

            //beta 
                double beta=0;
                branchName = "beta";
                TBranch *Br_beta = t->Branch(branchName.Data(),&beta,branchName+"/D");



            // ------------------- READING THE TREE ------------------------ //
            printf("Reading the tree ...\n");
            for (Int_t entry = 0; entry < t->GetEntries(); entry++){

                t->GetEntry(entry);
                gate = Cgate + Hod3T->at(0); // ns

                //HOD3 reconstruction WITHOUT captures
                Hod3_MaxEMyprompt_CopyNo=0; Hod3Hit_Myprompt=0;
                Hod3_MaxELate_CopyNo=0; Hod3Hit_Late=0;
                Hod3_MaxEMyprompt=0; Hod3_MaxELate=0;
                TOTHod3Prompt=0;

                Hod3_E_slabMyPrompt->resize(NScint, 0.);
                Hod3_E_slabLate->resize(NScint, 0.);

                //HOD2
                Hod2_MaxEMyprompt_CopyNo=0; Hod2Hit_Myprompt=0;
                Hod2_MaxELate_CopyNo=0; Hod2Hit_Late=0;
                Hod2_MaxEMyprompt=0; Hod2_MaxELate=0;
                TOTHod2Prompt=0;

                Hod2_E_slabMyPrompt->resize(NScint, 0.);
                Hod2_E_slabLate->resize(NScint, 0.);

                for(int i=0; i<NScint; i++){
                    Hod3_E_slabMyPrompt->at(i)=0;
                    Hod3_E_slabLate->at(i)=0;
                    Hod2_E_slabMyPrompt->at(i)=0;
                    Hod2_E_slabLate->at(i)=0;
                }
                
                //HECAL 
                HeCal_MaxEMyprompt_CopyNo=0; HeCal_MaxELate_CopyNo=0;
                HeCal_MaxEMyprompt=0; HeCal_MaxELate=0;
                TOTHeCalPrompt=0; TOTHeCalLate=0;

                HeCal_E_tankMyPrompt->resize(NTanks, 0.);
                HeCal_E_tankLate->resize(NTanks, 0.);

                for(int i=0; i<NTanks; i++){
                    HeCal_E_tankMyPrompt->at(i)=0;
                    HeCal_E_tankLate->at(i)=0;
                }


                // --- CYCLE ON THE CAPTURES --- //
                // entrance in the for cycle is granted by the CaptureT size that is always at least 1
                // In the absence of a capture the value inside CaptureT is 0
                // The reconstruction for anti particles is still made here, even if no capture happens
                
                for(int k=0; k < int(CaptureT->size()) ; k++){
                    
                    double time = CaptureT->at(k);  // = 0 if no capture
                    if(time==0 && CaptureT->size()>1) continue; // if one capture happens, skip the first entry in the vector
                    double NextCaptureTime = 0;
                    if( k+1 < int(CaptureT->size()) && time!=0 ) NextCaptureTime = CaptureT->at(k+1);

                    //HOD3 reconstruction WITHOUT captures
                        Hod3_MaxEMyprompt_CopyNo=0; Hod3Hit_Myprompt=0;
                        Hod3_MaxELate_CopyNo=0; Hod3Hit_Late=0;
                        Hod3_MaxEMyprompt=0; Hod3_MaxELate=0;
                        TOTHod3Prompt=0;

                    // ---- HOD3  prompt energy ---- //

                    //prompt energy and hit Hod3
                    for(int j=0; j<Hod3N; j++){
                        int copy3 = Hod3CopyNo->at(j);
                        if(Hod3T->at(j)<gate && Hod3E->at(j)!=0){
                            Hod3_E_slabMyPrompt->at(copy3) += Hod3E->at(j);
                            Hod3Hit_Myprompt++;
                        }
                    }

                    for(int hod=0; hod<NScint; hod++){
                        //adding a gaussian smearing
                        Hod3_E_slabMyPrompt->at(hod) = r3->Gaus(Hod3_E_slabMyPrompt->at(hod), 0.05*Hod3_E_slabMyPrompt->at(hod));

                        if(Hod3_E_slabMyPrompt->at(hod)>0.) TOTHod3Prompt += Hod3_E_slabMyPrompt->at(hod);
                        if(Hod3_E_slabMyPrompt->at(hod)>0 && Hod3_E_slabMyPrompt->at(hod)>Hod3_MaxEMyprompt){
                            Hod3_MaxEMyprompt = Hod3_E_slabMyPrompt->at(hod);
                            Hod3_MaxEMyprompt_CopyNo = hod;
                        }
                    }

                    Br_Hod3_E_slabMyPrompt->Fill();
                    Br_Hod3_MaxEMyprompt_CopyNo->Fill(); 
                    Br_Hod3_MaxEMyprompt->Fill();
                    Br_Hod3Hit_Myprompt->Fill();
                    Br_TOTHod3Prompt->Fill();

                    // ---- HOD3 capture delayed events ---- //
                    if(NextCaptureTime==0){ // no double capture
                        for(int j=0; j<Hod3N; j++){
                            int copy3 = Hod3CopyNo->at(j);
                            if(Hod3T->at(j)>gate && Hod3E->at(j)!=0){
                                Hod3_E_slabLate->at(copy3) += Hod3E->at(j);
                                Hod3Hit_Late++;
                            }
                        }
                        for(int hod=0; hod<NScint; hod++){
                            //adding a gaussian smearing
                            Hod3_E_slabLate->at(hod) = r3->Gaus(Hod3_E_slabLate->at(hod), 0.05*Hod3_E_slabLate->at(hod));
                            if(Hod3_E_slabLate->at(hod)>0 && Hod3_E_slabLate->at(hod)>Hod3_MaxELate) {
                                Hod3_MaxELate = Hod3_E_slabLate->at(hod);
                                Hod3_MaxELate_CopyNo = hod;
                            }
                        }
                    }

                    else{// double capture
                        for(int j=0; j<Hod3N; j++){
                            int copy3 = Hod3CopyNo->at(j);
                            if(Hod3T->at(j)>gate && Hod3E->at(j)!=0){
                                Hod3_E_slabLate->at(copy3) += Hod3E->at(j);
                                Hod3Hit_Late++;
                            }
                        }
                        for(int hod=0; hod<NScint; hod++){
                            //adding a gaussian smearing
                            Hod3_E_slabLate->at(hod) = r3->Gaus(Hod3_E_slabLate->at(hod), 0.05*Hod3_E_slabLate->at(hod));
                            if(Hod3_E_slabLate->at(hod)>0 && Hod3_E_slabLate->at(hod)>Hod3_MaxELate) {
                                Hod3_MaxELate = Hod3_E_slabLate->at(hod);
                                Hod3_MaxELate_CopyNo = hod;
                            }
                        }
                    }
                    Br_Hod3_E_slabLate->Fill();
                    Br_Hod3_MaxELate_CopyNo->Fill(); // Hod3 CopyNo for the MAX Delayed energy release 
                    Br_Hod3_MaxELate->Fill();
                    Br_Hod3Hit_Late->Fill();



                    //HOD2 reconstruction without captures
                        Hod2_MaxEMyprompt_CopyNo=0; Hod2Hit_Myprompt=0;
                        Hod2_MaxELate_CopyNo=0; Hod2Hit_Late=0;
                        Hod2_MaxEMyprompt=0; Hod2_MaxELate=0;
                        TOTHod2Prompt=0;

                    // ---- HOD2  prompt energy ---- //
                    for(int j=0; j<Hod2N; j++){
                        int copy2 = Hod2CopyNo->at(j);
                        if(Hod2T->at(j)<gate && Hod2E->at(j)!=0) {
                            Hod2_E_slabMyPrompt->at(copy2) += Hod2E->at(j);
                            Hod2Hit_Myprompt++;
                        }
                    }
                    for(int hod=0; hod<NScint; hod++){
                        //adding a gaussian smearing
                        Hod2_E_slabMyPrompt->at(hod) = r3->Gaus(Hod2_E_slabMyPrompt->at(hod), 0.05*Hod2_E_slabMyPrompt->at(hod));
                        if(Hod2_E_slabMyPrompt->at(hod)>0.) TOTHod2Prompt += Hod2_E_slabMyPrompt->at(hod);
                        if(Hod2_E_slabMyPrompt->at(hod)>0 && Hod2_E_slabMyPrompt->at(hod)>Hod2_MaxEMyprompt){
                            Hod2_MaxEMyprompt = Hod2_E_slabMyPrompt->at(hod);
                            Hod2_MaxEMyprompt_CopyNo = hod;
                        }
                    }

                    Br_Hod2_E_slabMyPrompt->Fill();
                    Br_Hod2_MaxEMyprompt_CopyNo->Fill(); 
                    Br_Hod2_MaxEMyprompt->Fill();
                    Br_Hod2Hit_Myprompt->Fill();
                    Br_TOTHod2Prompt->Fill();



                    
                    //HECAL reconstruction without captures
                        HeCal_MaxEMyprompt_CopyNo=0; HeCal_MaxELate_CopyNo=0;
                        HeCal_MaxEMyprompt=0; HeCal_MaxELate=0;
                        TOTHeCalPrompt=0; TOTHeCalLate=0;


                    // ---- HECAL prompt energy ---- //
                    for(int j=0; j<HGasN; j++){
                        int copy = HGasCopyNo->at(j);
                        if(HGasE->at(j)>0. && HGasT->at(j)<gate) HeCal_E_tankMyPrompt->at(copy) += HGasE->at(j);
                    }
                    for(int tank=0; tank<NTanks; tank++) {
                        //adding a gaussian smearing
                        HeCal_E_tankMyPrompt->at(tank) = r3->Gaus(HeCal_E_tankMyPrompt->at(tank), 0.1*HeCal_E_tankMyPrompt->at(tank));
                        if(HeCal_E_tankMyPrompt->at(tank) > 0.) TOTHeCalPrompt += HeCal_E_tankMyPrompt->at(tank);
                        if(HeCal_E_tankMyPrompt->at(tank)>HeCal_MaxEMyprompt) {
                            HeCal_MaxEMyprompt = HeCal_E_tankMyPrompt->at(tank);
                            HeCal_MaxEMyprompt_CopyNo = tank;
                        }
                    }

                    Br_HeCal_E_tankMyPrompt->Fill();
                    Br_HeCal_MaxEMyprompt_CopyNo->Fill();
                    Br_HeCal_MaxEMyprompt->Fill();
                    Br_TOTHeCalPrompt->Fill();


                    // ---- HECAL capture delayed events ---- //
                    if(NextCaptureTime==0){ // no double capture
                        for(int j=0; j<HGasN; j++){
                            int copy = HGasCopyNo->at(j);
                            //getting the energy released in each tank
                            if(HGasE->at(j)>0. && HGasT->at(j) > gate ) HeCal_E_tankLate->at(copy) += HGasE->at(j);
                        }
                        for(int tank=0; tank<NTanks; tank++) {
                            //adding a gaussian smearing
                            HeCal_E_tankLate->at(tank) = r3->Gaus(HeCal_E_tankLate->at(tank), 0.1*HeCal_E_tankLate->at(tank));
                            if(HeCal_E_tankLate->at(tank) > 0.) TOTHeCalLate += HeCal_E_tankLate->at(tank);
                            if(HeCal_E_tankLate->at(tank)>HeCal_MaxELate){
                                HeCal_MaxELate = HeCal_E_tankLate->at(tank);
                                HeCal_MaxELate_CopyNo = tank;
                            }
                        }
                    }

                    else{ // double capture
                        for(int j=0; j<HGasN; j++){
                            int copy = HGasCopyNo->at(j);
                            if(HGasE->at(j)>0. && HGasT->at(j) > gate ) HeCal_E_tankLate->at(copy) += HGasE->at(j);
                        }
                        for(int tank=0; tank<NTanks; tank++) {
                            //adding a gaussian smearing
                            HeCal_E_tankLate->at(tank) = r3->Gaus(HeCal_E_tankLate->at(tank), 0.1*HeCal_E_tankLate->at(tank));
                            if(HeCal_E_tankLate->at(tank) > 0.) TOTHeCalLate += HeCal_E_tankLate->at(tank);
                            if(HeCal_E_tankLate->at(tank)>HeCal_MaxELate){
                                HeCal_MaxELate = HeCal_E_tankLate->at(tank);
                                HeCal_MaxELate_CopyNo = tank;
                            }
                        }
                    }

                    Br_HeCal_E_tankLate->Fill();
                    Br_HeCal_MaxELate_CopyNo->Fill();
                    Br_HeCal_MaxELate->Fill();
                    Br_TOTHeCalLate->Fill();       



                    // ---- HOD2 capture delayed events ---- //
                    if(NextCaptureTime==0){ // no double capture
                        for(int j=0; j<Hod2N; j++){
                            int copy2 = Hod2CopyNo->at(j);
                            if(Hod2T->at(j)>gate && Hod2E->at(j)!=0) {
                                Hod2_E_slabLate->at(copy2) += Hod2E->at(j);
                                Hod2Hit_Late++;
                            }
                        }

                        for(int hod=0; hod<NScint; hod++){
                            //adding a gaussian smearing
                            Hod2_E_slabLate->at(hod) = r3->Gaus(Hod2_E_slabLate->at(hod), 0.05*Hod2_E_slabLate->at(hod));
                            if(Hod2_E_slabLate->at(hod)>0 && Hod2_E_slabLate->at(hod)>Hod2_MaxELate) {
                                Hod2_MaxELate = Hod2_E_slabLate->at(hod);
                                Hod2_MaxELate_CopyNo = hod;
                            }
                        }
                    }

                    else{ // double capture
                        for(int j=0; j<Hod2N; j++){
                            int copy2 = Hod2CopyNo->at(j);
                            if(Hod2T->at(j)>gate && Hod2E->at(j)!=0) {
                                Hod2_E_slabLate->at(copy2) += Hod2E->at(j);
                                Hod2Hit_Late++;
                            }
                        }

                        for(int hod=0; hod<NScint; hod++){
                            //adding a gaussian smearing
                            Hod2_E_slabLate->at(hod) = r3->Gaus(Hod2_E_slabLate->at(hod), 0.05*Hod2_E_slabLate->at(hod));
                            if(Hod2_E_slabLate->at(hod)>0 && Hod2_E_slabLate->at(hod)>Hod2_MaxELate) {
                                Hod2_MaxELate = Hod2_E_slabLate->at(hod);
                                Hod2_MaxELate_CopyNo = hod;
                            }
                        }
                    }
                    
                    Br_Hod2_E_slabLate->Fill();
                    Br_Hod2_MaxELate_CopyNo->Fill();
                    Br_Hod2_MaxELate->Fill();
                    Br_Hod2Hit_Late->Fill();



                    // --- Beta evaluatuion --- //
                    if(Hod3N>0 && Hod2N>0){
                        double dist = sqrt(pow((Hod2Y->at(0)-Hod3Y->at(0)),2) + pow((Hod2X->at(0)-Hod3X->at(0)),2) + pow((Hod2Z->at(0)-Hod3Z->at(0)),2));
                        beta = (1./299.792458)*dist/(Hod2T->at(0)-Hod3T->at(0));  
                        //adding a gaussian smearing
                        beta = r3->Gaus(beta, 0.05*beta);
                    }
                    Br_beta->Fill();
                    
                    break; // ONLY THE FIRST CAPTURE

                } // end of the cycle on the first captures
            }

        }




        //PARTICLES
        else{

            //HOD3  
                int Hod3_MaxEMyprompt_CopyNo=0, Hod3Hit_Myprompt=0,
                    Hod3_MaxEMyLate_CopyNo=0, Hod3Hit_MyLate=0;
                double Hod3_MaxEMyprompt=0, Hod3_MaxEMyLate=0;
                double TOTHod3Prompt=0;

                branchName = "Hod3_MaxEMyprompt_CopyNo";
                TBranch *Br_Hod3_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&Hod3_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "Hod3_MaxEMyprompt";
                TBranch *Br_Hod3_MaxEMyprompt = t->Branch(branchName.Data(),&Hod3_MaxEMyprompt,branchName+"/D");
                branchName = "Hod3Hit_Myprompt";
                TBranch *Br_Hod3Hit_Myprompt = t->Branch(branchName.Data(),&Hod3Hit_Myprompt,branchName+"/I");
                branchName = "TOTHod3Prompt";
                TBranch *Br_TOTHod3Prompt = t->Branch(branchName.Data(),&TOTHod3Prompt,branchName+"/D");
                
                branchName = "Hod3_MaxEMyLate_CopyNo";
                TBranch *Br_Hod3_MaxEMyLate_CopyNo = t->Branch(branchName.Data(),&Hod3_MaxEMyLate_CopyNo,branchName+"/I");
                branchName = "Hod3_MaxEMyLate";
                TBranch *Br_Hod3_MaxEMyLate = t->Branch(branchName.Data(),&Hod3_MaxEMyLate,branchName+"/D");
                branchName = "Hod3Hit_MyLate";
                TBranch *Br_Hod3Hit_MyLate = t->Branch(branchName.Data(),&Hod3Hit_MyLate,branchName+"/I");
                
                //Energy divided by slab
                vector<double> *Hod3_E_slabMyPrompt=0, *Hod3_E_slabMyLate=0;

                branchName = "Hod3_E_slabMyPrompt";
                TBranch *Br_Hod3_E_slabMyPrompt = t->Branch(branchName.Data(),&Hod3_E_slabMyPrompt);
                branchName = "Hod3_E_slabMyLate";
                TBranch *Br_Hod3_E_slabMyLate = t->Branch(branchName.Data(),&Hod3_E_slabMyLate);

            //HOD2
                int Hod2_MaxEMyprompt_CopyNo=0, Hod2Hit_Myprompt=0,
                Hod2_MaxEMyLate_CopyNo=0, Hod2Hit_MyLate=0;
                double Hod2_MaxEMyprompt=0, Hod2_MaxEMyLate=0;
                double TOTHod2Prompt=0;

                branchName = "Hod2_MaxEMyprompt_CopyNo";
                TBranch *Br_Hod2_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&Hod2_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "Hod2_MaxEMyprompt";
                TBranch *Br_Hod2_MaxEMyprompt = t->Branch(branchName.Data(),&Hod2_MaxEMyprompt,branchName+"/D");
                branchName = "Hod2Hit_Myprompt";
                TBranch *Br_Hod2Hit_Myprompt = t->Branch(branchName.Data(),&Hod2Hit_Myprompt,branchName+"/I");
                branchName = "TOTHod2Prompt";
                TBranch *Br_TOTHod2Prompt = t->Branch(branchName.Data(),&TOTHod2Prompt,branchName+"/D");
                branchName = "Hod2_MaxEMyLate_CopyNo";
                TBranch *Br_Hod2_MaxEMyLate_CopyNo = t->Branch(branchName.Data(),&Hod2_MaxEMyLate_CopyNo,branchName+"/I");
                branchName = "Hod2_MaxEMyLate";
                TBranch *Br_Hod2_MaxEMyLate = t->Branch(branchName.Data(),&Hod2_MaxEMyLate,branchName+"/D");
                branchName = "Hod2Hit_MyLate";
                TBranch *Br_Hod2Hit_MyLate = t->Branch(branchName.Data(),&Hod2Hit_MyLate,branchName+"/I");

                //Energy divided by slab
                vector<double> *Hod2_E_slabMyPrompt=0, *Hod2_E_slabMyLate=0;
                branchName = "Hod2_E_slabMyPrompt";
                TBranch *Br_Hod2_E_slabMyPrompt = t->Branch(branchName.Data(),&Hod2_E_slabMyPrompt);
                branchName = "Hod2_E_slabMyLate";
                TBranch *Br_Hod2_E_slabMyLate = t->Branch(branchName.Data(),&Hod2_E_slabMyLate);            


            //HECAL reconstruction 
                int HeCal_MaxEMyprompt_CopyNo=0, HeCal_MaxEMyLate_CopyNo=0;
                double HeCal_MaxEMyprompt=0, HeCal_MaxEMyLate=0;
                double TOTHeCalPrompt=0, TOTHeCalLate=0; 

                branchName = "HeCal_MaxEMyprompt_CopyNo";
                TBranch *Br_HeCal_MaxEMyprompt_CopyNo = t->Branch(branchName.Data(),&HeCal_MaxEMyprompt_CopyNo,branchName+"/I");
                branchName = "HeCal_MaxEMyprompt";
                TBranch *Br_HeCal_MaxEMyprompt = t->Branch(branchName.Data(),&HeCal_MaxEMyprompt,branchName+"/D");
                branchName = "HeCal_MaxEMyLate_CopyNo";
                TBranch *Br_HeCal_MaxEMyLate_CopyNo = t->Branch(branchName.Data(),&HeCal_MaxEMyLate_CopyNo,branchName+"/I");
                branchName = "HeCal_MaxEMyLate";
                TBranch *Br_HeCal_MaxEMyLate = t->Branch(branchName.Data(),&HeCal_MaxEMyLate,branchName+"/D");
                branchName = "TOTHeCalPrompt";
                TBranch *Br_TOTHeCalPrompt = t->Branch(branchName.Data(),&TOTHeCalPrompt,branchName+"/D");
                branchName = "TOTHeCalLate";
                TBranch *Br_TOTHeCalLate = t->Branch(branchName.Data(),&TOTHeCalLate,branchName+"/D");

                //Energy divided by tank
                vector<double> *HeCal_E_tankMyPrompt=0, *HeCal_E_tankMyLate=0;

                branchName = "HeCal_E_tankMyPrompt";
                TBranch *Br_HeCal_E_tankMyPrompt = t->Branch(branchName.Data(),&HeCal_E_tankMyPrompt);
                branchName = "HeCal_E_tankMyLate";
                TBranch *Br_HeCal_E_tankMyLate = t->Branch(branchName.Data(),&HeCal_E_tankMyLate);

                //beta 
                double beta=0;
                branchName = "beta";
                TBranch *Br_beta = t->Branch(branchName.Data(),&beta,branchName+"/D");


            // ------------------- READING THE TREE ------------------------ //
            printf("Reading the tree ...\n");
            for (Int_t entry = 0; entry < t->GetEntries(); entry++){

                t->GetEntry(entry);
                gate = gate + Hod3T->at(0); //ns


                // --- INITIALIZATION --- 
                //HOD3
                Hod3_MaxEMyprompt_CopyNo=0; Hod3Hit_Myprompt=0;
                Hod3_MaxEMyLate_CopyNo=0; Hod3Hit_MyLate=0;
                Hod3_MaxEMyprompt=0; Hod3_MaxEMyLate=0;
                TOTHod3Prompt=0;

                Hod3_E_slabMyPrompt->resize(NScint, 0.);
                Hod3_E_slabMyLate->resize(NScint, 0.);

                for(int i=0; i<NScint; i++){
                    Hod3_E_slabMyPrompt->at(i)=0;
                    Hod3_E_slabMyLate->at(i)=0;
                }

                // ---- HOD3 prompt energy ---- //
                for(int j=0; j<Hod3N; j++){
                    int copy3 = Hod3CopyNo->at(j);
                    if(Hod3T->at(j)<gate && Hod3E->at(j)!=0){
                    Hod3_E_slabMyPrompt->at(copy3) += Hod3E->at(j);
                    Hod3Hit_Myprompt++;
                    }
                }

                for(int hod=0; hod<NScint; hod++){
                    //adding a gaussian smearing
                    Hod3_E_slabMyPrompt->at(hod) = r3->Gaus(Hod3_E_slabMyPrompt->at(hod), 0.05*Hod3_E_slabMyPrompt->at(hod));

                    if(Hod3_E_slabMyPrompt->at(hod)>0.) TOTHod3Prompt += Hod3_E_slabMyPrompt->at(hod);


                    if(Hod3_E_slabMyPrompt->at(hod)>0 && Hod3_E_slabMyPrompt->at(hod)>Hod3_MaxEMyprompt){
                    Hod3_MaxEMyprompt = Hod3_E_slabMyPrompt->at(hod);
                    Hod3_MaxEMyprompt_CopyNo = hod;
                    }
                }
                Br_Hod3_MaxEMyprompt_CopyNo->Fill();
                Br_Hod3_MaxEMyprompt->Fill();
                Br_Hod3Hit_Myprompt->Fill();
                Br_Hod3_E_slabMyPrompt->Fill();
                Br_TOTHod3Prompt->Fill();

                // ---- HOD3 late events ---- //
                for(int j=0; j<Hod3N; j++){
                    int copy3 = Hod3CopyNo->at(j);
                    if(Hod3T->at(j)>gate && Hod3E->at(j)!=0){
                    Hod3_E_slabMyLate->at(copy3) += Hod3E->at(j);
                    Hod3Hit_MyLate++;
                    }
                }
                for(int hod=0; hod<NScint; hod++){
                    //adding a gaussian smearing
                    Hod3_E_slabMyLate->at(hod) = r3->Gaus(Hod3_E_slabMyLate->at(hod), 0.05*Hod3_E_slabMyLate->at(hod));

                    if(Hod3_E_slabMyLate->at(hod)>0 && Hod3_E_slabMyLate->at(hod)>Hod3_MaxEMyLate) {
                    Hod3_MaxEMyLate = Hod3_E_slabMyLate->at(hod);
                    Hod3_MaxEMyLate_CopyNo = hod;
                    }
                }
                Br_Hod3_MaxEMyLate_CopyNo->Fill();
                Br_Hod3_MaxEMyLate->Fill();
                Br_Hod3Hit_MyLate->Fill();
                Br_Hod3_E_slabMyLate->Fill(); 


                //HOD2
                Hod2_MaxEMyprompt_CopyNo=0; Hod2Hit_Myprompt=0;
                Hod2_MaxEMyLate_CopyNo=0; Hod2Hit_MyLate=0;
                Hod2_MaxEMyprompt=0; Hod2_MaxEMyLate=0;
                TOTHod2Prompt=0;

                Hod2_E_slabMyPrompt->resize(NScint, 0.);
                Hod2_E_slabMyLate->resize(NScint, 0.);

                for(int i=0; i<NScint; i++){
                    Hod2_E_slabMyPrompt->at(i)=0;
                    Hod2_E_slabMyLate->at(i)=0;
                }

                // ---- HOD2 prompt energy ---- //
                for(int j=0; j<Hod2N; j++){
                    int copy2 = Hod2CopyNo->at(j);
                    if(Hod2T->at(j)<gate && Hod2E->at(j)!=0) {
                    Hod2_E_slabMyPrompt->at(copy2) += Hod2E->at(j);
                    Hod2Hit_Myprompt++;
                    }
                }
                for(int hod=0; hod<NScint; hod++){
                    //adding a gaussian smearing
                    Hod2_E_slabMyPrompt->at(hod) = r3->Gaus(Hod2_E_slabMyPrompt->at(hod), 0.05*Hod2_E_slabMyPrompt->at(hod));

                    if(Hod2_E_slabMyPrompt->at(hod)>0.) TOTHod2Prompt += Hod2_E_slabMyPrompt->at(hod);


                    if(Hod2_E_slabMyPrompt->at(hod)>0 && Hod2_E_slabMyPrompt->at(hod)>Hod2_MaxEMyprompt){
                    Hod2_MaxEMyprompt = Hod2_E_slabMyPrompt->at(hod);
                    Hod2_MaxEMyprompt_CopyNo = hod;
                    }
                }
                Br_Hod2_MaxEMyprompt_CopyNo->Fill();
                Br_Hod2_MaxEMyprompt->Fill();
                Br_Hod2Hit_Myprompt->Fill();
                Br_Hod2_E_slabMyPrompt->Fill();
                Br_TOTHod2Prompt->Fill();


                // ---- HOD2 late events ---- //
                for(int j=0; j<Hod2N; j++){
                    int copy2 = Hod2CopyNo->at(j);
                    if(Hod2T->at(j)>gate && Hod2E->at(j)!=0) {
                    Hod2_E_slabMyLate->at(copy2) += Hod2E->at(j);
                    Hod2Hit_MyLate++;
                    }
                }
                for(int hod=0; hod<NScint; hod++){
                    //adding a gaussian smearing
                    Hod2_E_slabMyLate->at(hod) = r3->Gaus(Hod2_E_slabMyLate->at(hod), 0.05*Hod2_E_slabMyLate->at(hod));

                    if(Hod2_E_slabMyLate->at(hod)>0 && Hod2_E_slabMyLate->at(hod)>Hod2_MaxEMyLate) {
                    Hod2_MaxEMyLate = Hod2_E_slabMyLate->at(hod);
                    Hod2_MaxEMyLate_CopyNo = hod;
                    }
                }
                Br_Hod2_MaxEMyLate_CopyNo->Fill();
                Br_Hod2_MaxEMyLate->Fill();
                Br_Hod2Hit_MyLate->Fill();
                Br_Hod2_E_slabMyLate->Fill();

                //HECAL
                HeCal_MaxEMyprompt_CopyNo=0; HeCal_MaxEMyLate_CopyNo=0;
                HeCal_MaxEMyprompt=0; HeCal_MaxEMyLate=0;
                TOTHeCalPrompt=0; TOTHeCalLate=0;            

                HeCal_E_tankMyPrompt->resize(NTanks, 0.);
                HeCal_E_tankMyLate->resize(NTanks, 0.);

                for(int i=0; i<NTanks; i++){
                    HeCal_E_tankMyPrompt->at(i)=0;
                    HeCal_E_tankMyLate->at(i)=0;
                }

                // ---- HECAL prompt energy ---- //
                for(int j=0; j<HGasN; j++){
                    int copy = HGasCopyNo->at(j);
                    if(HGasE->at(j)>0. && HGasT->at(j)<gate) HeCal_E_tankMyPrompt->at(copy) += HGasE->at(j);
                }

                for(int tank=0; tank<NTanks; tank++) {
                    //adding a gaussian smearing
                    HeCal_E_tankMyPrompt->at(tank) = r3->Gaus(HeCal_E_tankMyPrompt->at(tank), 0.1*HeCal_E_tankMyPrompt->at(tank));
                    if(HeCal_E_tankMyPrompt->at(tank)>0.) TOTHeCalPrompt += HeCal_E_tankMyPrompt->at(tank);

                    if(HeCal_E_tankMyPrompt->at(tank)>HeCal_MaxEMyprompt) {
                        HeCal_MaxEMyprompt = HeCal_E_tankMyPrompt->at(tank);
                        HeCal_MaxEMyprompt_CopyNo = tank;
                    }
                }
                Br_HeCal_MaxEMyprompt_CopyNo->Fill();
                Br_HeCal_MaxEMyprompt->Fill();
                Br_HeCal_E_tankMyPrompt->Fill();
                Br_TOTHeCalPrompt->Fill();


                // ---- HECAL late events ---- //
                
                for(int j=0; j<HGasN; j++){
                    int copy = HGasCopyNo->at(j);
                    if(HGasE->at(j)>0. && HGasT->at(j)>gate) HeCal_E_tankMyLate->at(copy) += HGasE->at(j);
                }

                for(int tank=0; tank<NTanks; tank++) {
                    //adding a gaussian smearing
                    HeCal_E_tankMyLate->at(tank) = r3->Gaus(HeCal_E_tankMyLate->at(tank), 0.1*HeCal_E_tankMyLate->at(tank));
                    if(HeCal_E_tankMyLate->at(tank)>0.) TOTHeCalLate += HeCal_E_tankMyLate->at(tank);

                    if(HeCal_E_tankMyLate->at(tank)>HeCal_MaxEMyLate){
                        HeCal_MaxEMyLate = HeCal_E_tankMyLate->at(tank) ;
                        HeCal_MaxEMyLate_CopyNo = tank;
                    }
                }
                Br_HeCal_MaxEMyLate_CopyNo->Fill();
                Br_HeCal_MaxEMyLate->Fill();
                Br_HeCal_E_tankMyLate->Fill();
                Br_TOTHeCalLate->Fill();


                // --- Beta evaluatuion --- //

                if(Hod3N>0 && Hod2N>0){
                    double dist = sqrt(pow((Hod2Y->at(0)-Hod3Y->at(0)),2) + pow((Hod2X->at(0)-Hod3X->at(0)),2) + pow((Hod2Z->at(0)-Hod3Z->at(0)),2));
                    beta = (1./299.792458)*dist/(Hod2T->at(0)-Hod3T->at(0));  
                    //adding a gaussian smearing
                    beta = r3->Gaus(beta, 0.05*beta);
                }
                Br_beta->Fill();

                
                // evaluate incidence on ToF shell (F.Nozzoli)
                // Cosi[0]=1;
                // Cosi[1]=1;
                // Cosi[2]=1; // cosine WRT outher shell
                // double dir[3];
                // if(Hod3N>0 && Hod2N>0){
                //   dir[0] = Hod2X->at(0)-Hod3X->at(0);
                //   dir[1] = Hod2Y->at(0)-Hod3Y->at(0);
                //   dir[2] = Hod2Z->at(0)-Hod3Z->at(0); 
                // }
                // else {dir[0]=1;dir[1]=1;dir[2]=1;}
                // double len = sqrt(pow(dir[0],2) + pow(dir[1],2) + pow(dir[2],2)); 
                // for(int jj = 0;jj<3;jj++) dir[jj]=dir[jj]/len;
                // if(Hod3N>0){
                //   if(fabs(fabs(Hod3Z->at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[2]); // select SIDE
                //   if(fabs(fabs(Hod3Y->at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[1]);
                //   if(fabs(fabs(Hod3X->at(0))-MD->MyRS3)<0.5) Cosi[2]=fabs(dir[0]);
                // }
                // if(Hod2N>0){
                //   if(fabs(fabs(Hod2Z->at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[2]); // select SIDE
                //   if(fabs(fabs(Hod2Y->at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[1]);
                //   if(fabs(fabs(Hod2X->at(0))-MD->MyRS2)<0.5) Cosi[1]=fabs(dir[0]);
                // }
            } 

        } //END PARTICLES RECONSTRUCTION

        t->Write();
        outputFile->Write();
        
    } //CYCLE ON THE FILES
}
