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

void ADHD_doubleCapture()
{
    // gStyle->SetOptStat(0); // no statistics

    // ----------------- INPUT FILE ----------------- //
    const int nFiles = 1;
    // const int nFiles = 7;

    TString fInNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_STAT_gps_1raw.root"};

    
    // TString fInNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_raw.root", "../../sources/ROOT_files/pbar_10e7_gps_raw.root",
    //                           "../../sources/ROOT_files/proton_10e7_gps_rec.root", "../../sources/ROOT_files/Deuteron_10e7_gps_rec.root",
    //                           "../../sources/ROOT_files/He4_10e7_gps_rec.root", "../../sources/ROOT_files/C12_10e7_gps_rec.root",
    //                           "../../sources/ROOT_files/electron_10e7_gps_rec.root"};



    // Gate intervals
    const double Cgate = 50; // ns
    double gate = 0;
    TString Sgate = to_string(int(Cgate));


    // ----------------- OUTPUT FILE ----------------- //
    // TString fOutNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e6_gps_SM1rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM1rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM2rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM2rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM3rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM3rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM4rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM4rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM5rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM5rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM6rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM6rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM7rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM7rec_.root",
    //                              "../../sources/ROOT_files/Dbar_10e6_gps_SM8rec_.root", "../../sources/ROOT_files/pbar_10e6_gps_SM8rec_.root",
    //                              "../../sources/ROOT_files/proton_10e7_gps_SMrec_.root", "../../sources/ROOT_files/Deuteron_10e7_gps_SMrec_.root"};

    // TString fOutNames[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_rec_.root", "../../sources/ROOT_files/pbar_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/proton_10e7_gps_rec_.root", "../../sources/ROOT_files/Deuteron_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/He4_10e7_gps_rec_.root", "../../sources/ROOT_files/C12_10e7_gps_rec_.root",
    //                             "../../sources/ROOT_files/electron_10e7_gps_rec_.root"};

    // for(int i=0; i<nFiles; i++) fOutNames[i].ReplaceAll(".root", "_"+Sgate+"ns.root");

    TString branchName;

    // Constants
    const double MyDelay = 900000; // ns
    //Random generator
    TRandom3 *r3 = new TRandom3(0);


    TH1D *h1_mass1 = new TH1D("h1_mass1", "Particle captured mass before double capture; Mass [MeV/c^{2}]; Entries", 1100, 900, 1900);
    h1_mass1->SetLineWidth(3);
    h1_mass1->SetLineColor(kRed);
    TH1D *h1_mass2 = new TH1D("h1_mass2", "Particle captured mass in a double capture; Mass [MeV/c^{2}]; Entries", 100, 900, 1000);
    h1_mass2->SetLineWidth(3);


    // CYCLE ON THE FILES
    for (int i = 0; i < nFiles; i++){

        TFile *f = new TFile(fInNames[i], "OPEN");
        if (!f){
            printf("File %s not found\n", fInNames[i].Data());
            return;
        }

         // --- getting the tree and setting specific branches --- //
        //antiparticles
        // TTree *tr;
        TTree *t ;

        // if(fInNames[i].Contains("bar")) tr = (TTree *)f->Get("Hits_bar");
        if(fInNames[i].Contains("bar")) t = (TTree *)f->Get("Hits_bar");
        //particle
        else{
            std::printf("\nNot able to run on ordinary particle files: %s \n", fInNames[i].Data());
            return;
        }

        // std::printf("\nCopying tree from file %s \n", fInNames[i].Data());

        // gROOT->cd();
        // TFile* outputFile = TFile::Open(fOutNames[i], "RECREATE");
        // outputFile->cd();
        // TTree *t = (TTree*) tr->CloneTree();
        
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

        //HeCal
        //HeCal reconstruction for captures
        int HeCal_MaxEprompt_CopyNo=0, HeCal_MaxEAcap_CopyNo=0, HeCal_MaxEDelayed_CopyNo=0;
        double HeCal_MaxEprompt=0, HeCal_MaxEAcap=0, HeCal_MaxEDelayed=0;
   

        //Energy released for each HeCal tank
        vector<double> *HeCal_E_tankPrompt=0, *HeCal_E_tankAcap=0, *HeCal_E_tankDelayed=0;
        //Energy released for each HeCal tank WITHOUT captures
        vector<double> *HeCal_E_tankMyPrompt=0, *HeCal_E_tankLate=0;


        //HeCal reconstruction without captures
        int HeCal_MaxEMyprompt_CopyNo=0, HeCal_MaxELate_CopyNo=0;
        double HeCal_MaxEMyprompt=0, HeCal_MaxELate=0;
        double TOTHeCalPrompt=0, TOTHeCalLate=0; 


        //beta 
        double beta=0;
        branchName = "beta";


        // ------------------- READING THE TREE ------------------------ //
        printf("Reading the tree ...\n");
        for (Int_t entry = 0; entry < t->GetEntries(); entry++){

            t->GetEntry(entry);

            gate = Cgate + Hod3T->at(0); // ns


            // --- INITIALIZATION --- 


            //HeCal reconstruction for captures
            HeCal_MaxEprompt_CopyNo=0; HeCal_MaxEAcap_CopyNo=0; HeCal_MaxEDelayed_CopyNo=0;
            HeCal_MaxEprompt=0; HeCal_MaxEAcap=0; HeCal_MaxEDelayed=0;
                //reconstruction without captures
            HeCal_MaxEMyprompt_CopyNo=0; HeCal_MaxELate_CopyNo=0;
            HeCal_MaxEMyprompt=0; HeCal_MaxELate=0;
            TOTHeCalPrompt=0; TOTHeCalLate=0;

            //number of slabs = 64*6
            int NScint = 64*6;
            

            //number of tanks = 75
            int NTanks = 75;


            // --- CYCLE ON THE CAPTURES --- //
            // entrance in the for cycle is granted by the CaptureT size that is always at least 1
            // In the absence of a capture the value inside CaptureT is 0
            // The reconstruction for anti particles is still made here, even if no capture happens
            
            for(int k=0; k < int(CaptureT->size()) ; k++){
                
                double time = CaptureT->at(k);  // = 0 if no capture
                if(time==0 && CaptureT->size()>1) continue; // if one capture happens, skip the first entry in the vector
                double NextCaptureTime = 0;
                if( k+1 < int(CaptureT->size()) && time!=0 ) NextCaptureTime = CaptureT->at(k+1);

                //HeCal reconstruction for captures
                HeCal_MaxEprompt_CopyNo=0; HeCal_MaxEAcap_CopyNo=0; HeCal_MaxEDelayed_CopyNo=0;
                HeCal_MaxEprompt=0; HeCal_MaxEAcap=0; HeCal_MaxEDelayed=0;
                    //reconstruction without captures
                HeCal_MaxEMyprompt_CopyNo=0; HeCal_MaxELate_CopyNo=0;
                HeCal_MaxEMyprompt=0; HeCal_MaxELate=0;
                TOTHeCalPrompt=0; TOTHeCalLate=0;


                // ---- HECAL ---- //

            
                // capture delayed events

                if(NextCaptureTime!=0){ // double capture

                    h1_mass1->Fill(CaptureM->at(1));
                    h1_mass2->Fill(CaptureM->at(2));
                }

                
                break; // ONLY THE FIRST CAPTURE

            } // end of the cycle on the first captures
        }

        h1_mass1->Draw();
        TCanvas *c2 = new TCanvas("c2","c2", 800, 800);
        h1_mass2->Draw();

    

        // t->Write();
        // outputFile->Write();
        
    } //CYCLE ON THE FILES
}
