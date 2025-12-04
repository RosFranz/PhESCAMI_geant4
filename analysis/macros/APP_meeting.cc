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
Programma che richiede il valore del gate in ingresso (defualt e' 50 ns).
Legge in automatico dei file predefiniti e fornisce dei confronti tra diverse
distribuzioni; i.e. max E in HeCal vs MCEnergy
                    Max E (hod3)   vs MCEnergy
                    Max E (hod2)   vs MCEnergy
                   #hit (S1+S2)    vs MCEnergy.
DA USARE per file CON DIVERSE SPECIE DI PARTICELLE.
*/


// function to obtain log bin along X and Y (optional)
// axis of a TH2D histogram


TH2D* GenerateHistogram(TString Name, TString title, TString XTitle, TString YTitle, 
                        int NbinsX, double xmin, double xmax, int NbinsY, double ymin, double ymax, bool ylog)
{
    std::vector<double> xBins;
    std::vector<double> xBinsLog;

    int NEdgesX = NbinsX + 1;
    double xMaxLog = TMath::Log10(xmax);
    double xMinLog = TMath::Log10(xmin);

    double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;
    
    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog);
        xBins.push_back(TMath::Power(10, xBinsLog.at(i)));
    }

    std::vector<double> yBins;
    std::vector<double> yBinsLog;

    int NEdgesY = NbinsY + 1;
    double yMaxLog = TMath::Log10(ymax);
    double yMinLog = TMath::Log10(ymin);

    double yBinWidthLog = (yMaxLog - yMinLog) / NbinsY;
    double yBinWidth = (ymax - ymin) / NbinsY;

    for (int i = 0; i < NEdgesY; i++)
    {
        yBinsLog.push_back(yMinLog + i * yBinWidthLog);
        if(ylog) yBins.push_back(TMath::Power(10, yBinsLog.at(i)));
        else yBins.push_back(ymin + i * yBinWidth);
    }

    TString Title_tot =  title+";"+XTitle+";"+YTitle;

    TH2D *histogram = new TH2D(Name, Title_tot, NbinsX, &xBins[0], NbinsY, &yBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);

    return histogram;
}





//!!!!!!!!!!!!!!
// MAIN PROGRAM
//!!!!!!!!!!!!!!

void APP_meeting(int gate)
{
    gStyle->SetOptStat(0); // no statistics

    // ----------------- INPUT FILE ----------------- //
    const int nFiles = 5;
    Color_t colors[nFiles] = {kBlue, kRed, kRed - 3, kBlue - 5};
    //bugged on when starts the time counter
    // TString fNames_50[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_rec_corrected.root", "../../sources/ROOT_files/pbar_10e7_rec.root",
    //                           "../../sources/ROOT_files/proton_10e7_rec.root", "../../sources/ROOT_files/Deuteron_10e7_rec.root",
    //                           "../../sources/ROOT_files/He4_10e7_rec.root"};

    TString fNames_40[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_rec__40ns.root", "../../sources/ROOT_files/pbar_10e7_rec__40ns.root",
                              "../../sources/ROOT_files/proton_10e7_rec__40ns.root", "../../sources/ROOT_files/Deuteron_10e7_rec__40ns.root",
                              "../../sources/ROOT_files/He4_10e7_rec__40ns.root"};
    
    TString fNames_50[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/pbar_10e7_gps_rec__50ns.root",
                              "../../sources/ROOT_files/proton_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/Deuteron_10e7_gps_rec__50ns.root",
                              "../../sources/ROOT_files/He4_10e7_gps_rec__50ns.root"};
    
    TString fNames_60[nFiles] = {"../../sources/ROOT_files/Dbar_10e7_rec__60ns.root", "../../sources/ROOT_files/pbar_10e7_rec__60ns.root",
                              "../../sources/ROOT_files/proton_10e7_rec__60ns.root", "../../sources/ROOT_files/Deuteron_10e7_rec__60ns.root",
                              "../../sources/ROOT_files/He4_10e7_rec__60ns.root"};
    TString fNames[nFiles];
    TString outFileName;

    switch (gate){
        case 40:
            for(int i=0; i<nFiles; i++) fNames[i] = fNames_40[i];
            outFileName = "APP_meeting_40.root";
            printf("You selected a 40  ns gate\n");
            break;
        case 50:{
            for(int i=0; i<nFiles; i++) fNames[i] = fNames_50[i];
            outFileName = "APP_meeting_50.root";
            printf("You selected a 50  ns gate\n");
            break;
        }
        case 60:{
            for(int i=0; i<nFiles; i++) fNames[i] = fNames_60[i];
            outFileName = "APP_meeting_60.root";
            printf("You selected a 60  ns gate\n");
            break;
        }
        // in case the variable gate assumes another value, rise an error flag
        default:{
            printf("ERROR: gate value not valid\n Currently available values are 40, 50 and 60 ns\n");
            return;
        }
    }
    
    TString hNames[nFiles] = {"#bar{d}", "#bar{p}", "p", "d", "He4"};

    const int nPlot = 4;

    TString path = "/../../analysis/macros/";

    

    // ---- Histograms ---- //
    // MC TRUTH
    TH2D *h2_HeVSMCEnergy[nFiles];    // Max prompt energy in HeCal                               vs MCEnergy
    TH2D *h2_HeDelVSMCEnergy[nFiles]; // Delayed maximum energy release in HeCal                  vs MCEnergy

    TH2D *h2_Hod3VSMCEnergy[nFiles];    // Max prompt energy in the Hod3 (external scintillator)  vs MCEnergy
    TH2D *h2_Hod3DelVSMCEnergy[nFiles]; // Max delayed energy in the Hod3 (external scintillator) vs MCEnergy

    TH2D *h2_Hod2VSMCEnergy[nFiles];    // Max prompt energy in the Hod3 (external scintillator)  vs MCEnergy
    TH2D *h2_Hod2DelVSMCEnergy[nFiles]; // Max delayed energy in the Hod3 (external scintillator) vs MCEnergy

    TH2D *h2_Hod3NVSMCEnergy[nFiles];    // Number of HOD3 slab above PROMPT energy threshold     vs MCEnergy
    TH2D *h2_Hod3NDelVSMCEnergy[nFiles]; // Number of HOD3 slab above DLAYED energy threshold     vs MCEnergy

    TH2D *h2_Hod2NVSMCEnergy[nFiles];    // Number of HOD2 slab above PROMPT energy threshold     vs MCEnergy
    TH2D *h2_Hod2NDelVSMCEnergy[nFiles]; // Number of HOD2 slab above DLAYED energy threshold     vs MCEnergy


    // Thresholds
    const double thr_HeCal = 10.;    //Prompt and Delayed HeCal [MeV]
    const double thr_HodPrompt = 1.6;    // Prompt Scintillators [MeV]
    const double thr_HodDelayed = 0.8;   // Delayed Scintillators [MeV]
    const double thr_HodHitPrompt = 2.; // Prompt Hits in scintillators
    const double thr_HodHitDelayed = 2.;// Delayed Hits in scintillators

    printf("Thresholds:\n");
    printf("\t HeCal prompt > %.2f MeV\n", thr_HeCal);
    printf("\t HeCal delayed > %.2f MeV\n", thr_HeCal);
    printf("\t Scintillators prompt > %.2f MeV\n", thr_HodPrompt);
    printf("\t Scintillators delayed > %.2f MeV\n", thr_HodDelayed);
    printf("\t Hits prompt < %.2f\n", thr_HodHitPrompt);
    printf("\t Hits delayed > %.2f\n", thr_HodHitDelayed);


    // Constants
    const double MyDelay = 900000; // ns

    // CYCLE ON THE FILES
    for (int i = 0; i < nFiles; i++)
    {
        TFile *f = new TFile(fNames[i], "OPEN");
        if (!f){
            printf("File %s not found\n", fNames[i].Data());
            return;
        }

        // --- initializing histograms --- //
        h2_HeVSMCEnergy[i] = GenerateHistogram("HeCal_Max_MaxEMyprompt_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (HeCal prompt)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_HeVSMCEnergy[i]->Sumw2();

        h2_HeDelVSMCEnergy[i] = GenerateHistogram("HeCal_Max_MaxELate_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (HeCal delayed)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_HeDelVSMCEnergy[i]->Sumw2();

        h2_Hod3VSMCEnergy[i] = GenerateHistogram("Hod3_Max_MaxEMyprompt_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (Hod3 prompt)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_Hod3VSMCEnergy[i]->Sumw2();

        h2_Hod3DelVSMCEnergy[i] = GenerateHistogram("Hod3_Max_MaxELate_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (Hod3 delayed)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_Hod3DelVSMCEnergy[i]->Sumw2();

        h2_Hod2VSMCEnergy[i] = GenerateHistogram("Hod2_Max_MaxEMyprompt_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (Hod2 prompt)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_Hod2VSMCEnergy[i]->Sumw2();

        h2_Hod2DelVSMCEnergy[i] = GenerateHistogram("Hod2_Max_MaxELate_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "E_{max} (Hod2 delayed)", 400, 10, 10000, 100, 0.1, 200.1, true);
        h2_Hod2DelVSMCEnergy[i]->Sumw2();

        h2_Hod3NVSMCEnergy[i] = GenerateHistogram("N_Hod3_slab_prompt_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "N slab prompt Hod3", 400, 10, 10000, 20, -0.5, 19.5, false);
        h2_Hod3NVSMCEnergy[i]->Sumw2();

        h2_Hod2NVSMCEnergy[i] = GenerateHistogram("N_Hod2_slab_prompt_vs_" + hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "N slab prompt Hod2", 400, 10, 10000, 20, -0.5, 19.5, false);
        h2_Hod2NVSMCEnergy[i]->Sumw2();

        h2_Hod3NDelVSMCEnergy[i] = GenerateHistogram("N_Hod3_slab_late_vs_"+ hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "N slab late Hod3", 400, 10, 10000, 20, -0.5, 19.5, false);
        h2_Hod3NDelVSMCEnergy[i]->Sumw2();

        h2_Hod2NDelVSMCEnergy[i] = GenerateHistogram("N_Hod2_slab_late_vs_"+ hNames[i] + "_energy", hNames[i], "E_{k} [MeV]/N", "N slab late Hod2", 400, 10, 10000, 20, -0.5, 19.5, false);
        h2_Hod2NDelVSMCEnergy[i]->Sumw2();


        //variable for the branches
        //
        vector<double> *CaptureT=0;
        TBranch *b_CaptureT=0;

        // HeCal
        double HeCal_MaxEMyprompt, HeCal_MaxELate;
        // double HeCal_Acap_E;

        // Plastic scintillators
        double Hod3_MaxEMyprompt, Hod3_MaxELate;
        // double Hod3_Acap_E
        double Hod2_MaxEMyprompt, Hod2_MaxELate;
        // double Hod2_Acap_E;

        int Hod3Hit_Myprompt,  Hod2Hit_Myprompt,
            // Hod3Hit_Acap,    Hod2Hit_Acap,
            Hod3Hit_Late, Hod2Hit_Late;

        //Definition of the applied cuts
        TString cuts_bar = "HeCal_MaxEMyprompt_CopyNo==HeCal_MaxELate_CopyNo"; // antiparticles
        TString cuts = "HeCal_MaxEprompt_CopyNo==HeCal_MaxELate_CopyNo";       // particles


        // --- getting the tree and setting specific branches --- //
        if(fNames[i].Contains("bar")) { //antiparticles

            printf("\n------------- ANTIPARTICLE ------------- \n");

            //
            // NO SELECTION BEFORE TRIGGER DISTRIBUTIONS
            //
            // TTree *tr = (TTree *)f->Get("RHits_bar");
            // std::printf("\nCopying tree from file %s \n", fNames[i].Data());

            // //The problem is that you are copying the TTree into the (read-only) input file.
            // //If you intend on copying the TTree in memory you can do:
            // gROOT->cd();

            // TFile* outputFile = TFile::Open("CopiedTree.root", "RECREATE");
            // TTree *t = (TTree*) tr->CopyTree(cuts_bar.Data());

            TTree *t;
            if(fNames[i].Contains("__")) t = (TTree *)f->Get("Hits_bar");
            else t = (TTree *)f->Get("RHits_bar");

            // HeCal
            // t->SetBranchAddress("HeCal_Acap_E", &HeCal_Acap_E);
            t->SetBranchAddress("HeCal_MaxELate", &HeCal_MaxELate);

            // Plastic scintillators
            // t->SetBranchAddress("Hod3_Acap_E", &Hod3_Acap_E);
            // t->SetBranchAddress("Hod3Hit_Acap", &Hod3Hit_Acap);

            // t->SetBranchAddress("Hod2_Acap_E", &Hod2_Acap_E);
            // t->SetBranchAddress("Hod2Hit_Acap", &Hod2Hit_Acap);

            t->SetBranchAddress("Hod3_MaxELate", &Hod3_MaxELate);
            t->SetBranchAddress("Hod3Hit_Late", &Hod3Hit_Late); 

            t->SetBranchAddress("Hod2_MaxELate", &Hod2_MaxELate);
            t->SetBranchAddress("Hod2Hit_Late", &Hod2Hit_Late); 

            //clusterizzazione degli scintillatori: HOD3 e HOD2
            vector<double> *Hod3_E_slabMyPrompt=0, *Hod2_E_slabMyPrompt=0;
            TBranch *b_Hod3_E_slabMyPrompt=0, *b_Hod2_E_slabMyPrompt=0;
            t->SetBranchAddress("Hod3_E_slabMyPrompt",&Hod3_E_slabMyPrompt,&b_Hod3_E_slabMyPrompt);
            t->SetBranchAddress("Hod2_E_slabMyPrompt",&Hod2_E_slabMyPrompt,&b_Hod2_E_slabMyPrompt);

            vector<double> *Hod3_E_slabLate=0, *Hod2_E_slabLate=0;
            TBranch *b_Hod3_E_slabLate=0, *b_Hod2_E_slabLate=0;
            t->SetBranchAddress("Hod3_E_slabLate",&Hod3_E_slabLate,&b_Hod3_E_slabLate);
            t->SetBranchAddress("Hod2_E_slabLate",&Hod2_E_slabLate,&b_Hod2_E_slabLate);



            Long64_t entries = t->GetEntries();
            Long64_t captures = t->GetEntries("CaptureT>0 && CaptureT<9000");
            std::printf("\nReading file %s with %lld entries\n", fNames[i].Data(), entries);
            std::printf("Number of captures: %lld\n", captures);

            // setting other common branches
            t->SetBranchAddress("CaptureT", &CaptureT, &b_CaptureT);
                // HeCal
            t->SetBranchAddress("HeCal_MaxEMyprompt", &HeCal_MaxEMyprompt);

                // Plastic scintillators
            t->SetBranchAddress("Hod3_MaxEMyprompt", &Hod3_MaxEMyprompt);
            t->SetBranchAddress("Hod3Hit_Myprompt", &Hod3Hit_Myprompt);

            t->SetBranchAddress("Hod2_MaxEMyprompt", &Hod2_MaxEMyprompt);
            t->SetBranchAddress("Hod2Hit_Myprompt", &Hod2Hit_Myprompt);


                // MC truth
            double MCEnergy;
            int HGasN, Hod2N, Hod3N, EventID;

            t->SetBranchAddress("MCEnergy", &MCEnergy);
            t->SetBranchAddress("HGasN", &HGasN);
            t->SetBranchAddress("Hod2N", &Hod2N);
            t->SetBranchAddress("Hod3N", &Hod3N);
            t->SetBranchAddress("EventID", &EventID);

            // First beta measurements
            double beta;
            t->SetBranchAddress("beta", &beta);


            // ------------------- READING THE TREE ------------------------ //
            Long64_t counting = 0;
            int NHod3Prompt = 0, NHod3Late = 0,
                NHod2Prompt = 0, NHod2Late = 0;

            for (Int_t entry = 0; entry < entries; entry++){

                NHod3Prompt = 0; NHod3Late = 0;
                NHod2Prompt = 0; NHod2Late = 0;

                t->GetEntry(entry);

                //counting the number of scintillators above the PROMPT threshold
                // and the number of scintillators between [high, low] thresholds
                for(int i=0; i<Hod3_E_slabMyPrompt->size(); i++) {
                    if(Hod3_E_slabMyPrompt->at(i)>thr_HodPrompt) NHod3Prompt++;
                    if(Hod2_E_slabMyPrompt->at(i)>thr_HodPrompt) NHod2Prompt++;

                    if(thr_HodPrompt > Hod3_E_slabLate->at(i) && Hod3_E_slabLate->at(i) > thr_HodDelayed) NHod3Late++;
                    if(thr_HodPrompt > Hod2_E_slabLate->at(i) && Hod2_E_slabLate->at(i) > thr_HodDelayed) NHod2Late++;
                }

                if( (Hod3_MaxELate > thr_HodDelayed) && 
                    (Hod2_MaxELate > thr_HodDelayed) && 
                    (HeCal_MaxELate > thr_HeCal) &&
                    (NHod3Late > thr_HodHitDelayed) &&
                    (NHod2Late > thr_HodHitDelayed) 
                    ){
                        //filling the prompt histograms
                        h2_HeVSMCEnergy[i]->Fill(MCEnergy, HeCal_MaxEMyprompt);

                        h2_Hod3VSMCEnergy[i]->Fill(MCEnergy, Hod3_MaxEMyprompt);
                        h2_Hod2VSMCEnergy[i]->Fill(MCEnergy, Hod2_MaxEMyprompt);

                        h2_Hod3NVSMCEnergy[i]->Fill(MCEnergy, NHod3Prompt);
                        h2_Hod2NVSMCEnergy[i]->Fill(MCEnergy, NHod2Prompt);

                    }

                


                //filling the delayed histograms
                h2_Hod3NDelVSMCEnergy[i]->Fill(MCEnergy, NHod3Late);
                h2_Hod2NDelVSMCEnergy[i]->Fill(MCEnergy, NHod2Late);

                h2_HeDelVSMCEnergy[i]->Fill(MCEnergy, HeCal_MaxELate);

                h2_Hod3DelVSMCEnergy[i]->Fill(MCEnergy, Hod3_MaxELate);

                h2_Hod2DelVSMCEnergy[i]->Fill(MCEnergy, Hod2_MaxELate);
            }

            printf("Number of fillments: %.1f\n", h2_HeVSMCEnergy[i]->GetEntries());

            printf("\nNumber of captures: %lld\n", captures);
        }










        else {

            printf("\n------------- ORDINARY MATTER ------------- \n");


            TTree *t = (TTree *)f->Get("RHits");
            Long64_t entries = t->GetEntries();
            std::printf("\nReading file %s with %lld entries\n", fNames[i].Data(), entries);

            //clusterizzazione degli scintillatori: HOD3 e HOD2
            vector<double> *Hod3_E_slabMyPrompt=0, *Hod2_E_slabMyPrompt=0;
            TBranch *b_Hod3_E_slabMyPrompt=0, *b_Hod2_E_slabMyPrompt=0;
            vector<double> *Hod3_E_slabLate=0, *Hod2_E_slabLate=0;
            TBranch *b_Hod3_E_slabLate=0, *b_Hod2_E_slabLate=0;
            // MC truth
            double MCEnergy;
            int HGasN, Hod2N, Hod3N, EventID;
            // First beta measurements
            double beta;


            if(fNames[i].Contains("__")) { // reconstruction different gate respect to the default 50ns
                // HeCal
                t->SetBranchAddress("HeCal_MaxEMyLate", &HeCal_MaxELate);
                t->SetBranchAddress("HeCal_MaxEMyprompt", &HeCal_MaxEMyprompt);

                // Plastic scintillators
                t->SetBranchAddress("Hod3_MaxEMyLate", &Hod3_MaxELate);
                t->SetBranchAddress("Hod3Hit_MyLate", &Hod3Hit_Late);

                t->SetBranchAddress("Hod2_MaxEMyLate", &Hod2_MaxELate);
                t->SetBranchAddress("Hod2Hit_MyLate", &Hod2Hit_Late);

                
                t->SetBranchAddress("Hod3_E_slabMyPrompt",&Hod3_E_slabMyPrompt,&b_Hod3_E_slabMyPrompt);
                t->SetBranchAddress("Hod2_E_slabMyPrompt",&Hod2_E_slabMyPrompt,&b_Hod2_E_slabMyPrompt);

                t->SetBranchAddress("Hod3_E_slabMyLate",&Hod3_E_slabLate,&b_Hod3_E_slabLate);
                t->SetBranchAddress("Hod2_E_slabMyLate",&Hod2_E_slabLate,&b_Hod2_E_slabLate);

                // setting other common branches
                    // Plastic scintillators
                t->SetBranchAddress("Hod3_MaxEMyprompt", &Hod3_MaxEMyprompt);
                t->SetBranchAddress("Hod3Hit_Myprompt", &Hod3Hit_Myprompt);

                t->SetBranchAddress("Hod2_MaxEMyprompt", &Hod2_MaxEMyprompt);
                t->SetBranchAddress("Hod2Hit_Myprompt", &Hod2Hit_Myprompt);

                //MC truth
                t->SetBranchAddress("MCEnergy", &MCEnergy);
                t->SetBranchAddress("HGasN", &HGasN);
                t->SetBranchAddress("Hod2N", &Hod2N);
                t->SetBranchAddress("Hod3N", &Hod3N);
                t->SetBranchAddress("EventID", &EventID);
                
                t->SetBranchAddress("beta", &beta);
            }
            else{

                //
                // NO SELECTION BEFORE TRIGGER DISTRIBUTIONS
                //

                // TTree *tr = (TTree *)f->Get("RHits");

                // std::printf("\nCopying tree from file %s \n", fNames[i].Data());


                // //The problem is that you are copying the TTree into the (read-only) input file.
                // //If you intend on copying the TTree in memory you can do:
                // gROOT->cd();

                // TFile* outputFile = TFile::Open("CopiedTree.root", "RECREATE");

                // TTree *t = (TTree*) tr->CopyTree(cuts.Data());


                // HeCal
                t->SetBranchAddress("HeCal_MaxELate", &HeCal_MaxELate); 
                t->SetBranchAddress("HeCal_MaxEprompt", &HeCal_MaxEMyprompt);

                // Plastic scintillators
                t->SetBranchAddress("Hod3_MaxELate", &Hod3_MaxELate);
                t->SetBranchAddress("Hod3Hit_Late", &Hod3Hit_Late);

                t->SetBranchAddress("Hod2_MaxELate", &Hod2_MaxELate);
                t->SetBranchAddress("Hod2Hit_Late", &Hod2Hit_Late);

                
                t->SetBranchAddress("Hod3_E_slabPrompt",&Hod3_E_slabMyPrompt,&b_Hod3_E_slabMyPrompt);
                t->SetBranchAddress("Hod2_E_slabPrompt",&Hod2_E_slabMyPrompt,&b_Hod2_E_slabMyPrompt);

                t->SetBranchAddress("Hod3_E_slabLate",&Hod3_E_slabLate,&b_Hod3_E_slabLate);
                t->SetBranchAddress("Hod2_E_slabLate",&Hod2_E_slabLate,&b_Hod2_E_slabLate);

                // setting other common branches
                    // Plastic scintillators
                t->SetBranchAddress("Hod3_MaxEprompt", &Hod3_MaxEMyprompt);
                t->SetBranchAddress("Hod3Hit_prompt", &Hod3Hit_Myprompt);

                t->SetBranchAddress("Hod2_MaxEprompt", &Hod2_MaxEMyprompt);
                t->SetBranchAddress("Hod2Hit_prompt", &Hod2Hit_Myprompt);

                //MC truth
                t->SetBranchAddress("MCEnergy", &MCEnergy);
                t->SetBranchAddress("HGasN", &HGasN);
                t->SetBranchAddress("Hod2N", &Hod2N);
                t->SetBranchAddress("Hod3N", &Hod3N);
                t->SetBranchAddress("EventID", &EventID);
                
                t->SetBranchAddress("beta", &beta);
            }



                // ------------------- READING THE TREE ------------------------ //
            Long64_t counting = 0;
            int NHod3Prompt = 0, NHod3Late = 0,
                NHod2Prompt = 0, NHod2Late = 0;

            for (Int_t entry = 0; entry < entries; entry++)
            {
                NHod3Prompt = 0; NHod3Late = 0;
                NHod2Prompt = 0; NHod2Late = 0;

                t->GetEntry(entry);                

                //counting the number of scintillators above the PROMPT threshold
                // and the number of scintillators between [high, low] thresholds
                for(int i=0; i<Hod3_E_slabMyPrompt->size(); i++) {
                    if(Hod3_E_slabMyPrompt->at(i)>thr_HodPrompt) NHod3Prompt++;
                    if(Hod2_E_slabMyPrompt->at(i)>thr_HodPrompt) NHod2Prompt++;

                    // if(thr_HodPrompt > Hod3_E_slabLate->at(i) && Hod3_E_slabLate->at(i)>thr_HodDelayed) NHod3Late++;
                    // if(thr_HodPrompt > Hod2_E_slabLate->at(i) && Hod2_E_slabLate->at(i)>thr_HodDelayed) NHod2Late++;

                    if(thr_HodPrompt > Hod3_E_slabMyPrompt->at(i) && Hod3_E_slabMyPrompt->at(i) > thr_HodDelayed) NHod3Late++;
                    if(thr_HodPrompt > Hod2_E_slabMyPrompt->at(i) && Hod2_E_slabMyPrompt->at(i) > thr_HodDelayed) NHod2Late++;
                }

                if( (Hod3_MaxEMyprompt > thr_HodDelayed) && 
                    (Hod2_MaxEMyprompt > thr_HodDelayed) && 
                    (HeCal_MaxEMyprompt > thr_HeCal) &&
                    (NHod3Late > thr_HodHitDelayed) &&
                    (NHod2Late > thr_HodHitDelayed) 
                    ){
                        // filling the prompt histograms
                        h2_HeVSMCEnergy[i]->Fill(MCEnergy, HeCal_MaxEMyprompt);

                        h2_Hod3VSMCEnergy[i]->Fill(MCEnergy, Hod3_MaxEMyprompt);
                        h2_Hod2VSMCEnergy[i]->Fill(MCEnergy, Hod2_MaxEMyprompt);

                        h2_Hod3NVSMCEnergy[i]->Fill(MCEnergy, NHod3Prompt);
                        h2_Hod2NVSMCEnergy[i]->Fill(MCEnergy, NHod2Prompt);
                    }
                
                //filling the delayed histograms
                h2_Hod3NDelVSMCEnergy[i]->Fill(MCEnergy, NHod3Late);
                h2_Hod2NDelVSMCEnergy[i]->Fill(MCEnergy, NHod2Late);

                h2_HeDelVSMCEnergy[i]->Fill(MCEnergy, HeCal_MaxELate);

                h2_Hod3DelVSMCEnergy[i]->Fill(MCEnergy, Hod3_MaxELate);
                h2_Hod2DelVSMCEnergy[i]->Fill(MCEnergy, Hod2_MaxELate);
            }
            printf("Number of fillments: %f\n", h2_HeVSMCEnergy[i]->GetEntries());
        }
        
    }

    // ------------------- PLOTTING ------------------- //
    const double TopMargin = 0.025, RightMargin = 0.02, LeftMargin = 0.;
    const double labSize = 0.065, TitleOffSetX = -0.5, TitleOffSetY=0.6;
    const int lineWidth = 4;
    double MaxContent[3];
    double MaxContent2[2];
    TH2D * h2_dummy[3];
    TH2D *h2_dummy2[2];
    TFile *file = new TFile(outFileName, "RECREATE");
    file->cd();

    TCanvas *c_EMax = new TCanvas("c_EMax", "c_EMax", 1100, 800);
    c_EMax->Divide(2,3,0.001, 0.001); 
    TCanvas *c_2 = new TCanvas("c_2", "c_2", 1100, 800);
    c_2->Divide(2,2,0.001, 0.001); 

    for(int j=0; j<nPlot; j++){
        switch (j){
            case 0:{
                c_EMax->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetRightMargin(RightMargin);

                //PROMPT HOD3
                h2_Hod3VSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod3VSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod3VSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod3VSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod3VSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod3VSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod3VSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} Hod 2");
                h2_Hod3VSMCEnergy[0]->Write();

                
                TH2D* h2_copyHod3VSMCEnergy = (TH2D*)h2_Hod3VSMCEnergy[0]->Clone();
                h2_copyHod3VSMCEnergy->SetName("Hod3_MaxEMyprompt_Comulative");
                h2_copyHod3VSMCEnergy->SetTitle("");

                
                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod3VSMCEnergy->Add(h2_Hod3VSMCEnergy[i]);

                    h2_Hod3VSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod3VSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod3VSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod3VSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod3VSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod3VSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod3VSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} Hod 2");

                    h2_Hod3VSMCEnergy[i]->Write();
                }
                MaxContent[0] = h2_copyHod3VSMCEnergy->GetBinContent(h2_copyHod3VSMCEnergy->GetMaximumBin());
                printf("\nMaxContent[0] = %f\n", MaxContent[0]);
                h2_copyHod3VSMCEnergy->DrawCopy("col");
                h2_copyHod3VSMCEnergy->Write();
                TLine *thr_HodPrompt_line = new TLine(10, thr_HodPrompt, 10000, thr_HodPrompt);
                thr_HodPrompt_line->SetLineColor(kRed);
                // thr_HodPrompt_line->SetLineStyle(5);
                thr_HodPrompt_line->SetLineWidth(lineWidth);
                thr_HodPrompt_line->Draw("same");
                

                c_2->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetRightMargin(RightMargin);
                //PROMPT HOD HIT
                h2_Hod2NVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod2NVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod2NVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod2NVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod2NVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod2NVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod2NVSMCEnergy[0]->GetYaxis()->SetTitle("Number of slab Hod 1");
                h2_Hod2NVSMCEnergy[0]->GetYaxis()->SetNdivisions(210);
                h2_Hod2NVSMCEnergy[0]->Write();

                TH2D* h2_copyHod2NVSMCEnergy = (TH2D*)h2_Hod2NVSMCEnergy[0]->Clone();
                h2_copyHod2NVSMCEnergy->GetXaxis()->SetLabelSize(labSize-0.01);
                h2_copyHod2NVSMCEnergy->GetXaxis()->SetTitleSize(labSize-0.01);
                h2_copyHod2NVSMCEnergy->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_copyHod2NVSMCEnergy->GetYaxis()->SetLabelSize(labSize-0.01);
                h2_copyHod2NVSMCEnergy->GetYaxis()->SetTitleSize(labSize-0.01);
                h2_copyHod2NVSMCEnergy->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_copyHod2NVSMCEnergy->GetYaxis()->SetTitle("Number of slab Hod 1");
                h2_copyHod2NVSMCEnergy->SetName("N_Hod2_slab_prompt_Comulative");
                h2_copyHod2NVSMCEnergy->SetTitle("");

                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod2NVSMCEnergy->Add(h2_Hod3NVSMCEnergy[i]);

                    h2_Hod2NVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod2NVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod2NVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod2NVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod2NVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod2NVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod2NVSMCEnergy[i]->GetYaxis()->SetTitle("Number of slab Hod 1");
                    h2_Hod2NVSMCEnergy[i]->GetYaxis()->SetNdivisions(210);
                    

                    h2_Hod2NVSMCEnergy[i]->Write();
                }
                MaxContent2[0] = h2_copyHod2NVSMCEnergy->GetBinContent(h2_copyHod2NVSMCEnergy->GetMaximumBin());
                h2_copyHod2NVSMCEnergy->DrawCopy("col");
                h2_copyHod2NVSMCEnergy->Write();
                TLine *thr_HodHitPrompt_line = new TLine(10, thr_HodHitPrompt, 10000, thr_HodHitPrompt);
                thr_HodHitPrompt_line->SetLineColor(kRed);
                // thr_HodHitPrompt_line->SetLineStyle(5);
                thr_HodHitPrompt_line->SetLineWidth(lineWidth);
                thr_HodHitPrompt_line->Draw("same");
                break;
                
            }

            case 1:{
                c_EMax->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz(); 
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetLeftMargin(LeftMargin);
                //DELAYED HOD3

                h2_Hod3DelVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod3DelVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod3DelVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod3DelVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod3DelVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod3DelVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod3DelVSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} Hod 2");
                h2_Hod3DelVSMCEnergy[0]->Write();

                TH2D* h2_copyHod3DelVSMCEnergy = (TH2D*)h2_Hod3DelVSMCEnergy[0]->Clone();
                h2_copyHod3DelVSMCEnergy->SetName("Hod3_MaxELate_Comulative");
                h2_copyHod3DelVSMCEnergy->SetTitle("");
                h2_copyHod3DelVSMCEnergy->SetMaximum(MaxContent[0]);

                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod3DelVSMCEnergy->Add(h2_Hod3DelVSMCEnergy[i]);

                    h2_Hod3DelVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod3DelVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod3DelVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod3DelVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod3DelVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod3DelVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod3DelVSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} Hod 2");

                    h2_Hod3DelVSMCEnergy[i]->Write();
                }
                h2_copyHod3DelVSMCEnergy->DrawCopy("colz");
                h2_dummy[0] = (TH2D*)gPad->GetPrimitive("Hod3_MaxELate_Comulative_copy");
                h2_dummy[0]->SetMaximum(MaxContent[0]);

                h2_copyHod3DelVSMCEnergy->Write();
                TLine *thr_HodDelayed_line = new TLine(10, thr_HodDelayed, 10000, thr_HodDelayed);
                thr_HodDelayed_line->SetLineColor(kRed);
                // thr_HodDelayed_line->SetLineStyle(5);
                thr_HodDelayed_line->SetLineWidth(lineWidth);
                thr_HodDelayed_line->Draw("same");
                
                c_2->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetLeftMargin(LeftMargin);
                //DELAYED HOD HIT
                h2_Hod2NDelVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod2NDelVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod2NDelVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod2NDelVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod2NDelVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod2NDelVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod2NDelVSMCEnergy[0]->GetYaxis()->SetTitle("Number of slab Hod 1");
                h2_Hod2NDelVSMCEnergy[0]->GetYaxis()->SetNdivisions(210);
                h2_Hod2NDelVSMCEnergy[0]->Write();

                TH2D* h2_copyHod2NDelVSMCEnergy = (TH2D*)h2_Hod2NDelVSMCEnergy[0]->Clone();
                h2_copyHod2NDelVSMCEnergy->GetXaxis()->SetLabelSize(labSize-0.01);
                h2_copyHod2NDelVSMCEnergy->GetXaxis()->SetTitleSize(labSize-0.01);
                h2_copyHod2NDelVSMCEnergy->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_copyHod2NDelVSMCEnergy->GetYaxis()->SetTitle("Number of slab Hod 1");
                h2_copyHod2NDelVSMCEnergy->GetYaxis()->SetLabelSize(labSize-0.01);
                h2_copyHod2NDelVSMCEnergy->GetYaxis()->SetTitleSize(labSize-0.01);
                h2_copyHod2NDelVSMCEnergy->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_copyHod2NDelVSMCEnergy->SetName("N_Hod2_slab_late_Comulative");
                h2_copyHod2NDelVSMCEnergy->SetTitle("");
                h2_copyHod2NDelVSMCEnergy->SetMaximum(MaxContent2[0]);
                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod2NDelVSMCEnergy->Add(h2_Hod2NDelVSMCEnergy[i]);

                    h2_Hod2NDelVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod2NDelVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod2NDelVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod2NDelVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod2NDelVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod2NDelVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod2NDelVSMCEnergy[i]->GetYaxis()->SetTitle("Number of slab Hod 1");
                    h2_Hod2NDelVSMCEnergy[i]->GetYaxis()->SetNdivisions(210);

                    h2_Hod2NDelVSMCEnergy[i]->Write();
                }
                h2_copyHod2NDelVSMCEnergy->DrawCopy("colz");
                h2_dummy2[0] = (TH2D*)gPad->GetPrimitive("N_Hod2_slab_late_Comulative_copy");
                h2_dummy2[0]->SetMaximum(MaxContent2[0]);
                h2_copyHod2NDelVSMCEnergy->Write();

                TLine *thr_HodHitDelayed_line = new TLine(10, thr_HodHitDelayed, 10000, thr_HodHitDelayed);
                thr_HodHitDelayed_line->SetLineColor(kRed);
                // thr_HodHitDelayed_line->SetLineStyle(5);
                thr_HodHitDelayed_line->SetLineWidth(lineWidth);
                thr_HodHitDelayed_line->Draw("same");
                break;
                
            }

            case 2:{
                c_EMax->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetRightMargin(RightMargin);
                //PROMPT HOD2
                h2_Hod2VSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod2VSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod2VSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod2VSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod2VSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod2VSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod2VSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} Hod 1");
                h2_Hod2VSMCEnergy[0]->Write();

                TH2D* h2_copyHod2VSMCEnergy = (TH2D*)h2_Hod2VSMCEnergy[0]->Clone();
                h2_copyHod2VSMCEnergy->SetName("Hod2_MaxEMyprompt_Comulative");
                h2_copyHod2VSMCEnergy->SetTitle("");
                
                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod2VSMCEnergy->Add(h2_Hod2VSMCEnergy[i]);

                    h2_Hod2VSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod2VSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod2VSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod2VSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod2VSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod2VSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod2VSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} Hod 1");

                    h2_Hod2VSMCEnergy[i]->Write();
                }
                MaxContent[1] = h2_copyHod2VSMCEnergy->GetBinContent(h2_copyHod2VSMCEnergy->GetMaximumBin());
                printf("\nMaxContent[1] = %f\n", MaxContent[1]);
                h2_copyHod2VSMCEnergy->DrawCopy("col");
                h2_copyHod2VSMCEnergy->Write();
                TLine *thr_HodPrompt_line = new TLine(10, thr_HodPrompt, 10000, thr_HodPrompt);
                thr_HodPrompt_line->SetLineColor(kRed);
                // thr_HodPrompt_line->SetLineStyle(5);
                thr_HodPrompt_line->SetLineWidth(lineWidth);
                thr_HodPrompt_line->Draw("same");
                
                c_2->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetRightMargin(RightMargin);
                //PROMPT HOD HIT
                h2_Hod3NVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod3NVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod3NVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod3NVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod3NVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod3NVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod3NVSMCEnergy[0]->GetYaxis()->SetTitle("Number of slab Hod 2");
                h2_Hod3NVSMCEnergy[0]->GetYaxis()->SetNdivisions(210);
                h2_Hod3NVSMCEnergy[0]->Write();

                TH2D* h2_copyHodNVSMCEnergy = (TH2D*)h2_Hod3NVSMCEnergy[0]->Clone();
                h2_copyHodNVSMCEnergy->GetXaxis()->SetLabelSize(labSize-0.01);
                h2_copyHodNVSMCEnergy->GetXaxis()->SetTitleSize(labSize-0.01);
                h2_copyHodNVSMCEnergy->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_copyHodNVSMCEnergy->GetYaxis()->SetTitle("Number of slab Hod 2");
                h2_copyHodNVSMCEnergy->GetYaxis()->SetLabelSize(labSize-0.01);
                h2_copyHodNVSMCEnergy->GetYaxis()->SetTitleSize(labSize-0.01);
                h2_copyHodNVSMCEnergy->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_copyHodNVSMCEnergy->SetName("N_Hod3_slab_prompt_Comulative");
                h2_copyHodNVSMCEnergy->SetTitle("");
                
                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHodNVSMCEnergy->Add(h2_Hod3NVSMCEnergy[i]);

                    h2_Hod3NVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod3NVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod3NVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod3NVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod3NVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod3NVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod3NVSMCEnergy[i]->GetYaxis()->SetTitle("Number of slab Hod 2");
                    h2_Hod3NVSMCEnergy[i]->GetYaxis()->SetNdivisions(210);

                    h2_Hod3NVSMCEnergy[i]->Write();
                }

                MaxContent2[1] = h2_copyHodNVSMCEnergy->GetBinContent(h2_copyHodNVSMCEnergy->GetMaximumBin());
                h2_copyHodNVSMCEnergy->DrawCopy("col");
                h2_copyHodNVSMCEnergy->Write();
                TLine *thr_HodHitPrompt_line = new TLine(10, thr_HodHitPrompt, 10000, thr_HodHitPrompt);
                thr_HodHitPrompt_line->SetLineColor(kRed);
                // thr_HodHitPrompt_line->SetLineStyle(5);
                thr_HodHitPrompt_line->SetLineWidth(lineWidth);
                thr_HodHitPrompt_line->Draw("same");
                break;
            }

            case 3:{
                c_EMax->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetLeftMargin(LeftMargin);
                //DELAYED HOD2
                h2_Hod2DelVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod2DelVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod2DelVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod2DelVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod2DelVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod2DelVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod2DelVSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} Hod 1");
                h2_Hod2DelVSMCEnergy[0]->Write();

                TH2D* h2_copyHod2DelVSMCEnergy = (TH2D*)h2_Hod2DelVSMCEnergy[0]->Clone();
                h2_copyHod2DelVSMCEnergy->SetName("Hod2_MaxELate_Comulative");
                h2_copyHod2DelVSMCEnergy->SetTitle("");
                h2_copyHod2DelVSMCEnergy->SetMaximum(MaxContent[1]);

                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHod2DelVSMCEnergy->Add(h2_Hod2DelVSMCEnergy[i]);

                    h2_Hod2DelVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod2DelVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod2DelVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod2DelVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod2DelVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod2DelVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod2DelVSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} Hod 1");

                    h2_Hod2DelVSMCEnergy[i]->Write();
                }
                h2_copyHod2DelVSMCEnergy->DrawCopy("colz");
                h2_copyHod2DelVSMCEnergy->Write();
                h2_dummy[1] = (TH2D*)gPad->GetPrimitive("Hod2_MaxELate_Comulative_copy");
                h2_dummy[1]->SetMaximum(MaxContent[1]);

                TLine *thr_HodDelayed_line = new TLine(10, thr_HodDelayed, 10000, thr_HodDelayed);
                thr_HodDelayed_line->SetLineColor(kRed);
                // thr_HodDelayed_line->SetLineStyle(5);
                thr_HodDelayed_line->SetLineWidth(lineWidth);
                thr_HodDelayed_line->Draw("same");

                c_EMax->cd(j+2);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetRightMargin(RightMargin);
                //PROMPT HECAL
                h2_HeVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_HeVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_HeVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_HeVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_HeVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_HeVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_HeVSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} HeCal");
                h2_HeVSMCEnergy[0]->Write();

                TH2D *h2_copyHeVSMCEnergy = (TH2D*)h2_HeVSMCEnergy[0]->Clone();
                h2_copyHeVSMCEnergy->SetName("HeCal_MaxEMyprompt_Comulative");
                h2_copyHeVSMCEnergy->SetTitle("");
                
                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHeVSMCEnergy->Add(h2_HeVSMCEnergy[i]);

                    h2_HeVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_HeVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_HeVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_HeVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_HeVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_HeVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_HeVSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} HeCal");

                    h2_HeVSMCEnergy[i]->Write();
                }

                MaxContent[2] = h2_copyHeVSMCEnergy->GetBinContent(h2_copyHeVSMCEnergy->GetMaximumBin());
                printf("\nMaxContent[2] = %f\n", MaxContent[2]);
                h2_copyHeVSMCEnergy->DrawCopy("col");
                h2_copyHeVSMCEnergy->Write();
                TLine *thr_line = new TLine(10, thr_HeCal, 10000, thr_HeCal);
                thr_line->SetLineColor(kRed);
                // thr_line->SetLineStyle(5);
                thr_line->SetLineWidth(lineWidth);
                thr_line->Draw("same");
                

                c_EMax->cd(j+3);
                gPad->SetLogx();
                gPad->SetLogy();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetLeftMargin(LeftMargin);
                //DELAYED HECAL
                h2_HeDelVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_HeDelVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_HeDelVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_HeDelVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_HeDelVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_HeDelVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_HeDelVSMCEnergy[0]->GetYaxis()->SetTitle("Max E_{dep} HeCal");
                h2_HeDelVSMCEnergy[0]->Write();

                TH2D* h2_copyHeDelVSMCEnergy = (TH2D*)h2_HeDelVSMCEnergy[0]->Clone();
                h2_copyHeDelVSMCEnergy->SetName("HeCal_MaxELate_Comulative");
                h2_copyHeDelVSMCEnergy->SetTitle("");
                h2_copyHeDelVSMCEnergy->SetMaximum(MaxContent[2]);

                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHeDelVSMCEnergy->Add(h2_HeDelVSMCEnergy[i]);

                    h2_HeDelVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_HeDelVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_HeDelVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_HeDelVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_HeDelVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_HeDelVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_HeDelVSMCEnergy[i]->GetYaxis()->SetTitle("Max E_{dep} HeCal");

                    h2_HeDelVSMCEnergy[i]->Write();
                }
                h2_copyHeDelVSMCEnergy->DrawCopy("colz");
                h2_dummy[2] = (TH2D*)gPad->GetPrimitive("HeCal_MaxELate_Comulative_copy");
                h2_dummy[2]->SetMaximum(MaxContent[2]);
                h2_copyHeDelVSMCEnergy->Write();
                TLine *thr_HeCal_line = new TLine(10, thr_HeCal, 10000, thr_HeCal);
                thr_HeCal_line->SetLineColor(kRed);
                // thr_HeCal_line->SetLineStyle(5);
                thr_HeCal_line->SetLineWidth(lineWidth);
                thr_HeCal_line->Draw("same");
                
                
                c_2->cd(j+1);
                gPad->SetLogx();
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetTopMargin(TopMargin);
                gPad->SetLeftMargin(LeftMargin);
                //DELAYED HOD HIT
                h2_Hod3NDelVSMCEnergy[0]->GetXaxis()->SetLabelSize(labSize);
                h2_Hod3NDelVSMCEnergy[0]->GetXaxis()->SetTitleSize(labSize);
                h2_Hod3NDelVSMCEnergy[0]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_Hod3NDelVSMCEnergy[0]->GetYaxis()->SetLabelSize(labSize);
                h2_Hod3NDelVSMCEnergy[0]->GetYaxis()->SetTitleSize(labSize);
                h2_Hod3NDelVSMCEnergy[0]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_Hod3NDelVSMCEnergy[0]->GetYaxis()->SetTitle("Number of slab Hod 2");
                h2_Hod3NDelVSMCEnergy[0]->GetYaxis()->SetNdivisions(210);
                h2_Hod3NDelVSMCEnergy[0]->Write();

                TH2D* h2_copyHodNDelVSMCEnergy = (TH2D*)h2_Hod3NDelVSMCEnergy[0]->Clone();
                h2_copyHodNDelVSMCEnergy->GetXaxis()->SetLabelSize(labSize-0.01);
                h2_copyHodNDelVSMCEnergy->GetXaxis()->SetTitleSize(labSize-0.01);
                h2_copyHodNDelVSMCEnergy->GetXaxis()->SetTitleOffset(TitleOffSetX);
                h2_copyHodNDelVSMCEnergy->GetYaxis()->SetTitle("Number of slab Hod 2");
                h2_copyHodNDelVSMCEnergy->GetYaxis()->SetLabelSize(labSize-0.01);
                h2_copyHodNDelVSMCEnergy->GetYaxis()->SetTitleSize(labSize-0.01);
                h2_copyHodNDelVSMCEnergy->GetYaxis()->SetTitleOffset(TitleOffSetY);
                h2_copyHodNDelVSMCEnergy->SetName("N_Hod3_slab_late_Comulative");
                h2_copyHodNDelVSMCEnergy->SetTitle("");
                h2_copyHodNDelVSMCEnergy->SetMaximum(MaxContent2[1]);

                for (int i = 1; i < nFiles; i++)
                {
                    h2_copyHodNDelVSMCEnergy->Add(h2_Hod3NDelVSMCEnergy[i]);

                    h2_Hod3NDelVSMCEnergy[i]->GetXaxis()->SetLabelSize(labSize);
                    h2_Hod3NDelVSMCEnergy[i]->GetXaxis()->SetTitleSize(labSize);
                    h2_Hod3NDelVSMCEnergy[i]->GetXaxis()->SetTitleOffset(TitleOffSetX);
                    h2_Hod3NDelVSMCEnergy[i]->GetYaxis()->SetLabelSize(labSize);
                    h2_Hod3NDelVSMCEnergy[i]->GetYaxis()->SetTitleSize(labSize);
                    h2_Hod3NDelVSMCEnergy[i]->GetYaxis()->SetTitleOffset(TitleOffSetY);
                    h2_Hod3NDelVSMCEnergy[i]->GetYaxis()->SetTitle("Number of slab Hod 2");
                    h2_Hod3NDelVSMCEnergy[i]->GetYaxis()->SetNdivisions(210);

                    h2_Hod3NDelVSMCEnergy[i]->Write();
                }
                h2_copyHodNDelVSMCEnergy->DrawCopy("colz");
                h2_dummy2[1] = (TH2D*)gPad->GetPrimitive("N_Hod3_slab_late_Comulative_copy");
                h2_dummy2[1]->SetMaximum(MaxContent2[1]);
                h2_copyHodNDelVSMCEnergy->Write();
                TLine *thr_HodHitDelayed_line = new TLine(10, thr_HodHitDelayed, 10000, thr_HodHitDelayed);
                thr_HodHitDelayed_line->SetLineColor(kRed);
                // thr_HodHitDelayed_line->SetLineStyle(5);
                thr_HodHitDelayed_line->SetLineWidth(lineWidth);
                thr_HodHitDelayed_line->Draw("same");
                break;
            }
        }
    }

    file->Write();
    file->Close();
    delete file;
}
