#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <map>
#include <TH2.h>
#include <TLine.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <TLatex.h>
#include <iostream>
using namespace std;

/*
Programma che richiede il valore del gate in ingresso (defualt e' 50 ns).
Legge in automatico il file di output di APP_meeting.cc e fornisce 
per ciascuna variabile il contributo di una singola particella. 
Variabli;      i.e. max E in HeCal      vs MCEnergy
                    Max E (hod3)        vs MCEnergy
                    Max E (hod2)        vs MCEnergy
                    # prompt hit (Hod3) vs MCEnergy.
                    # prompt hit (Hod2) vs MCEnergy.
                    # late hit (Hod3)   vs MCEnergy.
                    # late hit (Hod2)   vs MCEnergy
*/



//!!!!!!!!!!!!!!
// MAIN PROGRAM
//!!!!!!!!!!!!!!

void single_contribution(int gate)
{
    gStyle->SetOptStat(0); // no statistics
    TString fName;

    // ----------------- INPUT FILE ----------------- //
    switch (gate){
            case 50:{
                fName = "APP_meeting_50.root";
                printf("You selected a 50  ns gate\n");
                break;
            }
            case 60:{
                fName = "APP_meeting_60.root";
                printf("You selected a 60  ns gate\n");
                break;
            }
            // in case the variable gate assumes another value, rise an error flag
            default:{
                printf("ERROR: gate value not valid\n Currently available values are 50 and 60 ns\n");
                return;
            }
    }
    TFile *f = new TFile(fName, "OPEN");


    const int nParticles = 5;
    TString pNames[nParticles] = {"#bar{d}", "#bar{p}", "d", "p", "He4"};

    const int nVariables = 10;
    TString vNames[nVariables] = {"Hod3_Max_MaxEMyprompt_vs_", "Hod3_Max_MaxELate_vs_",
                                  "Hod2_Max_MaxEMyprompt_vs_", "Hod2_Max_MaxELate_vs_", 
                                  "HeCal_Max_MaxEMyprompt_vs_","HeCal_Max_MaxELate_vs_",
                                  "N_Hod3_slab_prompt_vs_", "N_Hod3_slab_late_vs_",
                                  "N_Hod2_slab_prompt_vs_", "N_Hod2_slab_late_vs_"};
    //histograms
    TH2D *h2_Dbar[nVariables];  
    TH2D *h2_Pbar[nVariables]; 
    TH2D *h2_P[nVariables];
    TH2D *h2_D[nVariables];
    TH2D *h2_He4[nVariables];

    
    
    TString hNames[nVariables*nParticles];
    for(int i=0; i<nVariables; i++){
        for(int j=0; j<nParticles; j++){
            hNames[i*nParticles+j] = vNames[i] + pNames[j] + "_energy";
            if(pNames[j].Contains("#bar{d}")) h2_Dbar[i] = (TH2D*)f->Get(hNames[i*nParticles+j]);
            if(pNames[j].Contains("#bar{p}")) h2_Pbar[i] = (TH2D*)f->Get(hNames[i*nParticles+j]);
            if(pNames[j].Contains("p")) h2_P[i] = (TH2D*)f->Get(hNames[i*nParticles+j]);
            if(pNames[j].Contains("d")) h2_D[i] = (TH2D*)f->Get(hNames[i*nParticles+j]);
            if(pNames[j].Contains("He4")) h2_He4[i] = (TH2D*)f->Get(hNames[i*nParticles+j]);
        }
    }

    //canvas names
    TCanvas *c[nVariables];
    TString cNames[nVariables] = {"Hod3_prompt", "Hod3_late", "Hod2_prompt", "Hod2_late",
                                 "HeCal_prompt", "HeCal_late", "NHod3_prompt", "NHod3_late",
                                 "NHod2_prompt", "NHod2_late"}; 

    // Thresholds
    const double thr_HeCalPrompt = 10.; // Prompt HeCal [MeV]
    const double thr_HeCalAfter = 10.;  // Delayed HeCal [MeV]
    const double thr_HodPrompt = 1.6;   // Prompt Scintillators [MeV]
    const double thr_HodDelayed = 0.8;  // Delayed Scintillators [MeV]
    const double thr_HodHitPrompt = 2.; // Prompt Hits in scintillators
    const double thr_HodHitDelayed = 2.;// Delayed Hits in scintillators

    const int lineWidth = 4;
    const int lineStyle = 2; // kDashed

    TLine *l_thresholds[nVariables];
    for(int i=0; i<nVariables; i++){
        switch(i){
            case 0: case 2:
                l_thresholds[i] = new TLine(10, thr_HodPrompt, 10000, thr_HodPrompt);
                break;
            case 1: case 3:
                l_thresholds[i] = new TLine(10, thr_HodDelayed, 10000, thr_HodDelayed);
                break;
            case 4:
                l_thresholds[i] = new TLine(10, thr_HeCalPrompt, 10000, thr_HeCalPrompt);
                break;
            case 5:
                l_thresholds[i] = new TLine(10, thr_HeCalAfter, 10000, thr_HeCalAfter);
                break;
            case 6: case 8:
                l_thresholds[i] = new TLine(10, thr_HodHitPrompt, 10000, thr_HodHitPrompt);
                break;
            case 7: case 9:
                l_thresholds[i] = new TLine(10, thr_HodHitDelayed, 10000, thr_HodHitDelayed);
                break;
        }
        l_thresholds[i]->SetLineColor(kRed);
        l_thresholds[i]->SetLineWidth(lineWidth);
        l_thresholds[i]->SetLineStyle(lineStyle);
    }

    

    printf("Thresholds:\n");
    printf("\t HeCal prompt > %.2f MeV\n", thr_HeCalPrompt);
    printf("\t HeCal delayed > %.2f MeV\n", thr_HeCalAfter);
    printf("\t Scintillators prompt > %.2f MeV\n", thr_HodPrompt);
    printf("\t Scintillators delayed > %.2f MeV\n", thr_HodDelayed);
    printf("\t Hits prompt < %.2f\n", thr_HodHitPrompt);
    printf("\t Hits delayed > %.2f\n", thr_HodHitDelayed);



    // ------------------- PLOTTING ------------------- //
    // TFile *file = new TFile("APP_meeting.root", "RECREATE");
    // file->cd();
    const double TopMargin = 0.025, RightMargin = 0.02, LeftMargin = 0.;
    const double MaxContent = 2000;


    for(int i=0; i<nVariables; i++){
        c[i] = new TCanvas(cNames[i], cNames[i], 1100, 800);
        c[i]->Divide(2,3,0.001, 0.001);

        c[i]->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetTopMargin(TopMargin);
        gPad->SetRightMargin(RightMargin);


        if(i<6) gPad->SetLogy();
        gPad->SetLogz();
        h2_Dbar[i]->SetTitle("");
        h2_Dbar[i]->SetMaximum(MaxContent);
        h2_Dbar[i]->DrawCopy("col");
        l_thresholds[i]->Draw("same");


        c[i]->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetTopMargin(TopMargin);
        gPad->SetLeftMargin(LeftMargin);

        if(i<6) gPad->SetLogy();
        gPad->SetLogz();
        h2_Pbar[i]->SetTitle("");
        h2_Pbar[i]->SetMaximum(MaxContent);
        h2_Pbar[i]->DrawCopy("colz");
        l_thresholds[i]->Draw("same");

        c[i]->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetTopMargin(TopMargin);
        gPad->SetRightMargin(RightMargin);

        if(i<6) gPad->SetLogy();
        gPad->SetLogz();
        h2_D[i]->SetTitle("");
        h2_D[i]->SetMaximum(MaxContent);
        h2_D[i]->DrawCopy("col");
        l_thresholds[i]->Draw("same");


        c[i]->cd(4);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetTopMargin(TopMargin);
        gPad->SetLeftMargin(LeftMargin);

        if(i<6) gPad->SetLogy();
        h2_P[i]->SetTitle("");
        h2_P[i]->SetMaximum(MaxContent);
        h2_P[i]->DrawCopy("colz");
        l_thresholds[i]->Draw("same");
        gPad->SetLogz();

        c[i]->cd(5);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetTopMargin(TopMargin);
        gPad->SetRightMargin(RightMargin);

        if(i<6) gPad->SetLogy();
        gPad->SetLogz();
        h2_He4[i]->SetTitle("");
        h2_He4[i]->SetMaximum(MaxContent);
        h2_He4[i]->DrawCopy("col");
        l_thresholds[i]->Draw("same");
        
    }

    

    

    // file->Write();
    // file->Close();
    // delete file;
}
