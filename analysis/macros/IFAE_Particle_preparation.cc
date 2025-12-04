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
#include <TGraph.h>
#include <TGraphErrors.h>
using namespace std;





void IFAE_Particle_preparation(){

    // gStyle->SetOptStat(0); // no statistics

    const int NFiles = 2;
    const int NParticles = 2;
    TString particles[NParticles] = {"D", "P"};
    //, "Deuteron", "proton", "^{4}He", "^{12}C", "e^{-}"};
    TString particles_name[NParticles] = {"D", "P"};
    //, "Deuteron", "Proton", "He4", "C12", "electron"};

    // "Deuteron_10e7_gps_rec__50ns.root" , "proton_10e7_gps_rec__50ns.root"

    TString names[NFiles] = {"../../sources/ROOT_files/Deuteron_10e7_gps_SMrec__50ns.root", "../../sources/ROOT_files/proton_10e7_gps_SMrec__50ns.root"};

                            //  "../../sources/ROOT_files/Deuteron_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/proton_10e7_gps_rec__50ns.root", 
                            //  "../../sources/ROOT_files/He4_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/C12_10e7_gps_rec__50ns.root", 
                            //  "../../sources/ROOT_files/electron_10e7_gps_rec__50ns.root"}; 

    enum particle_order {D, P};
    // , Deuteron, Proton, He4, C12, Electron};

    //weight for antiP different statistics
    const double weight = 3;

    //HeCal binning prompt
    const int nBins_HeCalP = 90;
    const double HeCal_minP = 10, HeCal_maxP = 80; // MeV

    //HeCal binning prompt
    const int nBins_HeCalL = 250;
    const double HeCal_minL = 10, HeCal_maxL = 160; // MeV

    //Beta binning
    const int nBins_Beta = 150;
    const double Beta_min = 0.2, Beta_max = 0.5;

    //MIPs scaled binning
    const int nBins_MIP = 100;
    const double MIP_min = 1.0, MIP_max = 10.0;

    //Hod binning
    const int nBins_Hod = 260;
    const double Hod_min = 1.0, Hod_max = 14.0;

    //Hod slab binning
    const int nBins_HodSlab = 10;
    const double HodSlab_min = 0.5, HodSlab_max = 10.5;

    //HeCal tank binning
    const int nBins_HeCalPTank = 10;
    const double HeCalTank_min = 0.5, HeCalTank_max = 10.5;
    
                            
    TH2D *h2_Eprompt_VS_Beta[NParticles];
    TH2D *h2_TOTAL_Eprompt_VS_Beta = new TH2D("h2_TOTAL_Eprompt_VS_Beta", "HeCal prompt energy vs TOF #beta; HeCal Prompt Energy [MeV]; TOF #beta", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_Beta, Beta_min, Beta_max);

        h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_Eprompt_VS_Beta->GetXaxis()->SetTitleOffset(1.2);

    TH2D *h2_EHecal_VS_Hod2[NParticles];
    TH2D *h2_EHecal_VS_Hod3[NParticles];
    TH2D * h2_TOTAL_EHecal_VS_Hod2 = new TH2D("h2_TOTAL_EHecal_VS_Hod2", "HeCal prompt energy vs HOD2 prompt energy; HeCal Prompt Energy [MeV]; Hod2 Prompt Energy [MeV]", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_Hod, Hod_min, Hod_max);
    TH2D * h2_TOTAL_EHecal_VS_Hod3 = new TH2D("h2_TOTAL_EHecal_VS_Hod3", "HeCal prompt energy vs HOD3 prompt energy; HeCal Prompt Energy [MeV]; Hod3 Prompt Energy [MeV]", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_MIP, MIP_min, MIP_max);
    TH2D * h2_TOTAL_EHecal_VS_Hod = new TH2D("h2_TOTAL_EHecal_VS_Hod", "HeCal prompt energy vs HOD prompt energy; HeCal Prompt Energy [MeV]; Hod Prompt Energy [MeV]", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_Hod*2, Hod_min*2, Hod_max*2);

        h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetTitleOffset(1.0);

        h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetNdivisions(210);

        h2_TOTAL_EHecal_VS_Hod2->GetXaxis()->SetTitleOffset(1.2);
        h2_TOTAL_EHecal_VS_Hod3->GetXaxis()->SetTitleOffset(1.2);
        h2_TOTAL_EHecal_VS_Hod->GetXaxis()->SetTitleOffset(1.2);

    TH2D *h2_EHecalL_VS_Hod2SlabL[NParticles];
    TH2D *h2_EHecalL_VS_Hod3SlabL[NParticles];
    TH2D* h2_EHecalL_VS_HodSlabL[NParticles];
    TH2D * h2_TOTAL_EHecalL_VS_Hod2SlabL = new TH2D("h2_TOTAL_EHecalL_VS_Hod2SlabL", "HeCal late energy vs HOD2 slab number; HeCal Late Energy [MeV]; Hod2 Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab, HodSlab_min, HodSlab_max);
    TH2D * h2_TOTAL_EHecalL_VS_Hod3SlabL = new TH2D("h2_TOTAL_EHecalL_VS_Hod3SlabL", "HeCal late energy vs HOD3 slab number; HeCal Late Energy [MeV]; Hod3 Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab, HodSlab_min, HodSlab_max);
    TH2D * h2_TOTAL_EHecalL_VS_HodSlabL = new TH2D("h2_TOTAL_EHecalL_VS_HodSlabL", "HeCal late energy vs HOD slab number; HeCal Late Energy [MeV]; Hod Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab*2+1, HodSlab_min*2, HodSlab_max*2+1);

        h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetTitleOffset(1.0);

        h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetNdivisions(210);

        h2_TOTAL_EHecalL_VS_Hod2SlabL->GetXaxis()->SetTitleOffset(1.2);
        h2_TOTAL_EHecalL_VS_Hod3SlabL->GetXaxis()->SetTitleOffset(1.2);
        h2_TOTAL_EHecalL_VS_HodSlabL->GetXaxis()->SetTitleOffset(1.2);

    TH2D *h2_EHecalL_VS_HeCalTankL[NParticles];
    TH2D * h2_TOTAL_EHecalL_VS_HeCalTankL = new TH2D("h2_TOTAL_EHecalL_VS_HeCalTankL", "HeCal late energy vs HeCal tank number; HeCal Late Energy [MeV]; HeCal Tank Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HeCalPTank, HeCalTank_min, HeCalTank_max);

        h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecalL_VS_HeCalTankL->GetXaxis()->SetTitleOffset(1.2);

    TH2D *h2_TOTAL_EHecalL_VS_TotalHit_part[NParticles];

    TH2D *h2_TOTAL_EHecalL_VS_TotalHit = new TH2D("h2_TOTAL_EHecalL_VS_TotalHit", "HeCal late energy vs Total Hit Number; HeCal Late Energy [MeV]; Total Hit Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab*2+nBins_HeCalPTank, HodSlab_min*2+HeCalTank_min, HodSlab_max*2+HeCalTank_max);
    
            h2_TOTAL_EHecalL_VS_TotalHit->GetYaxis()->SetTitleOffset(1.0);
            h2_TOTAL_EHecalL_VS_TotalHit->GetYaxis()->SetNdivisions(210);
            h2_TOTAL_EHecalL_VS_TotalHit->GetXaxis()->SetTitleOffset(1.2);  

    const int LineWidth = 3;
    Color_t colors[NParticles] = {kRed, kBlue};
    //, kGreen+1, kMagenta, kOrange, kBlack, kCyan-3} ;


    for(int i=0; i<NParticles; i++) {
        h2_Eprompt_VS_Beta[i] = new TH2D("h2_Eprompt_VS_Beta_"+particles_name[i], "HeCal prompt energy vs TOF #beta_"+particles[i]+"; HeCal Prompt Energy [MeV]; TOF #beta_"+particles[i], nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_Beta, Beta_min, Beta_max);
        h2_EHecal_VS_Hod2[i] = new TH2D("h2_EHecal_VS_Hod2_"+particles_name[i], "HeCal prompt energy vs HOD2 prompt energy "+particles[i]+"; HeCal Prompt Energy [MeV]; Hod2 Prompt Energy [MeV]", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_Hod, Hod_min, Hod_max);
        h2_EHecal_VS_Hod3[i] = new TH2D("h2_EHecal_VS_Hod3_"+particles_name[i], "HeCal prompt energy vs HOD3 prompt energy "+particles[i]+"; HeCal Prompt Energy [MeV]; dE/dx [MIPs]", nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_MIP, MIP_min, MIP_max);

        // h2_Eprompt_VS_Beta[i] = GenerateHistogram("h2_Eprompt_VS_Beta_"+particles_name[i], "HeCal prompt energy vs TOF #beta_"+particles[i], "HeCal Prompt Energy [MeV]", "TOF #beta_"+particles[i], nBins_HeCalP, HeCal_minP, HeCal_maxP, nBins_HeCalP, 0.1, 1., false);
        // h2_EHecal_VS_Hod2[i] = GenerateHistogram("h2_EHecal_VS_Hod2_"+particles_name[i], "HeCal prompt energy vs HOD prompt energy "+particles[i], "HeCal Prompt Energy [MeV}", "Hod Prompt Energy [MeV]", nBins_HeCalP, HeCal_minP, HeCal_maxP);

        h2_EHecalL_VS_Hod2SlabL[i] = new TH2D("h2_EHecalL_VS_Hod2SlabL_"+particles_name[i], "HeCal late energy vs HOD2 slab number "+particles[i]+"; HeCal Late Energy [MeV]; Hod2 Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab, HodSlab_min, HodSlab_max);
        h2_EHecalL_VS_Hod3SlabL[i] = new TH2D("h2_EHecalL_VS_Hod3SlabL_"+particles_name[i], "HeCal late energy vs HOD3 slab number "+particles[i]+"; HeCal Late Energy [MeV]; Hod3 Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab, HodSlab_min, HodSlab_max);
        h2_EHecalL_VS_HodSlabL[i] = new TH2D("h2_EHecalL_VS_HodSlabL_"+particles_name[i], "HeCal late energy vs HOD slab number "+particles[i]+"; HeCal Late Energy [MeV]; Total Hod Slab Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab*2+1, HodSlab_min*2, HodSlab_max*2+1);

        h2_EHecalL_VS_HeCalTankL[i] = new TH2D("h2_EHecalL_VS_HeCalTankL_"+particles_name[i], "HeCal late energy vs HeCal tank number "+particles[i]+"; HeCal Late Energy [MeV]; HeCal Tank Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HeCalPTank, HeCalTank_min, HeCalTank_max);

        h2_TOTAL_EHecalL_VS_TotalHit_part[i] = new TH2D("h2_TOTAL_EHecalL_VS_TotalHit_"+particles_name[i], "HeCal late energy vs Total Hit Number; HeCal Late Energy [MeV]; Total Hit Number", nBins_HeCalL, HeCal_minL, HeCal_maxL, nBins_HodSlab*2+nBins_HeCalPTank, HodSlab_min*2+HeCalTank_min, HodSlab_max*2+HeCalTank_max);
                
        h2_Eprompt_VS_Beta[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_Eprompt_VS_Beta[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_EHecal_VS_Hod2[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_EHecal_VS_Hod2[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetNdivisions(210);
        h2_EHecalL_VS_Hod2SlabL[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetNdivisions(210);
        h2_EHecalL_VS_Hod3SlabL[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetNdivisions(210);
        h2_EHecalL_VS_HodSlabL[i]->GetXaxis()->SetTitleOffset(1.2);

        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetNdivisions(210);
        h2_EHecalL_VS_HeCalTankL[i]->GetXaxis()->SetTitleOffset(1.2);  

        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetTitleOffset(1.0);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetNdivisions(210);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetXaxis()->SetTitleOffset(1.2);
    }

    const double th_Hod_high = 1.6;
    const double th_Hod_low = 0.8;
    const double th_HeCal = 10.0;
    const double th_HodHit = 2.0;
    const double th_HodHitTot = 3.0;

    double MCEnergy;
    //prompt variables
    int NHod3Prompt, NHod2Prompt;
    double Hod3_MaxEMyprompt, Hod2_MaxEMyprompt, HeCal_MaxEMyprompt; 
    double TOTHeCalPrompt=0, TOTHod2Prompt=0, TOTHod3Prompt=0;
    double beta=0;

    vector<double> *Hod2X=0, *Hod2Y=0, *Hod2Z=0;
    vector<double> *Hod3X=0, *Hod3Y=0, *Hod3Z=0;
    TBranch *b_Hod2X=0, *b_Hod2Y=0, *b_Hod2Z=0;
    TBranch *b_Hod3X=0, *b_Hod3Y=0, *b_Hod3Z=0;
    const double Hod3Thick = 4.0 / 10.; // cm
    const double MIP_norm = 0.8 / Hod3Thick; // MeV/cm
    double PathLenght = 0;


    //delayed variables
    int NHod3Late, NHod2Late;
    double Hod3_MaxELate, Hod2_MaxELate, HeCal_MaxELate;
    double TOTHeCalLate=0, TOTHod2Slablate=0, TOTHod3SlabLate=0, TOTHeCalTankLate=0;
        //clusterizzazione HeCal
    vector<double> *HeCal_E_tankLate=0;
    TBranch *b_HeCal_E_tankLate=0;
        //clusterizzazione degli scintillatori: HOD3 e HOD2
    vector<double> *Hod3_E_slabLate=0, *Hod2_E_slabLate=0;
    TBranch *b_Hod3_E_slabLate=0, *b_Hod2_E_slabLate=0;

    const int NCut = 5;
    TString prompt_Sel[NCut] = {"Hod3_MaxEMyprompt", "Hod2_MaxEMyprompt", "HeCal_MaxEMyprompt", "NHod3Prompt", "NHod2Prompt"};
    TString prompt_Del[NCut] = {"Hod3_MaxELate", "Hod2_MaxELate", "HeCal_MaxELate", "NHod3Late", "NHod2Late"};

    printf("\n ------- CURRENT PROMPT SELECTION -------\n");
    for(int j=0;j<NCut;j++){
        switch(j){
            case 0: case 1:
                printf(" %s \t > %.2f \n", prompt_Sel[j].Data(), th_Hod_high);
                break;
            case 2:
                printf(" %s \t > %.2f \n", prompt_Sel[j].Data(), th_HeCal);
                break;
            case 3: case 4:
                printf(" %s \t < %.2f \n", prompt_Sel[j].Data(), th_HodHit);
                printf(" %s \t > %.2f \n", prompt_Sel[j].Data(), 0.0);
                if(j==4) {
                    printf(" %s + %s \t < %.2f \n", prompt_Sel[3].Data(), prompt_Sel[4].Data(), th_HodHitTot);
                }
                break;
        }
    }

    printf("\n ------- CURRENT DELAYED SELECTION -------\n");
    for(int j=0;j<NCut;j++){
        switch(j){
            case 0: case 1:
                printf(" %s \t > %.2f \n", prompt_Del[j].Data(), th_Hod_high);
                break;
            case 2:
                printf(" %s \t > %.2f \n", prompt_Del[j].Data(), th_HeCal);
                break;
            case 3: case 4:
                printf(" %s \t > %.2f \n", prompt_Del[j].Data(), th_HodHit);
                // if(j==4) {
                //     printf(" %s + %s \t < %.2f \n", prompt_Del[3].Data(), prompt_Del[4].Data(), th_HodHitTot);
                // }
                break;
        }     
    }


    for(int i=0; i<NFiles; i++){

        TFile *file = new TFile(names[i], "OPEN");
        printf("\r Reading file %d/%d: %s", i+1, NFiles, names[i].Data());
        fflush(stdout);

        TTree *t;
        if(names[i].Contains("bar")) { //antiparticles
            if(names[i].Contains("__")) t = (TTree *)file->Get("Hits_bar");
            else t = (TTree *)file->Get("RHits_bar");
        }
        else t = (TTree *)file->Get("RHits");       

        t->SetBranchAddress("MCEnergy", &MCEnergy);

        //prompt variables
        t->SetBranchAddress("Hod3_MaxEMyprompt", &Hod3_MaxEMyprompt);
        t->SetBranchAddress("Hod2_MaxEMyprompt", &Hod2_MaxEMyprompt);
        t->SetBranchAddress("HeCal_MaxEMyprompt", &HeCal_MaxEMyprompt);
        t->SetBranchAddress("beta", &beta);

        //scintillators position
        t->SetBranchAddress("Hod3X",&Hod3X,&b_Hod3X);
        t->SetBranchAddress("Hod3Y",&Hod3Y,&b_Hod3Y);
        t->SetBranchAddress("Hod3Z",&Hod3Z,&b_Hod3Z);
        t->SetBranchAddress("Hod2X",&Hod2X,&b_Hod2X);
        t->SetBranchAddress("Hod2Y",&Hod2Y,&b_Hod2Y);
        t->SetBranchAddress("Hod2Z",&Hod2Z,&b_Hod2Z);

        //clusterizzazione HeCal
        vector<double> *HeCal_E_tankMyPrompt=0;
        TBranch *b_HeCal_E_tankMyPrompt=0;
        t->SetBranchAddress("HeCal_E_tankMyPrompt",&HeCal_E_tankMyPrompt,&b_HeCal_E_tankMyPrompt);

        //clusterizzazione degli scintillatori: HOD3 e HOD2
        vector<double> *Hod3_E_slabMyPrompt=0, *Hod2_E_slabMyPrompt=0;
        TBranch *b_Hod3_E_slabMyPrompt=0, *b_Hod2_E_slabMyPrompt=0;
        t->SetBranchAddress("Hod3_E_slabMyPrompt",&Hod3_E_slabMyPrompt,&b_Hod3_E_slabMyPrompt);
        t->SetBranchAddress("Hod2_E_slabMyPrompt",&Hod2_E_slabMyPrompt,&b_Hod2_E_slabMyPrompt);


        //delayed variables
        if(names[i].Contains("bar")) { //antiparticles
            t->SetBranchAddress("Hod3_MaxELate", &Hod3_MaxELate);
            t->SetBranchAddress("Hod2_MaxELate", &Hod2_MaxELate);
            t->SetBranchAddress("HeCal_MaxELate", &HeCal_MaxELate);

            //clusterizzazione HeCal
            t->SetBranchAddress("HeCal_E_tankLate",&HeCal_E_tankLate,&b_HeCal_E_tankLate);

            //clusterizzazione degli scintillatori: HOD3 e HOD2
            t->SetBranchAddress("Hod3_E_slabLate",&Hod3_E_slabLate,&b_Hod3_E_slabLate);
            t->SetBranchAddress("Hod2_E_slabLate",&Hod2_E_slabLate,&b_Hod2_E_slabLate);
        }
        


        for(int entry=0; entry<t->GetEntries(); entry++){
            t->GetEntry(entry);

            //Path length in Hod 3
            PathLenght = 0;

            //Number of slab Prompt
            NHod3Prompt = 0; NHod2Prompt = 0;
            //Number of slab Late
            NHod3Late = 0; NHod2Late = 0;
            //Number of tank Late
            TOTHeCalTankLate = 0;
            
            //Total Prompt Energy
            TOTHeCalPrompt = 0;
            TOTHod2Prompt = 0;
            TOTHod3Prompt = 0;

            //Total Late Energy
            TOTHeCalLate = 0;


            //counting the number of scintillators above the PROMPT threshold
            for(int iscint=0; iscint<Hod3_E_slabMyPrompt->size(); iscint++) {
                if(Hod3_E_slabMyPrompt->at(iscint)>th_Hod_high) NHod3Prompt++;
                if(Hod2_E_slabMyPrompt->at(iscint)>th_Hod_high) NHod2Prompt++;
            }
            


            // PROMPT SELECTION
            if( (Hod3_MaxEMyprompt > th_Hod_high) && 
                (Hod2_MaxEMyprompt > th_Hod_high) && 
                (HeCal_MaxEMyprompt > th_HeCal) && 
                (NHod3Prompt < th_HodHit) && 
                (NHod2Prompt < th_HodHit) &&  
                (NHod3Prompt > 0) && 
                (NHod2Prompt > 0) &&
                (NHod3Prompt+NHod2Prompt) < th_HodHitTot ){

                    for(int tank=0; tank<HeCal_E_tankMyPrompt->size(); tank++) {
                        if(HeCal_E_tankMyPrompt->at(tank)>0.) TOTHeCalPrompt += HeCal_E_tankMyPrompt->at(tank);
                    }
                    for(int slab=0; slab<Hod3_E_slabMyPrompt->size(); slab++) {
                        if(Hod3_E_slabMyPrompt->at(slab)>0.) TOTHod3Prompt += Hod3_E_slabMyPrompt->at(slab);
                        if(Hod2_E_slabMyPrompt->at(slab)>0.) TOTHod2Prompt += Hod2_E_slabMyPrompt->at(slab);
                    }


                    if(Hod3X->size()>0 && Hod2X->size()>0){
                        double dist = sqrt(pow((Hod2Y->at(0)-Hod3Y->at(0)),2) + pow((Hod2X->at(0)-Hod3X->at(0)),2) + pow((Hod2Z->at(0)-Hod3Z->at(0)),2));
                        double deltaZ = abs(Hod2Z->at(0)-Hod3Z->at(0));
                        double cosTheta = deltaZ/dist; 
                        PathLenght = Hod3Thick/(cosTheta); //cm
                    }

                    if(i%2 == 0){
                        h2_Eprompt_VS_Beta[D]->Fill(TOTHeCalPrompt, beta);
                        h2_EHecal_VS_Hod2[D]->Fill(TOTHeCalPrompt, TOTHod2Prompt);
                        h2_EHecal_VS_Hod3[D]->Fill(TOTHeCalPrompt, TOTHod3Prompt/(PathLenght*MIP_norm));

                        // h2_EHecalL_VS_Hod2SlabL[D]->Fill(TOTHeCalLate, NHod2Late);
                        // h2_EHecalL_VS_Hod3SlabL[D]->Fill(TOTHeCalLate, NHod3Late);
                        // h2_EHecalL_VS_HodSlabL[D]->Fill(TOTHeCalLate, NHod2Late+NHod3Late);
                        // h2_EHecalL_VS_HeCalTankL[D]->Fill(TOTHeCalLate, TOTHeCalTankLate);

                        // h2_TOTAL_EHecalL_VS_TotalHit_part[D]->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate);

                        h2_TOTAL_EHecal_VS_Hod->Fill(TOTHeCalPrompt, TOTHod2Prompt+TOTHod3Prompt);
                        // h2_TOTAL_EHecalL_VS_HodSlabL->Fill(TOTHeCalLate, NHod2Late+NHod3Late);
                        // h2_TOTAL_EHecalL_VS_TotalHit->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate);
                    }
                    else{
                        h2_Eprompt_VS_Beta[P]->Fill(TOTHeCalPrompt, beta, weight);
                        h2_EHecal_VS_Hod2[P]->Fill(TOTHeCalPrompt, TOTHod2Prompt, weight);
                        h2_EHecal_VS_Hod3[P]->Fill(TOTHeCalPrompt, TOTHod3Prompt/(PathLenght*MIP_norm), weight);

                        // h2_EHecalL_VS_Hod2SlabL[P]->Fill(TOTHeCalLate, NHod2Late, weight);
                        // h2_EHecalL_VS_Hod3SlabL[P]->Fill(TOTHeCalLate, NHod3Late, weight);
                        // h2_EHecalL_VS_HodSlabL[P]->Fill(TOTHeCalLate, NHod2Late+NHod3Late, weight);
                        // h2_EHecalL_VS_HeCalTankL[P]->Fill(TOTHeCalLate, TOTHeCalTankLate, weight);

                        // h2_TOTAL_EHecalL_VS_TotalHit_part[P]->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate, weight);

                        h2_TOTAL_EHecal_VS_Hod->Fill(TOTHeCalPrompt, TOTHod2Prompt+TOTHod3Prompt, weight);
                        // h2_TOTAL_EHecalL_VS_HodSlabL->Fill(TOTHeCalLate, NHod2Late+NHod3Late, weight);
                        // h2_TOTAL_EHecalL_VS_TotalHit->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate, weight);
                    }
                } // prompt selection


                NHod3Prompt = 0; NHod2Prompt = 0;

                //counting the number of scintillators above the DELAYED threshold
                for(int iscint=0; iscint<Hod3_E_slabMyPrompt->size(); iscint++){
                    if(th_Hod_high > Hod3_E_slabMyPrompt->at(iscint) && Hod3_E_slabMyPrompt->at(iscint)>th_Hod_low) NHod3Prompt++;
                    if(th_Hod_high > Hod3_E_slabMyPrompt->at(iscint) && Hod3_E_slabMyPrompt->at(iscint)>th_Hod_low) NHod2Prompt++;
                }

                // // DELAYED SELECTION
                // if( (Hod3_MaxEMyprompt > th_Hod_low) && 
                //     (Hod2_MaxEMyprompt > th_Hod_low) && 
                //     (HeCal_MaxEMyprompt > th_HeCal) && 
                //     (NHod3Prompt > th_HodHit) && 
                //     (NHod2Prompt > th_HodHit) ){


                //     for(int tank=0; tank<HeCal_E_tankMyPrompt->size(); tank++) {
                //         if(HeCal_E_tankMyPrompt->at(tank)>0.) TOTHeCalLate += HeCal_E_tankMyPrompt->at(tank);
                //         if(HeCal_E_tankMyPrompt->at(tank) > 2) TOTHeCalTankLate++;
                //     }
                //     for(int slab=0; slab<Hod3_E_slabMyPrompt->size(); slab++) {
                //         if(Hod3_E_slabMyPrompt->at(slab)>0.) TOTHod3Prompt += Hod3_E_slabMyPrompt->at(slab);
                //         if(Hod2_E_slabMyPrompt->at(slab)>0.) TOTHod2Prompt += Hod2_E_slabMyPrompt->at(slab);
                //     }

                //     //counting the number of scintillators above the DELAYED threshold
                //     NHod3Late = 0; NHod2Late = 0;
                //     for(int iscint=0; iscint<Hod3_E_slabLate->size(); iscint++){
                //         if(Hod3_E_slabLate->at(iscint)>th_Hod_low+0.2) NHod3Late++;
                //         if(Hod2_E_slabLate->at(iscint)>th_Hod_low+0.2) NHod2Late++;
                //     }

                //     if(Hod3X->size()>0 && Hod2X->size()>0){
                //         double dist = sqrt(pow((Hod2Y->at(0)-Hod3Y->at(0)),2) + pow((Hod2X->at(0)-Hod3X->at(0)),2) + pow((Hod2Z->at(0)-Hod3Z->at(0)),2));
                //         double deltaZ = abs(Hod2Z->at(0)-Hod3Z->at(0));
                //         double cosTheta = deltaZ/dist; 
                //         PathLenght = Hod3Thick/(cosTheta); //cm
                //     }

                //     if(i%2 == 0){
                //         h2_Eprompt_VS_Beta[D]->Fill(TOTHeCalPrompt, beta);
                //         h2_EHecal_VS_Hod2[D]->Fill(TOTHeCalPrompt, TOTHod2Prompt);
                //         h2_EHecal_VS_Hod3[D]->Fill(TOTHeCalPrompt, TOTHod3Prompt/(PathLenght*MIP_norm));

                //         h2_EHecalL_VS_Hod2SlabL[D]->Fill(TOTHeCalLate, NHod2Late);
                //         h2_EHecalL_VS_Hod3SlabL[D]->Fill(TOTHeCalLate, NHod3Late);
                //         h2_EHecalL_VS_HodSlabL[D]->Fill(TOTHeCalLate, NHod2Late+NHod3Late);
                //         h2_EHecalL_VS_HeCalTankL[D]->Fill(TOTHeCalLate, TOTHeCalTankLate);

                //         h2_TOTAL_EHecalL_VS_TotalHit_part[D]->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate);

                //         h2_TOTAL_EHecal_VS_Hod->Fill(TOTHeCalPrompt, TOTHod2Prompt+TOTHod3Prompt);
                //         h2_TOTAL_EHecalL_VS_HodSlabL->Fill(TOTHeCalLate, NHod2Late+NHod3Late);
                //         h2_TOTAL_EHecalL_VS_TotalHit->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate);
                //     }
                //     else{
                //         h2_Eprompt_VS_Beta[P]->Fill(TOTHeCalPrompt, beta, weight);
                //         h2_EHecal_VS_Hod2[P]->Fill(TOTHeCalPrompt, TOTHod2Prompt, weight);
                //         h2_EHecal_VS_Hod3[P]->Fill(TOTHeCalPrompt, TOTHod3Prompt/(PathLenght*MIP_norm), weight);

                //         h2_EHecalL_VS_Hod2SlabL[P]->Fill(TOTHeCalLate, NHod2Late, weight);
                //         h2_EHecalL_VS_Hod3SlabL[P]->Fill(TOTHeCalLate, NHod3Late, weight);
                //         h2_EHecalL_VS_HodSlabL[P]->Fill(TOTHeCalLate, NHod2Late+NHod3Late, weight);
                //         h2_EHecalL_VS_HeCalTankL[P]->Fill(TOTHeCalLate, TOTHeCalTankLate, weight);

                //         h2_TOTAL_EHecalL_VS_TotalHit_part[P]->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate, weight);

                //         h2_TOTAL_EHecal_VS_Hod->Fill(TOTHeCalPrompt, TOTHod2Prompt+TOTHod3Prompt, weight);
                //         h2_TOTAL_EHecalL_VS_HodSlabL->Fill(TOTHeCalLate, NHod2Late+NHod3Late, weight);
                //         h2_TOTAL_EHecalL_VS_TotalHit->Fill(TOTHeCalLate, NHod2Late+NHod3Late+TOTHeCalTankLate, weight);
                //     }
                // } //Delayed selection
        } //end loop on the events
    } //end loop on the files       

    //Saving the histograms 
    TFile *f_out = new TFile ("plots_florence_particles.root", "RECREATE");
    f_out->cd();
    for(int i_hist=0; i_hist<NParticles; i_hist++) {
        //Set all null bins to -1 in order to be displayed
        for(int i=1; i<=h2_Eprompt_VS_Beta[i_hist]->GetNbinsX(); i++){
            for(int j=1; j<=h2_Eprompt_VS_Beta[i_hist]->GetNbinsY(); j++){
                if(h2_Eprompt_VS_Beta[i_hist]->GetBinContent(i,j)==0) h2_Eprompt_VS_Beta[i_hist]->SetBinContent(i,j,-1);
            }
        }
        h2_Eprompt_VS_Beta[i_hist]->Write();


        for(int i=1; i<=h2_EHecal_VS_Hod2[i_hist]->GetNbinsX(); i++){
            for(int j=1; j<=h2_EHecal_VS_Hod2[i_hist]->GetNbinsY(); j++){
                if(h2_EHecal_VS_Hod2[i_hist]->GetBinContent(i,j)==0) h2_EHecal_VS_Hod2[i_hist]->SetBinContent(i,j,-1);
            }
        }
        h2_EHecal_VS_Hod2[i_hist]->Write();

        for(int i=1; i<=h2_EHecal_VS_Hod3[i_hist]->GetNbinsX(); i++){
            for(int j=1; j<=h2_EHecal_VS_Hod3[i_hist]->GetNbinsY(); j++){
                if(h2_EHecal_VS_Hod3[i_hist]->GetBinContent(i,j)==0) h2_EHecal_VS_Hod3[i_hist]->SetBinContent(i,j,-1);
            }
        }
        h2_EHecal_VS_Hod3[i_hist]->Write();

        // h2_EHecalL_VS_Hod2SlabL[i_hist]->Write();
        // h2_EHecalL_VS_Hod3SlabL[i_hist]->Write();
        // h2_EHecalL_VS_HodSlabL[i_hist]->Write();
        // h2_EHecalL_VS_HeCalTankL[i_hist]->Write();

        // h2_TOTAL_EHecalL_VS_TotalHit_part[i_hist]->Write();

        h2_TOTAL_Eprompt_VS_Beta->Add(h2_Eprompt_VS_Beta[i_hist]);

        h2_TOTAL_EHecal_VS_Hod2->Add(h2_EHecal_VS_Hod2[i_hist]);
        h2_TOTAL_EHecal_VS_Hod3->Add(h2_EHecal_VS_Hod3[i_hist]);

        // h2_TOTAL_EHecalL_VS_Hod2SlabL->Add(h2_EHecalL_VS_Hod2SlabL[i_hist]);
        // h2_TOTAL_EHecalL_VS_Hod3SlabL->Add(h2_EHecalL_VS_Hod3SlabL[i_hist]);

        // h2_TOTAL_EHecalL_VS_HeCalTankL->Add(h2_EHecalL_VS_HeCalTankL[i_hist]);
    }

    for(int i=1; i<=h2_TOTAL_Eprompt_VS_Beta->GetNbinsX(); i++){
        for(int j=1; j<=h2_TOTAL_Eprompt_VS_Beta->GetNbinsY(); j++){
            if(h2_TOTAL_Eprompt_VS_Beta->GetBinContent(i,j)==0) h2_TOTAL_Eprompt_VS_Beta->SetBinContent(i,j,-1);
        }
    }
    h2_TOTAL_Eprompt_VS_Beta->Write();

    for(int i=1; i<=h2_TOTAL_EHecal_VS_Hod->GetNbinsX(); i++){
        for(int j=1; j<=h2_TOTAL_EHecal_VS_Hod->GetNbinsY(); j++){
            if(h2_TOTAL_EHecal_VS_Hod->GetBinContent(i,j)==0) h2_TOTAL_EHecal_VS_Hod->SetBinContent(i,j,-1);
        }
    }
    h2_TOTAL_EHecal_VS_Hod->Write();

    for(int i=1; i<=h2_TOTAL_EHecal_VS_Hod2->GetNbinsX(); i++){
        for(int j=1; j<=h2_TOTAL_EHecal_VS_Hod2->GetNbinsY(); j++){
            if(h2_TOTAL_EHecal_VS_Hod2->GetBinContent(i,j)==0) h2_TOTAL_EHecal_VS_Hod2->SetBinContent(i,j,-1);
        }
    }
    h2_TOTAL_EHecal_VS_Hod2->Write();

    for(int i=1; i<=h2_TOTAL_EHecal_VS_Hod3->GetNbinsX(); i++){
        for(int j=1; j<=h2_TOTAL_EHecal_VS_Hod3->GetNbinsY(); j++){
            if(h2_TOTAL_EHecal_VS_Hod3->GetBinContent(i,j)==0) h2_TOTAL_EHecal_VS_Hod3->SetBinContent(i,j,-1);
        }
    }
    h2_TOTAL_EHecal_VS_Hod3->Write();
    // h2_TOTAL_EHecalL_VS_Hod2SlabL->Write();
    // h2_TOTAL_EHecalL_VS_Hod3SlabL->Write();
    // h2_TOTAL_EHecalL_VS_HodSlabL->Write();
    // h2_TOTAL_EHecalL_VS_HeCalTankL->Write();
    // h2_TOTAL_EHecalL_VS_TotalHit->Write();

    f_out->Close();

    return;
}







