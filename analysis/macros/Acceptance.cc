#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <map>
#include <TH2.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <iostream>
using namespace std;

/*
Programma senza richiesta di argomenti in ingresso.
Legge in automatico dei file predefiniti e fornisce le accettanze
*/

void Acceptance()
{
    //gStyle->SetOptStat(0); // no statistics
    
    // ----------------- INPUT FILE ----------------- //
    const int nFiles = 4;
    Color_t colors[nFiles] = {kBlue,kRed,kRed-3,kBlack};
    TString fNames[nFiles] = {"../analysis_data/dbar_10e8_filtred.root", "../analysis_data/pbar_10e8_filtred.root", "../analysis_data/proton_10e7_REC.root", "../analysis_data/deuteron_10e7_REC.root"};
    TString hNames[nFiles] = {"#bar{D}","#bar{p}","p","d"};

    // ---- Histograms ---- //
    // MC TRUTH
    TH1D* h_MCEnergy[nFiles]; // energy of the primary
    TH1D* h_Acceptance[nFiles]; // acceptance

    // TH1D* h_MCEnergy_c[nFiles-1]; // energy of the capture
    
    //He histos
    const int NTanks = 75;
    double hitE[NTanks] ;//= {};
    double hitE_A[NTanks] ;//= {};
    double maxHeCal=0, maxHeCal_A=0;

    //Plastic scintillators
    const int NScint = 64;
    double hitScint2[NScint], hitScint3[NScint];
    double hitHod2E, hitHod3E;

    //Thresholds
    double thr = 20; // HeCal [MeV]
    double thrS = 6.0 ; // Scintillators [MeV]
    double thrH = 2.5; // Hits in scintillators
    

    // --- Acceptance constants --- //
    
    const float pi = acos(-1.);
    const float Area_mm = 2.3104e+08; // mm^2 2.3104e+08
    float Area = (Area_mm*1e-6); // m^2
    cout << "Area: " << Area << " m^2" << endl;
    const float PlaneA = pi*Area;
    
    // CYCLE ON THE FILES 
    for(int i = 0; i<nFiles; i++){
        TFile *f = new TFile(fNames[i],"READ");
        if (!f) return;

        // --- initializing histograms --- //
        if(i==0||i==1) h_MCEnergy[i] = new TH1D(hNames[i]+"_energy",hNames[i]+"; E_{k} [MeV]/N; Entries",180,10,1000); 
        else h_MCEnergy[i] = new TH1D(hNames[i]+"_energy",hNames[i]+"; E_{k} [MeV]/N; Entries",1980,10,10000); 
        h_MCEnergy[i]->Sumw2(); //h_MCEnergy[i]->SetMarkerStyle(8);  h_MCEnergy[i]->SetMarkerColor(colors[i]); h_MCEnergy[i]->SetMarkerSize(0.3);
        h_MCEnergy[i]->SetLineColor(colors[i]); h_MCEnergy[i]->SetLineWidth(2); //h_MCEnergy[i]->SetFillColor(colors[i]);

        // h_MCEnergy_c[i] =  new TH1D(hNames[i]+"_captured",hNames[i]+"; E_{k} [MeV]/N; Entries",990,10,10000);
        // h_MCEnergy_c[i]->Sumw2(); //h_MCEnergy_c[i]->SetMarkerStyle(8);  h_MCEnergy_c[i]->SetMarkerColor(colors[i]); h_MCEnergy_c[i]->SetMarkerSize(0.3);
        // h_MCEnergy_c[i]->SetLineColor(colors[i]); h_MCEnergy_c[i]->SetLineWidth(2); h_MCEnergy_c[i]->SetFillColor(colors[i]);

        if(i==0||i==1) h_Acceptance[i] = new TH1D(hNames[i]+"_acceptance",hNames[i]+"; E_{k} [MeV]/N; [m^{2} sr]",180,10,1000);
        else h_Acceptance[i] = new TH1D(hNames[i]+"_acceptance",hNames[i]+"; E_{k} [MeV]/N; [m^{2} sr]",1980,10,10000);
        h_Acceptance[i]->Sumw2(); //h_Acceptance[i]->SetMarkerStyle(8);  h_Acceptance[i]->SetMarkerColor(colors[i]); h_Acceptance[i]->SetMarkerSize(0.3);
        h_Acceptance[i]->SetLineColor(colors[i]); h_Acceptance[i]->SetLineWidth(2); //h_Acceptance[i]->SetFillColor(colors[i]);

        // --- getting the tree --- //
        TTree *t = (TTree*)f->Get("Hits");
        Long64_t entries = t->GetEntries();
        std::printf("Reading file %s with %lld entries\n",fNames[i].Data(),entries);

        double MCEnergy;
        int HGasN, Hod2N, Hod3N, EventID;

            // setting the branches
        t->SetBranchAddress("MCEnergy",&MCEnergy);
        t->SetBranchAddress("HGasN",&HGasN);
        t->SetBranchAddress("Hod2N",&Hod2N);
        t->SetBranchAddress("Hod3N",&Hod3N);
        t->SetBranchAddress("EventID",&EventID);

            // Inner scintillator
        vector<double> *Hod2E=0, *Hod2T=0;
        vector<int> *Hod2CopyNo=0;
        TBranch *b_Hod2E=0, *b_Hod2T=0, *b_Hod2CopyNo=0;
        t->SetBranchAddress("Hod2E",&Hod2E,&b_Hod2E);
        t->SetBranchAddress("Hod2T",&Hod2T,&b_Hod2T);
        t->SetBranchAddress("Hod2CopyNo",&Hod2CopyNo,&b_Hod2CopyNo);

             // Outer scintillator
        vector<double> *Hod3E=0, *Hod3T=0;
        vector<int> *Hod3CopyNo=0;
        TBranch *b_Hod3E=0, *b_Hod3T=0, *b_Hod3CopyNo=0;
        t->SetBranchAddress("Hod3E",&Hod3E,&b_Hod3E);
        t->SetBranchAddress("Hod3T",&Hod3T,&b_Hod3T);
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
        

        // ------------------- READING THE TREE ------------------------ //
        int over=0, overS=0, belowH=0;
        int NHit3=0, NHit2=0;

        for(Int_t entry = 0; entry<entries; entry++){

            for(int tank=0; tank<NTanks; tank++) {
                hitE[tank] = 0; hitE_A[tank] = 0;
            }
            for(int hod=0; hod<NScint; hod++) {
                hitScint2[hod] = 0; hitScint3[hod] = 0;
            }
            maxHeCal = 0; maxHeCal_A = 0;
            hitHod3E = 0; hitHod2E = 0;
            NHit3=0, NHit2=0;

            t->GetEntry(entry);

            //Dbar Pbar
            if(i==0 || i==1){
                for(UInt_t k=0; k<CaptureT->size(); k++){
                    double time = CaptureT->at(k);
                    if(time >0.){
                       for(UInt_t j=0; j<HGasN; j++){
                            int copy = HGasCopyNo->at(j);
                            if(HGasE->at(j)>0.){
                                if(HGasT->at(j)<time) hitE[copy] = hitE[copy] + HGasE->at(j);
                            }
                        }
                        for(int tank=0; tank<NTanks; tank++){
                            if(hitE[tank]>maxHeCal) maxHeCal = hitE[tank];
                        }
                        if(maxHeCal>thr){
                            for(UInt_t j=0; j<Hod3N; j++){
                                double Hod3Copy = Hod3CopyNo->at(j)-NScint;
                                if(Hod3Copy>=0 && Hod3Copy<NScint){
                                    int copy3 = Hod3Copy;
                                    if(Hod3T->at(j)<2 && Hod3E->at(j)!=0){
                                        hitScint3[copy3] = hitScint3[copy3] + Hod3E->at(j);
                                        NHit3++;
                                    }
                                }
                                if(Hod3T->at(j)>2 && Hod3T->at(j)<100 && Hod3E->at(j)!=0) NHit3++;
                            }
                            for(UInt_t j=0; j<Hod2N; j++){
                                double Hod2Copy = Hod2CopyNo->at(j)-NScint;
                                if(Hod2Copy>=0 && Hod2Copy<NScint){
                                    int copy2 = Hod2Copy;
                                    if(Hod2T->at(j)<7 && Hod2E->at(j)!=0) {
                                        hitScint2[copy2] = hitScint2[copy2] + Hod2E->at(j);
                                        NHit2++;
                                    }
                                }
                                if(Hod2T->at(j)>7 && Hod2T->at(j)<100 && Hod2E->at(j)!=0) NHit2++;
                            }
                            for(int hod=0; hod<NScint; hod++){
                                if(hitScint3[hod]>0 && hitScint3[hod]>hitHod3E) hitHod3E = hitScint3[hod];
                                if(hitScint2[hod]>0 && hitScint2[hod]>hitHod2E) hitHod2E = hitScint2[hod];
                            }
                            double sum = hitHod2E+hitHod3E;
                            if(sum>thrS){
                                int HitSum = NHit2+NHit3;
                                if(HitSum<thrH) h_MCEnergy[i]->Fill(MCEnergy);
                            }
                        }
                        break;
                    }
                }
            }
            
            
            else{
                for(UInt_t j=0; j<HGasN; j++){
                    int copy = HGasCopyNo->at(j);
                    if(HGasE->at(j)>0.){
                        if(HGasT->at(j)<100) hitE[copy] = hitE[copy] + HGasE->at(j);
                    }
                }
                for(int tank=0; tank<NTanks; tank++) {
                    if(hitE[tank]>maxHeCal) maxHeCal = hitE[tank];
                }
                if(maxHeCal>thr) {
                    for(UInt_t j=0; j<Hod3N; j++){
                        double Hod3Copy = Hod3CopyNo->at(j)-NScint;
                        if(Hod3Copy>=0 && Hod3Copy<NScint){
                            int copy3 = Hod3Copy;
                            if(Hod3T->at(j)<2. && Hod3E->at(j)!=0.){
                                hitScint3[copy3] = hitScint3[copy3] + Hod3E->at(j);
                                NHit3++;
                            }
                        }
                        if(Hod3T->at(j)>2 && Hod3T->at(j)<100 && Hod3E->at(j)!=0) NHit3++;
                    }
                    for(UInt_t j=0; j<Hod2N; j++){
                        double Hod2Copy = Hod2CopyNo->at(j)-NScint;
                        if(Hod2Copy>=0 && Hod2Copy<NScint){
                            int copy2 = Hod2Copy;
                            if(Hod2T->at(j)<7. && Hod2E->at(j)!=0.){
                                hitScint2[copy2] = hitScint2[copy2] + Hod2E->at(j);
                                NHit2++;
                            }
                        }
                        if(Hod2T->at(j)>7 && Hod2T->at(j)<100 && Hod2E->at(j)!=0) NHit2++;
                    }
                    for(int hod=0; hod<NScint; hod++){
                        if(hitScint3[hod]>hitHod3E) hitHod3E = hitScint3[hod];
                        if(hitScint2[hod]>hitHod2E) hitHod2E = hitScint2[hod];
                    }
                    double sum = hitHod2E+hitHod3E;
                    if(sum>thrS) {
                        int HitSum = NHit2+NHit3;
                        if(HitSum<thrH) h_MCEnergy[i]->Fill(MCEnergy);
                    }
                }
            }   
        }

        printf("\t sel. events (%f) \n", h_MCEnergy[i]->GetEntries());

        // --- ACCEPTANCE --- //
    
        //uniform in log10(Ek) from [10,1'000]
        float nGen = 50000000;  
        float GenDensity = nGen/log10(1000./10.);
        float upperLimit = 1000;
        std::printf("GenDensity: %f\n", GenDensity);

        double HighEdge, LowEdge;
        double sel, selErr;

        //uniform in log10(Ek) from [2000,10'000] MeV
        if(i==2 || i==3) {
            nGen = 10000000;
            GenDensity = nGen/log10(10000./10.);
            upperLimit = 10000;
        }
        
        for(int b=1; b<=h_Acceptance[i]->GetNbinsX(); b++){

            HighEdge = h_MCEnergy[i]->GetBinLowEdge(b+1);
            LowEdge = h_MCEnergy[i]->GetBinLowEdge(b);
            
            if(HighEdge <= upperLimit) {
                double log10BinWidth = log10(HighEdge/LowEdge);
                double Gen = GenDensity*log10BinWidth;
 
                sel = h_MCEnergy[i]->GetBinContent(b);
                selErr = h_MCEnergy[i]->GetBinError(b);
                
                double acc = (sel/Gen)*PlaneA;
                double accErr = (selErr/Gen)*PlaneA;
                
                h_Acceptance[i]->SetBinContent(b,acc);
                h_Acceptance[i]->SetBinError(b,accErr);
                nGen = nGen - Gen;
            }
            else{
                h_Acceptance[i]->SetBinContent(b,0);
                h_Acceptance[i]->SetBinError(b,0);
            }
        }
        std::printf("CrossCheck ... %f\n", nGen);
    }
    
    
    
    
    // ------------------- PLOTTING ------------------- //
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->SetLogx(); c1->SetLogy();
    for(int i=0; i<nFiles; i++){
        if(i==0) {
            h_MCEnergy[i+nFiles-1]->DrawNormalized("hist");
        }
        else h_MCEnergy[i-1]->DrawNormalized("hist same");
        printf("number of captured particles: %f\n", h_MCEnergy[i]->GetEntries());
    }
    c1->BuildLegend(); 
    
    
    TCanvas *c2 = new TCanvas("c2","c2",800,800);
    c2->SetLogx(); c2->SetLogy(); 
    for(int i=0; i<nFiles; i++){
        //patch SIF 2023
        if(i==0) h_Acceptance[i+1]->Draw("hist e");
        if(i==1) h_Acceptance[i-1]->Draw("hist e same");

        //codice precedente
        // if(i==0){
        //     //h_Acceptance[i+nFiles-1]->GetYaxis()->SetRangeUser(1e-3,7e0);
        //     h_Acceptance[i+nFiles-1]->Draw("hist e");
        // } 
        // else h_Acceptance[i-1]->Draw("hist e same");
    }
    c2->BuildLegend();
    
}