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
Legge in automatico dei file predefiniti e fornisce dei confronti tra diverse 
distribuzioni; i.e. accettanze, energie nell'He e distribuzioni temporali.
DA USARE per file CON DIVERSE SPECIE DI PARTICELLE.
*/

void comparison_debug()
{
    //gStyle->SetOptStat(0); // no statistics
    
    // ----------------- INPUT FILE ----------------- //
    const int nFiles = 3;
    Color_t colors[nFiles] = {kBlue,kRed,kBlack};
    TString fNames[nFiles] = {"dbar_5_10e7.root","pbar_5_10e7.root","proton_5_10e7.root"};
    // TString fNames[nFiles] = {"dbar_5_10e6_debug.root","pbar_5_10e6_debug.root","proton_5_10e6_debug.root"};
    //proton_10e5_debug.root
    //TString fNames[nFiles] = {"dbar_10e6.root","pbar_10e6.root","proton_10e6.root"};
    TString hNames[nFiles] = {"#bar{d}","#bar{p}","p"};

    // ---- Histograms ---- //
    // MC TRUTH
    TH1D* h_MCEnergy[nFiles]; // energy of the primary
    TH1D* h_Acceptance[nFiles]; // acceptance

    TH1D* h_MCEnergy_c[nFiles-1]; // energy of the capture
    // TH1D* h_MCEnergy_cP[nFiles-1]; // energy of the pbar capture
    // TH1D* h_MCEnergy_DC[nFiles-1]; // energy of the double capture
    //TH1D* h_MCEnergy_cD; // energy of the dbar capture
    
    //He histos
    const int NTanks = 75;
    double hitE[NTanks] ;//= {};
    double hitE_A[NTanks] ;//= {};
    TH1D* h_HeE_tanks[nFiles][NTanks];
    TH1D* h_HeE_tanksA[nFiles-1][NTanks];

    double hitHod2E[NTanks], hitHod3E[NTanks];
    TH1D* h_Hod2E_tanks[nFiles], *h_Hod3E_tanks[nFiles];
    

    // --- Acceptance constants --- //
    
    const float pi = acos(-1.);
    const float Area_mm = 2.3104e+08; // mm^2 (divided by 4) (OLD VALUE 2.3104e+06) 2.3104e+08
    float Area = 4*(Area_mm*1e-6); // m^2
    cout << "Area: " << Area << " m^2" << endl;
    const float PlaneA = pi*Area;
    

    //TFile *f = new TFile("dbar_10e5_debug.root","READ");
    for(int i = 0; i<nFiles; i++){
        TFile *f = new TFile(fNames[i],"READ");
        if (!f) return;

        // --- getting the tree --- //
        TTree *t = (TTree*)f->Get("Hits");
        Long64_t entries = t->GetEntries();
        std::printf("Reading file %s with %lld entries\n",fNames[i].Data(),entries);

        double MCEnergy;
        int HGasN, Hod2N, Hod3N;

            // setting the branches
        t->SetBranchAddress("MCEnergy",&MCEnergy);
        t->SetBranchAddress("HGasN",&HGasN);
        t->SetBranchAddress("Hod2N",&Hod2N);
        t->SetBranchAddress("Hod3N",&Hod3N);

            // Inner scintillator
        vector<double> *Hod2E=0, *Hod2T=0;
        TBranch *b_Hod2E=0, *b_Hod2T=0;
        t->SetBranchAddress("Hod2E",&Hod2E,&b_Hod2E);
        t->SetBranchAddress("Hod2T",&Hod2T,&b_Hod2T);

            // Outer scintillator
        vector<double> *Hod3E=0, *Hod3T=0;
        TBranch *b_Hod3E=0, *b_Hod3T=0;
        t->SetBranchAddress("Hod3E",&Hod3E,&b_Hod3E);
        t->SetBranchAddress("Hod3T",&Hod3T,&b_Hod3T);

            // array of the He scintillator
        vector<double> *HGasE=0, *HGasT=0;
        TBranch *b_HGasE=0, *b_HGasT=0;
        vector<int> *HGasCopyNo=0;
        TBranch *b_HGasCopyNo=0;
        t->SetBranchAddress("HGasE",&HGasE,&b_HGasE);
        t->SetBranchAddress("HGasT",&HGasT,&b_HGasT);
        t->SetBranchAddress("HGasCopyNo",&HGasCopyNo,&b_HGasCopyNo);



        // --- initializing the histograms --- //
        if (i==0 || i==1){
            h_MCEnergy[i] = new TH1D("Primary_"+hNames[i],"Primary "+hNames[i]+"; E_{k} [MeV]/N; Entries",99,10,1000);  // energy of the primary
            h_MCEnergy[i]->Sumw2(); h_MCEnergy[i]->SetLineWidth(2);  h_MCEnergy[i]->SetLineColor(colors[i]);

            h_Acceptance[i] = new TH1D("Acceptance_"+hNames[i],"Acceptance "+hNames[i]+"; E_{k} [MeV]/N; [m^{2} sr]",99,10,1000); // acceptance
            h_Acceptance[i]->Sumw2(); h_Acceptance[i]->SetLineWidth(2);  h_Acceptance[i]->SetLineColor(colors[i]);

            h_MCEnergy_c[i] =  new TH1D("Capture_"+hNames[i],"Capture "+hNames[i]+"; E_{k} [MeV]/N; Entries",99,10,1000);      // energy of the capture
            h_MCEnergy_c[i]->Sumw2(); h_MCEnergy_c[i]->SetLineWidth(2);  h_MCEnergy_c[i]->SetLineColor(colors[i]);

            // h_MCEnergy_cP[i] = new TH1D("#bar{p} captured ("+hNames[i]+")","#bar{p} captured ("+hNames[i]+"); E_{k} [MeV]/N; Entries",99,10,1000);  // energy of the pbar capture
            // h_MCEnergy_cP[i]->Sumw2(); h_MCEnergy_cP[i]->SetLineWidth(2);  h_MCEnergy_cP[i]->SetLineColor(colors[i]);

            // h_MCEnergy_DC[i] = new TH1D("Double capture ("+hNames[i]+")","Double capture ("+hNames[i]+"); E_{k} [MeV]/N; Entries",99,10,1000); // energy of the double capture
            // h_MCEnergy_DC[i]->Sumw2(); h_MCEnergy_DC[i]->SetLineWidth(2);  h_MCEnergy_DC[i]->SetLineColor(kViolet+1);

            // deuteron captures
            if(i==0){
                // h_MCEnergy_cD = new TH1D("#bar{d} captured","#bar{d}; E_{k} [MeV]/N; Entries",99,10,1000);
                // h_MCEnergy_cD->Sumw2(); h_MCEnergy_cD->SetLineWidth(2);  h_MCEnergy_cD->SetLineColor(colors[i]);
            }

        }
        else {
            h_MCEnergy[i] = new TH1D("Primary_"+hNames[i],"Primary "+hNames[i]+"; E_{k} [MeV]/N; Entries",990,2000,10000);  // energy of the primary
            h_MCEnergy[i]->Sumw2(); h_MCEnergy[i]->SetLineWidth(2);  h_MCEnergy[i]->SetLineColor(colors[i]);

            h_Acceptance[i] = new TH1D("Acceptance_"+hNames[i],"Acceptance "+hNames[i]+"; E_{k} [MeV]/N; [m^{2} sr]",990,2000,10000);   // acceptance
            h_Acceptance[i]->Sumw2(); h_Acceptance[i]->SetLineWidth(2);  h_Acceptance[i]->SetLineColor(colors[i]);
        }


        TString HODname = "Hod2E_"+hNames[i]+Form("_%d",i);
        h_Hod2E_tanks[i] = new TH1D(HODname,HODname+"; E_{He} [MeV]; Entries",75,0,150);
        h_Hod2E_tanks[i]->SetLineColor(colors[i]); h_Hod2E_tanks[i]->SetLineWidth(2);

        HODname = "Hod3E_"+hNames[i]+Form("_%d",i);
        h_Hod3E_tanks[i] = new TH1D(HODname,HODname+"; E_{He} [MeV]; Entries",75,0,150);
        h_Hod3E_tanks[i]->SetLineColor(colors[i]); h_Hod3E_tanks[i]->SetLineWidth(2);
            
        for(int j=0; j<NTanks; j++) {
            TString name = "HeE_"+hNames[i]+Form("_%d",j);
            h_HeE_tanks[i][j] = new TH1D(name,name+"; E_{He} [MeV]; Entries",75,0,150);
            h_HeE_tanks[i][j]->SetLineColor(colors[i]); h_HeE_tanks[i][j]->SetLineWidth(2);
            if(i==0 || i==1){
                name = "HeE_"+hNames[i]+Form("_%d_A",j);
                h_HeE_tanksA[i][j] = new TH1D(name,name+"; E_{He} [MeV]; Entries",75,0,250);
                h_HeE_tanksA[i][j]->SetLineColor(colors[i]); h_HeE_tanksA[i][j]->SetLineStyle(2);
                h_HeE_tanksA[i][j]->SetLineWidth(2);
            }
        }
            // array of the capture
        vector<double> *CaptureT=0, *CaptureM=0;
        TBranch *b_CaptureT=0, *b_CaptureM=0;
        t->SetBranchAddress("CaptureT",&CaptureT,&b_CaptureT);
        t->SetBranchAddress("CaptureM",&CaptureM,&b_CaptureM);

        

        // ------------------- READING THE TREE ------------------------ //

        Long64_t counting = 0;
        int doubleCapture = 0;

        for(Int_t entry = 0; entry<entries; entry++){
            for(int tank=0; tank<NTanks; tank++) {
                hitE[tank] = 0;
                hitE_A[tank] = 0;
            }

            t->GetEntry(entry);

            //summing the hits Dbar Pbar
            if(i==0 || i==1){
                for(UInt_t k=0; k<CaptureT->size(); k++){
                    double time = CaptureT->at(k);
                    if(time >0.){
                        for(UInt_t j=0; j<HGasN; j++){
                            int copy = HGasCopyNo->at(j);
                            //printf("primo for su copy\t j: %lld\t HGas: %f \t HGasE: %f\t copy: %d\t Entry: %lld\t CaptureT: %f \n",j,HGasT.At(j), HGasE.At(j),copy, entry, CaptureT.At(j));
                            if(HGasE->at(j)>0.){
                                if(HGasT->at(j)<time) hitE[copy] = hitE[copy] + HGasE->at(j);
                                else hitE_A[copy] = hitE_A[copy] + HGasE->at(j);
                            }
                        }
                        for(int tank=0; tank<NTanks; tank++) {
                            if(hitE[tank]>0) h_HeE_tanks[i][tank]->Fill(hitE[tank]);
                            if(hitE_A[tank]>0) h_HeE_tanksA[i][tank]->Fill(hitE_A[tank]);
                        }
                        counting++;
                    }
                }
            }
            
            
            else{
                for(UInt_t j=0; j<HGasN; j++){
                    int copy = HGasCopyNo->at(j);
                    //printf("primo for su copy\t j: %lld\t HGas: %f \t HGasE: %f\t copy: %d\t Entry: %lld\t CaptureT: %f \n",j,HGasT.At(j), HGasE.At(j),copy, entry, CaptureT.At(j));
                    if(HGasE->at(j)>0.){
                        if(HGasT->at(j)<10000) hitE[copy] = hitE[copy] + HGasE->at(j);
                        else hitE_A[copy] = hitE_A[copy] + HGasE->at(j);
                    }
                }
                for(int tank=0; tank<NTanks; tank++) {
                    if(hitE[tank]>0) h_HeE_tanks[i][tank]->Fill(hitE[tank]);
                    if(hitE_A[tank]>0) h_HeE_tanksA[i][tank]->Fill(hitE_A[tank]);
                }
                counting++;
            }
            

        
            h_MCEnergy[i]->Fill(MCEnergy);
            if(i==0 || i==1){
                for(UInt_t j=0; j<CaptureT->size(); j++){
                    double capT = CaptureT->at(j);
                    if(capT>0.){
                        h_MCEnergy_c[i]->Fill(MCEnergy);
                        // if(CaptureM->at(j)>1000.) h_MCEnergy_cD->Fill(MCEnergy);
                        // else h_MCEnergy_cP[i]->Fill(MCEnergy);
                        if(capT>900000.){
                            //if(CaptureM_v.at(dummy) >= 938.28){
                                //std::printf("\n!!!! --- (entry = %lld) --- !!!!\n", entry);
                                //std::printf("X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",CaptureX.At(dummy),CaptureY.At(dummy),CaptureZ.At(dummy),CaptureE.At(dummy),k,CaptureM_v.at(dummy),CaptureF.At(dummy));
                                //printf("Primary particle charge: %f ; mass: %f\n", MCCharge, MCMass);
                            //}
                            doubleCapture++;
                            // h_MCEnergy_DC[i]->Fill(MCEnergy);
                        }
                    }
                }
            }
                
        }

        std::printf("\n\nNumber of fillments: %lld\n", counting);
        std::printf("number of double capture: %d\n", doubleCapture);



        // --- ACCEPTANCE --- //
    
        //uniform in log10(Ek) from [10,1'000]
        float nGen = 1000000;  
        float GenDensity = nGen/log10(1000./10.);
        float upperLimit = 1000;
        std::printf("GenDensity: %f\n", GenDensity);

        double HighEdge, LowEdge;
        double sel, selErr;

        //uniform in log10(Ek) from [2000,10'000] MeV
        if(i==2) {
            GenDensity = nGen/log10(10000./2000.);
            upperLimit = 10000;
        }
        
        for(int b=1; b<=h_Acceptance[i]->GetNbinsX(); b++){
            if(i==0 || i==1){
                HighEdge = h_MCEnergy_c[i]->GetBinLowEdge(b+1);
                LowEdge = h_MCEnergy_c[i]->GetBinLowEdge(b);
            }
            else{
                HighEdge = h_MCEnergy[i]->GetBinLowEdge(b+1);
                LowEdge = h_MCEnergy[i]->GetBinLowEdge(b);
            }

            if(HighEdge < upperLimit) {
                double log10BinWidth = log10(HighEdge/LowEdge);
                double Gen = GenDensity*log10BinWidth;
                if(i==0 || i==1){
                    sel = h_MCEnergy_c[i]->GetBinContent(b);
                    selErr = h_MCEnergy_c[i]->GetBinError(b);
                }
                else{
                    sel = h_MCEnergy[i]->GetBinContent(b);
                    selErr = h_MCEnergy[i]->GetBinError(b);
                }
                
                //printf("log10 Width: %f\t Gen: %f\t Sel: %f\t LowEdge: %f\t HighEdge: %f\n",log10BinWidth, Gen, sel/Gen, LowEdge, HighEdge);
                
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
    for(int i=0; i<2; i++){
        if(i==0) {
            h_MCEnergy[i]->Draw("e");
            h_MCEnergy_c[i]->Draw("same");
        }
        else {
            h_MCEnergy[i]->Draw("e same");
            h_MCEnergy_c[i]->Draw("same");   
        }
    }
    c1->BuildLegend(); 
    
    printf("number of captured particles: %f\n", h_MCEnergy_c[0]->GetEntries());
    
    
    
    //h_MCEnergy_DC[i]->Draw("same");
    //h_MCEnergy_cD->Draw("same");
    // h_MCEnergy_cP[i]->Draw("same");
    
    
    TCanvas *c2 = new TCanvas("c2","c2",800,800);
    for(int i=0; i<3; i++){
        if(i==0) h_Acceptance[i]->Draw("hist");
        else h_Acceptance[i]->Draw("hist same");
    }
    c2->BuildLegend();

    
    TCanvas* c_tanks[NTanks];
    for(int j=0; j<NTanks; j++){
        TString name = +Form("c_tanks_%d",j);
        c_tanks[j] = new TCanvas(name,name,800,800);
        c_tanks[j]->SetLogy();
        for(int i=2; i>=0; i--){
            if(i==2) {
                h_HeE_tanks[i][j]->DrawNormalized("hist");
            }
            else {
                h_HeE_tanks[i][j]->DrawNormalized("hist same");
                if(i==0 || i ==1) h_HeE_tanksA[i][j]->DrawNormalized("hist same");
            }
        }
        //c_tanks[j]->SetLogy();
        c_tanks[j]->BuildLegend();
        c_tanks[j]->SaveAs("plots/"+name+".png");
        c_tanks[j]->Close();
    }
    
}