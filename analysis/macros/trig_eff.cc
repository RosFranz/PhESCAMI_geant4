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








TH1D* GenerateHistogram(TString Name, TString title, TString XTitle, TString YTitle, 
                        int NbinsX, double xmin, double xmax)
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

    TString Title_tot =  title+";"+XTitle+";"+YTitle;

    TH1D *histogram = new TH1D(Name, Title_tot, NbinsX, &xBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
    histogram->Sumw2();

    return histogram;
}






void trig_eff(){

    // gStyle->SetOptStat(0); // no statistics

    const int NAntiParticles = 2;
    TH1D *h1_acceptance[NAntiParticles];


    const int NFiles = 7;
    TString particles[NFiles] = {"#bar{D}", "#bar{p}", "Deuteron", "proton", "^{4}He", "^{12}C", "e^{-}"};
    TString particles_name[NFiles] = {"antiD", "antiP", "Deuteron", "Proton", "He4", "C12", "electron"};
    // TString names[NFiles] = {"../../sources/ROOT_files/Dbar_10e6_gps_SM1rec_TESTdir__50ns.root", "../../sources/ROOT_files/pbar_10e6_gps_SM1rec_TESTdir__50ns.root", 
    TString names[NFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/pbar_10e7_gps_rec__50ns.root", 
                             "../../sources/ROOT_files/Deuteron_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/proton_10e7_gps_rec__50ns.root", 
                             "../../sources/ROOT_files/He4_10e7_gps_rec__50ns.root", "../../sources/ROOT_files/C12_10e7_gps_rec__50ns.root", 
                             "../../sources/ROOT_files/electron_10e7_gps_rec__50ns.root"}; 

                             //{"../../sources/ROOT_files/Dbar_10e7_rec__50ns.root", "../../sources/ROOT_files/pbar_10e7_rec__50ns.root", 
    //                         "../../sources/ROOT_files/Deuteron_10e7_rec__50ns.root", "../../sources/ROOT_files/proton_10e7_rec__50ns.root", 
    //                         "../../sources/ROOT_files/He4_10e7_rec__50ns.root"}; //old He and proton file with GunParticle

    
    // TString names_GenSpectra[NFiles] = {"../../sources/ROOT_files/Dbar_10e6_STAT_gps_1rec.root", "../../sources/ROOT_files/pbar_10e6_STAT_gps_1rec.root", 
    TString names_GenSpectra[NFiles] = {"../../sources/ROOT_files/Dbar_10e7_gps_rec.root", "../../sources/ROOT_files/pbar_10e7_gps_rec.root", 
                             "../../sources/ROOT_files/Deuteron_10e7_gps_rec.root", "../../sources/ROOT_files/proton_10e7_gps_rec.root", 
                             "../../sources/ROOT_files/He4_10e7_gps_rec.root", "../../sources/ROOT_files/C12_10e7_gps_rec.root", 
                             "../../sources/ROOT_files/electron_10e7_gps_rec.root"};

    enum particle_order {Dbar, Pbar, Deuteron, Proton, He4, C12, Electron};
    
                            
    TH1D *h1_effPrompt[NFiles];
    TH1D *h1_effDelayed[NFiles];
    const int LineWidth = 3;
    Color_t colors[NFiles] = {kRed, kBlue, kGreen+1, kMagenta, kOrange, kBlue-9, kCyan-9} ;
    TPaveText* title;
    const int nbins = 400;
    const double Ek_min = 10, Ek_max = 100000; // MeV
    const double BinWidthLog = (TMath::Log10(Ek_max) - TMath::Log10(Ek_min)) / nbins; // bin width in log scale

    const double GenArea = 231.04; // m^2
    const double PlaneAcceptance = TMath::Pi() * GenArea; // m^2 sr
    // const double NGen = 1e7; // number of generated events (flat in log scale)
    // const double NGenDensity = NGen / log10(Ek_max/Ek_min); // density of generated events in log scale
    // printf("NGenDensity: %.9e\n", NGenDensity);

    TH1D *h1_genEvents[NFiles];
    // TH1D *h1_genEvents = GenerateHistogram("h1_genEvents", "Generated events in each bin", "E_{k}/N [MeV/c]", "Entries", nbins, Ek_min, Ek_max);
    // for(int i=1; i<=h1_genEvents->GetNbinsX(); i++){
    //     h1_genEvents->SetBinContent(i, NGenDensity * log10(h1_genEvents->GetBinLowEdge(i+1)/h1_genEvents->GetBinLowEdge(i)));
    //     h1_genEvents->SetBinError(i, sqrt(h1_genEvents->GetBinContent(i)));
    // }
    


    for(int i=0; i<NFiles; i++) {
        h1_effPrompt[i] = GenerateHistogram("h1_effPrompt_"+particles_name[i], "Prompt_efficiency_"+particles[i], "E_{k}/N [MeV/c]", "#epsilon_{prompt}", nbins, Ek_min, Ek_max);
        h1_effDelayed[i] = GenerateHistogram("h1_effDelayed_"+particles_name[i], "Delayed_efficiency_"+particles[i], "E_{k}/N [MeV/c]", "#epsilon_{delayed}", nbins, Ek_min, Ek_max);
        
        h1_effPrompt[i]->SetLineWidth(LineWidth);
        h1_effPrompt[i]->SetLineColor(colors[i]);
        h1_effPrompt[i]->GetYaxis()->SetTitleOffset(1.0);
        h1_effPrompt[i]->GetXaxis()->SetTitleOffset(1.2);

        h1_effDelayed[i]->SetLineWidth(LineWidth);
        h1_effDelayed[i]->SetLineColor(colors[i]);
        h1_effDelayed[i]->GetYaxis()->SetTitleOffset(1.0);
        h1_effDelayed[i]->GetXaxis()->SetTitleOffset(1.2);

        if(i>NAntiParticles) continue;
        else{
            h1_acceptance[i] = GenerateHistogram("h1_acceptance_"+particles_name[i], "Acceptance_"+particles[i], "E_{k}/N [MeV/c]", "Acceptance [m^{2} sr]", nbins, Ek_min, Ek_max);
            h1_acceptance[i]->SetLineWidth(LineWidth);
            h1_acceptance[i]->SetLineColor(colors[i]);
            h1_acceptance[i]->GetYaxis()->SetTitleOffset(1.0);
            h1_acceptance[i]->GetXaxis()->SetTitleOffset(1.2);
        }
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
    //delayed variables
    int NHod3Late, NHod2Late;
    double Hod3_MaxELate, Hod2_MaxELate, HeCal_MaxELate;
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

        TFile *file_GenSpectra = new TFile(names_GenSpectra[i], "OPEN");
        h1_genEvents[i] = (TH1D*)file_GenSpectra->Get("Gen_spectra");
        

        TFile *file = new TFile(names[i], "OPEN");
        printf("\n Reading file %d/%d: %s\n", i+1, NFiles, names[i].Data());

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
            //clusterizzazione degli scintillatori: HOD3 e HOD2
            t->SetBranchAddress("Hod3_E_slabLate",&Hod3_E_slabLate,&b_Hod3_E_slabLate);
            t->SetBranchAddress("Hod2_E_slabLate",&Hod2_E_slabLate,&b_Hod2_E_slabLate);
        }
        


        for(int entry=0; entry<t->GetEntries(); entry++){
            t->GetEntry(entry);
            NHod3Prompt = 0; NHod2Prompt = 0;
            NHod3Late = 0; NHod2Late = 0;


            //counting the number of scintillators above the PROMPT threshold
            for(int iscint=0; iscint<Hod3_E_slabMyPrompt->size(); iscint++) {
                if(Hod3_E_slabMyPrompt->at(iscint)>th_Hod_high) NHod3Prompt++;
                if(Hod2_E_slabMyPrompt->at(iscint)>th_Hod_high) NHod2Prompt++;
            }
            


            // PROMPT SELECTION
            if( (Hod3_MaxEMyprompt > th_Hod_high) && 
                (Hod2_MaxEMyprompt > th_Hod_high) && 
                (HeCal_MaxEMyprompt > th_HeCal) && 
                (NHod3Prompt <= th_HodHit) && 
                (NHod2Prompt <= th_HodHit) &&  
                (NHod3Prompt > 0) && 
                (NHod2Prompt > 0) &&
                (NHod3Prompt+NHod2Prompt) <= th_HodHitTot ){
                    h1_effPrompt[i]->Fill(MCEnergy);

                    if(i==Dbar || i==Pbar){
                        //counting the number of scintillators above the DELAYED threshold
                        for(int iscint=0; iscint<Hod3_E_slabLate->size(); iscint++){
                            if(th_Hod_high > Hod3_E_slabLate->at(iscint) && Hod3_E_slabLate->at(iscint)>th_Hod_low) NHod3Late++;
                            if(th_Hod_high > Hod2_E_slabLate->at(iscint) && Hod2_E_slabLate->at(iscint)>th_Hod_low) NHod2Late++;
                        }



                        // DELAYED SELECTION
                        if( (Hod3_MaxELate > th_Hod_low) && 
                            (Hod2_MaxELate > th_Hod_low) && 
                            (HeCal_MaxELate > th_HeCal) && 
                            (NHod3Late > th_HodHit) && 
                            (NHod2Late > th_HodHit) ){
                            // && (NHod3Late+NHod2Late) < th_HodHitTot ){
                                h1_effDelayed[i]->Fill(MCEnergy);

                                h1_acceptance[i]->Fill(MCEnergy);
                        }
                    }
                }
            
            if(i==Deuteron || i==Proton || i==He4 || i==C12 || i==Electron) {   // only for particles

                NHod3Prompt = 0;
                NHod2Prompt = 0;

                 //counting the number of scintillators above the DELAYED threshold
                for(int iscint=0; iscint<Hod3_E_slabMyPrompt->size(); iscint++) {
                    if(th_Hod_high > Hod3_E_slabMyPrompt->at(iscint) && Hod3_E_slabMyPrompt->at(iscint) > th_Hod_low) NHod3Prompt++;
                    if(th_Hod_high > Hod2_E_slabMyPrompt->at(iscint) && Hod2_E_slabMyPrompt->at(iscint) > th_Hod_low) NHod2Prompt++;
                }

                // DELAYED SELECTION
                if( (Hod3_MaxEMyprompt > th_Hod_low) && 
                (Hod2_MaxEMyprompt > th_Hod_low) &&     
                (HeCal_MaxEMyprompt > th_HeCal) && 
                (NHod3Prompt > th_HodHit) && 
                (NHod2Prompt > th_HodHit) ){
                    // && (NHod3Late+NHod2Late) < th_HodHitTot ){
                        h1_effDelayed[i]->Fill(MCEnergy);
                }
            }
            
        } //end loop on the events

        h1_effPrompt[i]->Divide(h1_genEvents[i]);
        h1_effDelayed[i]->Divide(h1_genEvents[i]);
        if(i<NAntiParticles) {
            h1_acceptance[i]->Divide(h1_genEvents[i]);
            h1_acceptance[i]->Scale(PlaneAcceptance);
        }
    } //end loop on the files

    TCanvas *c_test = new TCanvas("c_test", "c_test", 800, 600);
    c_test->SetLogx();
    for(int i_hist = 0; i_hist<NFiles; i_hist++){
        if(i_hist == 0) h1_genEvents[0]->Draw("e");
        else h1_genEvents[1]->Draw("e same");
    }

   

    TCanvas *c = new TCanvas("c", "prompt_efficiency", 800, 800);
    c->SetLogx();
    c->SetLogy();
    c->SetGrid();
    for(int i=0; i<NFiles; i++) {
        // h1_effDelayed[i]->GetYaxis()->SetRangeUser(9e-8,3e-4);
        h1_effPrompt[i]->GetYaxis()->SetRangeUser(2.e-6,2e-1);
        h1_effPrompt[i]->SetTitle(particles[i].Data());

        if(i==0) h1_effPrompt[i]->Draw("hist");
        else h1_effPrompt[i]->Draw("hist same");

        if(i!=NFiles-1) continue;
        c->BuildLegend();
        h1_effPrompt[0]->SetTitle("Prompt efficiency");
    }

    TCanvas *cDel = new TCanvas("cDel", "delayed_efficiency", 800, 800);
    cDel->SetLogx();
    cDel->SetLogy();
    cDel->SetGrid();
    for(int ifile=0; ifile<NFiles; ifile++) {
        h1_effDelayed[ifile]->GetYaxis()->SetRangeUser(2.e-6,2.e-1);
        h1_effDelayed[ifile]->SetTitle(particles[ifile].Data());

        if(ifile==0) h1_effDelayed[ifile]->Draw("hist");
        else h1_effDelayed[ifile]->Draw("hist same");

        if(ifile!=NFiles-1) continue;
        cDel->BuildLegend();
        h1_effDelayed[0]->SetTitle("Delayed efficiency");
    }



    // // ----- RATE ESTIMATION -----

    //histo for rates
    const int NRates = 8;
    TH1D *h_rate[NRates];

    TString rate_names[NRates] = {"h_pRate", "h_He4Rate", "h_C12Rate", "h_eRate", "h_pRateD", "h_He4RateD", "h_C12RateD", "h_eRateD"};
    TString rate_titles[NRates] = {"proton PROMPT", "^{4}He PROMPT", "^{12}C PROMPT", "e^{-} PROMPT",  "proton DELAYED", "^{4}He DELAYED", "^{12}C DELAYED", "e^{-} DELAYED"};
    Color_t rate_colors[NRates] = {colors[Proton], colors[He4], colors[C12], colors[Electron], static_cast<Color_t>( colors[Proton]-2 ), static_cast<Color_t>( colors[He4]+2 ), static_cast<Color_t>( colors[C12]+3 ), static_cast<Color_t>( colors[Electron]+3 )};
    const int LineStyle = 2;

    for(int i_rate=0; i_rate<NRates; i_rate++){
        h_rate[i_rate] = GenerateHistogram(rate_names[i_rate], rate_titles[i_rate], "E_{k}/N [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min, Ek_max);
        h_rate[i_rate]->SetLineColor(rate_colors[i_rate]);
        h_rate[i_rate]->SetLineWidth(LineWidth);
        if(rate_names[i_rate].Contains("D")) h_rate[i_rate]->SetLineStyle(LineStyle);
    }

    //getting the flux TGraph
    const int NFluxes = NRates/2;
    TGraph *g_fluxes[NFluxes];
    TString txt_names[NFluxes] = {"protoni_100GeV.txt", "He_TOTAL_100GeV.txt", "C_TOTAL_100GeV.txt", "elettroni_100GeV.txt"};
    TString flux_titles[NFluxes] = {"p", "He", "C", "e^{-}"};
    TString X_axis = "E_{k}/N [MeV/c]";
    TString Y_axis = "Flux [s^{-1} m^{-2} sr^{-1} (GeV)^{-1} ]";
    const double MarkerSize = 1.2;
    const int MarkerStyle = 20;

    for(int i_flux=0; i_flux<NFluxes; i_flux++){
        g_fluxes[i_flux] = new TGraph(txt_names[i_flux], "%lg %lg");
        g_fluxes[i_flux]->Sort();
        g_fluxes[i_flux]->SetLineWidth(LineWidth);
        g_fluxes[i_flux]->SetLineStyle(LineStyle);
        g_fluxes[i_flux]->SetMarkerColor(rate_colors[i_flux]);
        g_fluxes[i_flux]->SetMarkerStyle(MarkerStyle);
        g_fluxes[i_flux]->SetMarkerSize(MarkerSize);

        g_fluxes[i_flux]->SetTitle(flux_titles[i_flux]);
        g_fluxes[i_flux]->GetXaxis()->SetTitle(X_axis);
        g_fluxes[i_flux]->GetYaxis()->SetTitle(Y_axis);
    }

    TGraph *g_He3flux = new TGraph("He3.txt", "%lg %lg");
    TGraph *g_He4flux = new TGraph("He4.txt", "%lg %lg");

    g_He3flux->Sort();
    g_He4flux->Sort();

    g_He3flux->SetLineWidth(LineWidth);
    g_He4flux->SetLineWidth(LineWidth);

    g_He3flux->SetLineColor(colors[He4]-3);
    g_He4flux->SetLineColor(colors[He4]+3);

    g_He3flux->SetMarkerColor(colors[He4]-3);
    g_He4flux->SetMarkerColor(colors[He4]+3);

    g_He3flux->SetLineStyle(2);
    g_He4flux->SetLineStyle(2);

    g_He3flux->SetMarkerStyle(20);
    g_He4flux->SetMarkerStyle(20);

    g_He3flux->SetMarkerSize(1.2);
    g_He4flux->SetMarkerSize(1.2);

    

    TCanvas *c2 = new TCanvas("c2", "Fluxes", 800, 600);
    c2->SetLogx();
    c2->SetLogy();
    g_He3flux->SetTitle("^{3}He");
    g_He4flux->SetTitle("^{4}He");
    g_fluxes[0]->GetYaxis()->SetRangeUser(1.e-4,8e+3);
    for(int i_flux=0; i_flux<NFluxes; i_flux++){
        if(i_flux==0) g_fluxes[i_flux]->Draw("apl");
        else g_fluxes[i_flux]->Draw("pl same");
    }
    g_He3flux->Draw("pl same");
    g_He4flux->Draw("pl same");
    c2->BuildLegend();


    // //fit function 
    // TF1 *f_fitP = new TF1("f_fitP", "[0]*pow(x, [1])", 4000, 10000);
    // TF1 *f_fitHe4 = new TF1("f_fitH4", "[0]*pow(x, [1])", 4000, 10000);
    // TF1 *f_fitP_D = new TF1("f_fitP_He4_D", "[0]*pow(x, [1])", 4000, 10000);
    // TF1 *f_fitHe4_D = new TF1("f_fitP_He4_D", "[0]*pow(x, [1])", 4000, 10000);
    // f_fitP->SetParameters(4.77761e+07, -1.65490e+00);
    // f_fitHe4->SetParameters( 1.78256e+08, -2.19514e+00);
    // double integralP = 0, integralHe4 = 0;

    int eff_counter;
    eff_counter = (Proton > He4) ? He4 : Proton;

    const int rate_bias = 4;


    //cycle on the different fluxes
    for(int i_flux=0; i_flux<NFluxes; i_flux++){

        //cycle on the bins
        for(int i_bin=1; i_bin<=nbins; i_bin++){

            double binCenter = h1_effPrompt[eff_counter]->GetBinCenter(i_bin)/1000.;        //GeV
            double binWidth = h1_effPrompt[eff_counter]->GetBinWidth(i_bin)/1000.;    //GeV
            double fluxValue = g_fluxes[i_flux]->Eval(binCenter);
            double effValue = h1_effPrompt[eff_counter]->GetBinContent(i_bin);
            double rate = (binWidth * PlaneAcceptance * effValue) * fluxValue;
            // printf("Rate: %.9e ; binCenter: %.2e ; binWidth: %.2e ; efficiency: %.9e ; AcceptancePlane: %.2e ; flux: %.2e  \n", rate, binCenter, binWidth, effValue, PlaneAcceptance, fluxValue);
            h_rate[i_flux]->SetBinContent(i_bin, rate);
            h_rate[i_flux]->SetBinError(i_bin, 0.05*rate); // 5% error (completely wrong but just to put some number)

            //delayed proton
            binCenter = h1_effDelayed[eff_counter]->GetBinCenter(i_bin)/1000.;  //GeV
            binWidth = h1_effDelayed[eff_counter]->GetBinWidth(i_bin)/1000.;    //GeV
            fluxValue = g_fluxes[i_flux]->Eval(binCenter);
            effValue = h1_effDelayed[eff_counter]->GetBinContent(i_bin);
            rate = (binWidth * PlaneAcceptance * effValue) * fluxValue;
            h_rate[i_flux+rate_bias]->SetBinContent(i_bin, rate);
            h_rate[i_flux+rate_bias]->SetBinError(i_bin, 0.05*rate); // 5% error (completely wrong but just to put some number)

            //how to extrapolate the fit
            // if(i<=h1_effPrompt[0]->GetNbinsX()){
            // }
            // else{
                // if(i==h1_effPrompt[0]->GetNbinsX()+1){
                //     integralP = h_pRate->Integral();
                //     integralHe4 = h_He4Rate->Integral();
                //     h_pRate->Fit("f_fitP", "QR");
                //     h_He4Rate->Fit("f_fitH4", "QR");
                // } 
                // double binCenter = h_pRate->GetBinCenter(i);        //MeV
                // double binValue = f_fitP->Eval(binCenter);
                // h_pRate->SetBinContent(i, binValue);
                // h_pRate->SetBinError(i, 0.05*binValue); // 5% error (completely wrong but just to put some number)

                // binCenter = h_He4Rate->GetBinCenter(i);        //MeV
                // binValue = f_fitHe4->Eval(binCenter);
                // h_He4Rate->SetBinContent(i, binValue);
                // h_He4Rate->SetBinError(i, 0.05*binValue); // 5% error (completely wrong but just to put some number)
            // }
        } // end of the bin cycle

        eff_counter++;

    } // end of the flux cycle



    //Plotting the rates
    TCanvas *c3 = new TCanvas("c3", "Proton rate", 800, 600);
    c3->cd();
    c3->SetLogx();
    c3->SetLogy();
    c3->SetGridx(); 
    c3->SetGridy();

    TString integrals_names[NRates] = {"PROTON prompt ", "He4 prompt    ", "C12 prompt    ", "e- prompt     ",
                                       "PROTON delayed", "He4 delayed   ", "C12 delayed   ", "e- delayed    "};

    for(int i_plot=0; i_plot<NRates; i_plot++){
        if(i_plot==0) h_rate[i_plot]->Draw("hist");
        else h_rate[i_plot]->Draw("hist same");
        printf(integrals_names[i_plot]+"  ----- Total rate: %.5e \n", h_rate[i_plot]->Integral());
    }
    c3->BuildLegend();
    h_rate[0]->SetTitle("Expected rates");


    TCanvas *c4 = new TCanvas("c4", "Antiparticle acceptances", 800, 600);
    c4->cd();
    c4->SetLogx();
    c4->SetLogy();
    c4->SetGridx();
    c4->SetGridy();
    for(int i_hist=0; i_hist<NAntiParticles; i_hist++){
        if(i_hist==0) h1_acceptance[i_hist]->Draw("hist");
        else h1_acceptance[i_hist]->Draw("hist same");
    }
    c4->BuildLegend();
    h1_acceptance[0]->SetTitle("Acceptances");

    //Saving the histograms 
    TFile *f_out = new TFile ("rate_acceptance.root", "RECREATE");
    f_out->cd();
    for(int i_hist=0; i_hist<NRates; i_hist++) h_rate[i_hist]->Write();
    for(int i_hist=0; i_hist<NFiles; i_hist++) {
        h1_effPrompt[i_hist]->Write();
        h1_effDelayed[i_hist]->Write();
    }
    for(int i_g=0; i_g<NFluxes; i_g++) g_fluxes[i_g]->Write();
    for(int i_hist=0; i_hist<NAntiParticles; i_hist++) h1_acceptance[i_hist]->Write();
    f_out->Close();

    return;

}







