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


void trigger_plotter(){

    const int NAntiParticles = 2;
    const int NRates = 8;

    TH1D * h_rate[NRates], *h1_acceptance[NAntiParticles];
    TString rate_names[NRates] = {"h_pRate", "h_He4Rate", "h_C12Rate", "h_eRate", "h_pRateD", "h_He4RateD", "h_C12RateD", "h_eRateD"};
    TString rate_titles[NRates] = {"p prompt", "^{4}He prompt", "^{12}C prompt", "e^{-} prompt",  "p delayed", "^{4}He delayed", "^{12}C delayed", "e^{-} delayed"};

    TString particles_name[NAntiParticles] = {"antiD", "antiP"};
    TString particles[NAntiParticles] = {"#bar{D}", "#bar{p}"};

    enum particle_order {Proton, He4, C12, Electron};
    enum anti_particle_order {Dbar, Pbar};

    const double ProbDelCapture = 0.033;
    const double AnnAfter50nsPbar = 1;//0.64;
    const double AnnAfter50nsDbar = 1;//0.5;

    const double gate = 3.95e-6; // s

    const double Days = 105;
    const double SecondsInADay = 86400;
    const double FlightDuration = Days*SecondsInADay;

    const int nbins = 400;
    const double Ek_min = 10, Ek_max = 100000; // MeV



    const double LabelSize = 0.05;
    const double TitleSize = 0.05;
    const double TitleOffset = 1.2;
    
    TFile *f_in = new TFile ("rate_acceptance.root", "OPEN");


    for(int i_hist=0; i_hist<NRates; i_hist++) {
        h_rate[i_hist] = (TH1D*)f_in->Get(rate_names[i_hist]);
        h_rate[i_hist]->SetLineWidth(2);
        if(i_hist == Proton) h_rate[i_hist]->SetLineColor(kMagenta-2);
        if(i_hist == He4) h_rate[i_hist]->SetLineColor(kCyan-6);
        if(i_hist == C12) h_rate[i_hist]->SetLineColor(kAzure-4);
        if(i_hist == Electron) h_rate[i_hist]->SetLineColor(kGreen+3);

        if(i_hist == NRates-4) h_rate[i_hist]->SetLineColor(kMagenta);
        if(i_hist == NRates-3) h_rate[i_hist]->SetLineColor(kCyan-9);
        if(i_hist == NRates-2) h_rate[i_hist]->SetLineColor(kAzure-7);
        if(i_hist == NRates-1) h_rate[i_hist]->SetLineColor(kGreen+1);

        if(i_hist>Electron) h_rate[i_hist]->SetFillColorAlpha(h_rate[i_hist]->GetLineColor(), 1-(i_hist*0.1));
        h_rate[i_hist]->SetTitle(rate_titles[i_hist]);
        h_rate[i_hist]->Rebin(4);
        if(rate_names[i_hist].Contains("D")) h_rate[i_hist]->SetLineStyle(1);

    }
    for(int i_hist=0; i_hist<NAntiParticles; i_hist++) {
        h1_acceptance[i_hist] = (TH1D*)f_in->Get("h1_acceptance_"+particles_name[i_hist]);
         h1_acceptance[i_hist]->SetLineWidth(2);
        if(i_hist == Dbar) {
            h1_acceptance[i_hist]->SetLineColor(kBlue);
            h1_acceptance[i_hist]->Scale(ProbDelCapture*AnnAfter50nsDbar);
        }
        if(i_hist == Pbar){
            for(int i_bin=0; i_bin<=nbins; i_bin++){
                if(h1_acceptance[i_hist]->GetBinCenter(i_bin)>0.) continue;
                h1_acceptance[i_hist]->SetBinContent(i_bin, 0.);
            }
            h1_acceptance[i_hist]->SetLineColor(kRed);
            h1_acceptance[i_hist]->Scale(ProbDelCapture*AnnAfter50nsPbar);
            h1_acceptance[i_hist]->GetXaxis()->SetRangeUser(10, 1.e+03);
        }

        h1_acceptance[i_hist]->SetFillColorAlpha(h1_acceptance[i_hist]->GetLineColor(), 1 - ( (i_hist+3) *0.2));
        h1_acceptance[i_hist]->SetTitle(particles[i_hist]);
        h1_acceptance[i_hist]->Rebin(2);
    }


    TH1D *h1_DbarSensitivity = GenerateHistogram("h1_Dbar_SenibleFlux", "SenibleFlux_"+particles[Dbar], "E_{k}/N [MeV/c]", "#phi [m^{2} sr s]^{-1}", nbins, Ek_min, Ek_max);
    TH1D *h1_mean = GenerateHistogram("h1_mean", "Mean"+particles[Dbar], "E_{k}/N [MeV/c]", "mean #phi [m^{2} sr s]^{-1}", nbins, Ek_min, Ek_max);

    for(int i_bin=1; i_bin<=nbins; i_bin++){
        h1_DbarSensitivity->SetBinContent(i_bin, 1.);
        h1_DbarSensitivity->SetBinError(i_bin, 1.);
    }
    h1_DbarSensitivity->Rebin(2);
    h1_mean->Rebin(2);
    h1_DbarSensitivity->Scale(1./2.);
    h1_DbarSensitivity->Divide(h1_acceptance[Dbar]);
    h1_DbarSensitivity->Scale(1./FlightDuration);

    gStyle->SetOptStat(0); // no statistics

    const double TopMargin = 0.05, RightMargin = 0.14, LeftMargin = 0.15, BottomMargin = 0.12;

    TCanvas *cCombined = new TCanvas("cCombined", "Acceptances and background rates", 1500, 563);
    cCombined->Divide(2,1, 0.001, 0.001); 
    

    



    //Plotting the rates
    TCanvas *c3 = new TCanvas("c3", "Proton rate", 800, 600);
    c3->cd();
    c3->SetLogx();
    c3->SetLogy();
    c3->SetGridx(); 
    c3->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);

    TString integrals_names[NRates] = {"PROTON prompt ", "He4 prompt    ", "C12 prompt    ", "e- prompt     ",
                                       "PROTON delayed", "He4 delayed   ", "C12 delayed   ", "e- delayed    "};

    double Total_PromptRate = 0, Total_DelayedRate = 0;
    TH1 *h_copy, *h_copyAcc;

    for(int i_plot=0; i_plot<NRates; i_plot++){
        if(i_plot==0) {
            h_rate[i_plot]->GetYaxis()->SetRangeUser(0.15, 1.e+3);
            h_rate[i_plot]->GetXaxis()->SetLabelSize(LabelSize);
            h_rate[i_plot]->GetYaxis()->SetLabelSize(LabelSize);
            h_rate[i_plot]->GetXaxis()->SetTitleSize(TitleSize);
            h_rate[i_plot]->GetYaxis()->SetTitleSize(TitleSize);
            h_rate[i_plot]->GetXaxis()->SetTitleOffset(TitleOffset);
            h_rate[i_plot]->GetYaxis()->SetTitleOffset(TitleOffset);
            h_copy = (TH1*)h_rate[i_plot]->DrawClone("hist ][");
            h_copy->GetXaxis()->SetTitle("E_{k}/N [MeV]");
            cCombined->cd(1);
            h_copy->Draw("hist ][");
            c3->cd();
            // h_rate[i_plot]->GetYaxis()->SetRangeUser(0.1, 2000);
            // h_rate[i_plot]->GetXaxis()->SetRangeUser(20, 1.e+05);
        }
        else {
            h_rate[i_plot]->Draw("hist ][ same");            
            cCombined->cd(1);
            h_rate[i_plot]->Draw("hist ][ same");
            c3->cd();
        }
        printf(integrals_names[i_plot]+"  ----- Total rate: %.5e \n", h_rate[i_plot]->Integral());
        if(i_plot<=Electron) Total_PromptRate+=h_rate[i_plot]->Integral();
        else Total_DelayedRate+=h_rate[i_plot]->Integral();
    }

    printf("Total expected acquired data: %.5e \n", Total_PromptRate*Total_DelayedRate*gate);
    TLegend *leg_rate = c3->BuildLegend();
    h_copy->SetTitle("");
    leg_rate->SetNColumns(4);
    leg_rate->SetFillStyle(0);
    leg_rate->SetLineWidth(0);
    gPad->RedrawAxis();
    
    cCombined->cd(1);
    leg_rate->Draw("same");
    gPad->RedrawAxis();
    c3->SaveAs("MatterRates.pdf");
    

    TCanvas *c4 = new TCanvas("c4", "Antiparticle acceptances", 800, 600);
    c4->cd();
    c4->SetLogx();
    c4->SetLogy();
    c4->SetGridx();
    c4->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);
    for(int i_hist=0; i_hist<NAntiParticles; i_hist++){
        h1_acceptance[i_hist]->GetYaxis()->SetRangeUser(4.e-3, 1.4);
        h1_acceptance[i_hist]->SetMinimum(5.e-3);
        if(i_hist==0) {
            h1_acceptance[i_hist]->GetXaxis()->SetRangeUser(30, 1.e+03);
            h1_acceptance[i_hist]->GetXaxis()->SetLabelSize(LabelSize);
            h1_acceptance[i_hist]->GetYaxis()->SetLabelSize(LabelSize);
            h1_acceptance[i_hist]->GetXaxis()->SetTitleSize(TitleSize);
            h1_acceptance[i_hist]->GetYaxis()->SetTitleSize(TitleSize);
            h1_acceptance[i_hist]->GetXaxis()->SetTitleOffset(TitleOffset);
            h1_acceptance[i_hist]->GetYaxis()->SetTitleOffset(TitleOffset);
            h_copyAcc = (TH1*)h1_acceptance[i_hist]->DrawClone("hist ][");
            h_copyAcc->GetXaxis()->SetTitle("E_{k}/N [MeV]");
            cCombined->cd(2);
            h_copyAcc->Draw("hist ][");
            h_copyAcc->SetTitle("#bar{d}");
            c4->cd();
        }
        else {
            h1_acceptance[i_hist]->Draw("hist ][ same");
            cCombined->cd(2);
            h1_acceptance[i_hist]->Draw("hist ][ same");
            c4->cd();
        }
    }
    TLegend *leg_acc = c4->BuildLegend();
    h_copyAcc->SetTitle("");
    leg_acc->SetMargin(0.5);
    leg_acc->SetFillStyle(0);
    leg_acc->SetLineWidth(0);
    gPad->RedrawAxis();
    cCombined->cd(2);
    leg_acc->Draw("same");
    c4->SaveAs("AntiAcceptances.pdf");



    cCombined->cd(1);
    // gPad->SetGrid();
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);
    gPad->SetTicks(1, 1);
    gPad->SetFrameLineWidth(3);


    cCombined->cd(2);
    // gPad->SetGrid();
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);
    gPad->SetTicks(1, 1);
    gPad->SetFrameLineWidth(3);
    gPad->Update();
    gPad->RedrawAxis();


    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.SetLineWidth(4);
    l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
    // TLine l2;
    // l2.SetLineWidth(3);
    // l2.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
    // l2.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());






    TCanvas *c_sensitivity = new TCanvas("c_sensitivity", "Sensitivity", 800, 600);
    c_sensitivity->cd();
    c_sensitivity->SetLogx();
    c_sensitivity->SetLogy();
    c_sensitivity->SetGridx();
    c_sensitivity->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);
    h1_DbarSensitivity->SetLineWidth(3);
    h1_DbarSensitivity->GetXaxis()->SetRangeUser(30, 1.e+03);
    h1_DbarSensitivity->GetXaxis()->SetLabelSize(LabelSize);
    h1_DbarSensitivity->GetYaxis()->SetLabelSize(LabelSize);
    h1_DbarSensitivity->GetXaxis()->SetTitleSize(TitleSize);
    h1_DbarSensitivity->GetYaxis()->SetTitleSize(TitleSize);
    h1_DbarSensitivity->GetXaxis()->SetTitleOffset(TitleOffset);
    h1_DbarSensitivity->GetYaxis()->SetTitleOffset(TitleOffset);
    h1_DbarSensitivity->Draw("e");


    // TCanvas *c_mean = new TCanvas("c_mean", "mean", 800, 600);
    // c_mean->cd();
    // c_mean->SetLogx();
    // c_mean->SetLogy();
    // c_mean->SetGridx();
    // c_mean->SetGridy();

    // const double mean_min = 30, mean_max = 1.e+03;
    // // const int binNumber = (mean_max-mean_min)/h1_DbarSensitivity->GetBinNumber();

    // h1_DbarSensitivity->SetLineWidth(3);
    // h1_DbarSensitivity->GetXaxis()->SetRangeUser(30, 1.e+03);
    // h1_DbarSensitivity->GetXaxis()->SetLabelSize(LabelSize);
    // h1_DbarSensitivity->GetYaxis()->SetLabelSize(LabelSize);
    // h1_DbarSensitivity->GetXaxis()->SetTitleSize(TitleSize);
    // h1_DbarSensitivity->GetYaxis()->SetTitleSize(TitleSize);
    // h1_DbarSensitivity->GetXaxis()->SetTitleOffset(TitleOffset);
    // h1_DbarSensitivity->GetYaxis()->SetTitleOffset(TitleOffset);
    // h1_DbarSensitivity->Draw("e");

    return;

}







