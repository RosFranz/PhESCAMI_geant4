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







void IFAE_plot(){

    gStyle->SetOptStat(0); // no statistics
    gStyle->SetFrameLineWidth(4);

    const int NParticles = 2;
    TString particles[NParticles] = {"#bar{D}", "#bar{p}"};
    TString particles_name[NParticles] = {"antiD", "antiP"};

    enum particle_order {Dbar, Pbar};

    // TFile *f_in = new TFile ("plots_florence_test.root", "OPEN");
    TFile *f_in = new TFile ("plots_florence_testdir.root", "OPEN");

    const double label_size = 0.05;
    const double offset = 1.2, Loffset = 1.e-2;

                            
    TH2D *h2_TOTAL_Eprompt_VS_Beta = (TH2D*)f_in->Get("h2_TOTAL_Eprompt_VS_Beta");
    h2_TOTAL_Eprompt_VS_Beta->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_Eprompt_VS_Beta->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_Eprompt_VS_Beta->Smooth(1);


    TH2D * h2_TOTAL_EHecal_VS_Hod2 = (TH2D*)f_in->Get("h2_TOTAL_EHecal_VS_Hod2"); 
    h2_TOTAL_EHecal_VS_Hod2->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecal_VS_Hod2->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod2->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecal_VS_Hod2->Smooth(1);


    TH2D * h2_TOTAL_EHecal_VS_Hod3 = (TH2D*)f_in->Get("h2_TOTAL_EHecal_VS_Hod3"); 
    h2_TOTAL_EHecal_VS_Hod3->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetTitle("dE/dx [MIPs]");
    h2_TOTAL_EHecal_VS_Hod3->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecal_VS_Hod3->Smooth(1);



    TH2D * h2_TOTAL_EHecal_VS_Hod = (TH2D*)f_in->Get("h2_TOTAL_EHecal_VS_Hod");
    h2_TOTAL_EHecal_VS_Hod->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecal_VS_Hod->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecal_VS_Hod->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecal_VS_Hod->Smooth(1);



    TH2D * h2_TOTAL_EHecalL_VS_Hod2SlabL = (TH2D*)f_in->Get("h2_TOTAL_EHecalL_VS_Hod2SlabL");
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecalL_VS_Hod2SlabL->Smooth(1);


    TH2D * h2_TOTAL_EHecalL_VS_Hod3SlabL = (TH2D*)f_in->Get("h2_TOTAL_EHecalL_VS_Hod3SlabL");
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecalL_VS_Hod3SlabL->Smooth(1);


    TH2D * h2_TOTAL_EHecalL_VS_HodSlabL = (TH2D*)f_in->Get("h2_TOTAL_EHecalL_VS_HodSlabL");
    h2_TOTAL_EHecalL_VS_HodSlabL->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HodSlabL->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecalL_VS_HodSlabL->Smooth(1);


    TH2D * h2_TOTAL_EHecalL_VS_HeCalTankL = (TH2D*)f_in->Get("h2_TOTAL_EHecalL_VS_HeCalTankL");
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetXaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetXaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetTitleOffset(offset);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetZaxis()->SetLabelSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetZaxis()->SetTitleSize(label_size);
    h2_TOTAL_EHecalL_VS_HeCalTankL->GetYaxis()->SetLabelOffset(Loffset);
    h2_TOTAL_EHecalL_VS_HeCalTankL->Smooth(1);


    TH2D *h2_Eprompt_VS_Beta[NParticles];

    TH2D *h2_EHecal_VS_Hod2[NParticles];
    TH2D *h2_EHecal_VS_Hod3[NParticles];

    TH2D *h2_EHecalL_VS_Hod2SlabL[NParticles];
    TH2D *h2_EHecalL_VS_Hod3SlabL[NParticles];
    TH2D* h2_EHecalL_VS_HodSlabL[NParticles];

    TH2D *h2_EHecalL_VS_HeCalTankL[NParticles];

    TH2D* h2_TOTAL_EHecalL_VS_TotalHit_part[NParticles];

    for(int i=0; i<NParticles; i++) {
        h2_Eprompt_VS_Beta[i] = (TH2D*)f_in->Get("h2_Eprompt_VS_Beta_"+particles_name[i]);
        h2_Eprompt_VS_Beta[i]->GetXaxis()->SetLabelSize(label_size);
        h2_Eprompt_VS_Beta[i]->GetXaxis()->SetTitleSize(label_size);
        h2_Eprompt_VS_Beta[i]->GetYaxis()->SetLabelSize(label_size);
        h2_Eprompt_VS_Beta[i]->GetYaxis()->SetTitleSize(label_size);
        h2_Eprompt_VS_Beta[i]->GetYaxis()->SetTitleOffset(offset);
        h2_Eprompt_VS_Beta[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_Eprompt_VS_Beta[i]->GetZaxis()->SetLabelSize(label_size);
        h2_Eprompt_VS_Beta[i]->GetZaxis()->SetTitleSize(label_size);
        h2_Eprompt_VS_Beta[i]->Smooth(1);

        h2_EHecal_VS_Hod2[i] = (TH2D*)f_in->Get("h2_EHecal_VS_Hod2_"+particles_name[i]);
        h2_EHecal_VS_Hod2[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod2[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod2[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod2[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod2[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecal_VS_Hod2[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecal_VS_Hod2[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod2[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod2[i]->Smooth(1);

        
        h2_EHecal_VS_Hod3[i] = (TH2D*)f_in->Get("h2_EHecal_VS_Hod3_"+particles_name[i]);
        h2_EHecal_VS_Hod3[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod3[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod3[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod3[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod3[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecal_VS_Hod3[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecal_VS_Hod3[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecal_VS_Hod3[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecal_VS_Hod3[i]->Smooth(1);



        h2_EHecalL_VS_Hod2SlabL[i] = (TH2D*)f_in->Get("h2_EHecalL_VS_Hod2SlabL_"+particles_name[i]);
        h2_EHecalL_VS_Hod2SlabL[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecalL_VS_Hod2SlabL[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecalL_VS_Hod2SlabL[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod2SlabL[i]->Smooth(1);


        h2_EHecalL_VS_Hod3SlabL[i] = (TH2D*)f_in->Get("h2_EHecalL_VS_Hod3SlabL_"+particles_name[i]);
        h2_EHecalL_VS_Hod3SlabL[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecalL_VS_Hod3SlabL[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecalL_VS_Hod3SlabL[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_Hod3SlabL[i]->Smooth(1);
        

        h2_EHecalL_VS_HodSlabL[i] = (TH2D*)f_in->Get("h2_EHecalL_VS_HodSlabL_"+particles_name[i]);
        h2_EHecalL_VS_HodSlabL[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecalL_VS_HodSlabL[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecalL_VS_HodSlabL[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HodSlabL[i]->Smooth(1);


        h2_EHecalL_VS_HeCalTankL[i] = (TH2D*)f_in->Get("h2_EHecalL_VS_HeCalTankL_"+particles_name[i]);
        h2_EHecalL_VS_HeCalTankL[i]->GetXaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->GetXaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetTitleOffset(offset);
        h2_EHecalL_VS_HeCalTankL[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_EHecalL_VS_HeCalTankL[i]->GetZaxis()->SetLabelSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->GetZaxis()->SetTitleSize(label_size);
        h2_EHecalL_VS_HeCalTankL[i]->Smooth(1);


        h2_TOTAL_EHecalL_VS_TotalHit_part[i] = (TH2D*)f_in->Get("h2_TOTAL_EHecalL_VS_TotalHit_"+particles_name[i]);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetXaxis()->SetLabelSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetXaxis()->SetTitleSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetLabelSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetTitleSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetTitleOffset(offset);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetYaxis()->SetLabelOffset(Loffset);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetZaxis()->SetLabelSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->GetZaxis()->SetTitleSize(label_size);
        h2_TOTAL_EHecalL_VS_TotalHit_part[i]->Smooth(1);

    }

    const double th_Hod_high = 1.6;
    const double th_Hod_low = 0.8;
    const double th_HeCal = 10.0;
    const double th_HodHit = 2.0;
    const double th_HodHitTot = 3.0;

    //Plotting the histos 
    TLatex TLdbar, TLpbar;
    TLdbar.SetTextSize(2*label_size);
    TLdbar.SetTextColor(kWhite);
    TLpbar.SetTextSize(2*label_size);
    TLpbar.SetTextColor(kWhite);


    TCanvas *c_TOF = new TCanvas("c_TOF", "c_TOF", 1500, 563);
    c_TOF->Divide(2,1, 0.001, 0.001); 
    const double TopMargin = 0.05, RightMargin = 0.14, LeftMargin = 0.15, BottomMargin = 0.12;
    const double TopMargin2 = 0.05, BottomMargin2 = 0.12;

    c_TOF->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);

    h2_TOTAL_Eprompt_VS_Beta->SetTitle("");
    h2_TOTAL_Eprompt_VS_Beta->GetYaxis()->SetRangeUser(0.21, 0.5);
    h2_TOTAL_Eprompt_VS_Beta->Draw("colz");
    TLdbar.DrawLatex(65, 0.32, "#bar{d}");
    TLpbar.DrawLatex(50, 0.38, "#bar{p}");


    c_TOF->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin2);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin2);

    h2_TOTAL_EHecal_VS_Hod3->SetTitle("");
    h2_TOTAL_EHecal_VS_Hod3->GetYaxis()->SetRangeUser(0.5, 10.5);
    h2_TOTAL_EHecal_VS_Hod3->Draw("colz");
    TLdbar.DrawLatex(60, 6.5, "#bar{d}");
    TLpbar.DrawLatex(50, 4., "#bar{p}");

    TLine *l_MIP = new TLine(10, 1, 80, 1);
    l_MIP->SetLineColor(kRed);
    l_MIP->SetLineWidth(4);
    l_MIP->SetLineStyle(9);
    l_MIP->Draw("same");

    c_TOF->Update();
    c_TOF->SaveAs("/home/franz/Desktop/IFAE_2024/Proceeding/PLOTS/HeCalPromptE_TOF.eps");



    TCanvas *c_hodHIT = new TCanvas("c_hodHIT", "c_hodHIT", 1500, 563);
    c_hodHIT->Divide(2,1, 0.001, 0.001); 

    c_hodHIT->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);

    h2_EHecalL_VS_HodSlabL[Dbar]->SetTitle("");
    h2_EHecalL_VS_HodSlabL[Dbar]->GetYaxis()->SetRangeUser(1,21);
    h2_EHecalL_VS_HodSlabL[Dbar]->Draw("colz");
    TLdbar.DrawLatex(140, 18., "#bar{d}");
    

    TLine *l_pbarHitE = new TLine(50, 1.5, 50, 21);
    l_pbarHitE->SetLineColor(kRed);
    l_pbarHitE->SetLineWidth(4);
    l_pbarHitE->SetLineStyle(9);
    l_pbarHitE->Draw("same");

    TLine *l_pbarHit = new TLine(10, 8, 160, 8);
    l_pbarHit->SetLineColor(kRed);
    l_pbarHit->SetLineWidth(4);
    l_pbarHit->SetLineStyle(9);
    l_pbarHit->Draw("same");


    c_hodHIT->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin2);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin2);

    h2_EHecalL_VS_HodSlabL[Pbar]->SetTitle("");
    h2_EHecalL_VS_HodSlabL[Pbar]->GetYaxis()->SetRangeUser(1.5, 21.);
    h2_EHecalL_VS_HodSlabL[Pbar]->Draw("colz");
    TLpbar.DrawLatex(140, 18., "#bar{p}");

    l_pbarHitE->Draw("same");
    l_pbarHit->Draw("same");

    c_hodHIT->Update();
    c_hodHIT->SaveAs("/home/franz/Desktop/IFAE_2024/Proceeding/PLOTS/HeCalPromptE_TOFhit.eps");






    TCanvas *c_HeCalHIT = new TCanvas("c_HeCalHIT", "c_HeCalHIT", 1500, 563);
    c_HeCalHIT->Divide(2, 1, 0.001, 0.001); 

    c_HeCalHIT->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);

    h2_EHecalL_VS_HeCalTankL[Dbar]->SetTitle("");
    h2_EHecalL_VS_HeCalTankL[Dbar]->Draw("colz");
    TLdbar.DrawLatex(140, 9., "#bar{d}");


    c_HeCalHIT->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin2);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin2);

    h2_EHecalL_VS_HeCalTankL[Pbar]->SetTitle("");
    h2_EHecalL_VS_HeCalTankL[Pbar]->Draw("colz");
    TLpbar.DrawLatex(140, 9., "#bar{p}");


    // c_HeCalHIT->Update();
    // c_HeCalHIT->SaveAs("/home/franz/Desktop/IFAE_2024/Proceeding/PLOTS/HeCalPromptE_HeCalHit.eps");


    TCanvas *c_TotalHIT = new TCanvas("c_TotalHIT", "c_TotalHIT", 1500, 563);
    c_TotalHIT->Divide(2, 1, 0.001, 0.001); 

    c_TotalHIT->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin);

    h2_TOTAL_EHecalL_VS_TotalHit_part[Dbar]->SetTitle("");
    h2_TOTAL_EHecalL_VS_TotalHit_part[Dbar]->GetYaxis()->SetNdivisions(220);
    h2_TOTAL_EHecalL_VS_TotalHit_part[Dbar]->Draw("colz");
    TLdbar.DrawLatex(140, 27., "#bar{d}");


    TLine *l_pbarTOTHitE = new TLine(50, 1.5, 50, 31);
    l_pbarTOTHitE->SetLineColor(kRed);
    l_pbarTOTHitE->SetLineWidth(4);
    l_pbarTOTHitE->SetLineStyle(9);
    l_pbarTOTHitE->Draw("same");

    TLine *l_pbarTOTHit = new TLine(10, 12, 160, 12);
    l_pbarTOTHit->SetLineColor(kRed);
    l_pbarTOTHit->SetLineWidth(4);
    l_pbarTOTHit->SetLineStyle(9);
    l_pbarTOTHit->Draw("same");


    c_TotalHIT->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(TopMargin2);
    gPad->SetRightMargin(RightMargin);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(BottomMargin2);

    h2_TOTAL_EHecalL_VS_TotalHit_part[Pbar]->SetTitle("");
    h2_TOTAL_EHecalL_VS_TotalHit_part[Pbar]->GetYaxis()->SetNdivisions(220);
    h2_TOTAL_EHecalL_VS_TotalHit_part[Pbar]->Draw("colz");
    TLpbar.DrawLatex(140, 27., "#bar{p}");


    l_pbarTOTHitE->Draw("same");
    l_pbarTOTHit->Draw("same");

    c_TotalHIT->Update();
    gPad->Update();
    c_TotalHIT->SaveAs("/home/franz/Desktop/IFAE_2024/Proceeding/PLOTS/HeCalPromptE_TOTALhit.eps");

    

    return;

}







