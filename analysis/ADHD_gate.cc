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
Programma senza argomenti in ingresso. Usa dei percorsi predefiniti ai file root "rec" per antiparticelle.
Recupera il tree dai file, applica selezioni specifiche per Dbar (doppie catture) e pbar, 
per ottenere una distribuzione temporale delle hit nell'HeCal per la parte "prompt" dell'evento.
Questa distribuzione viene pesata con l'energia delle hit stesse.
Nel caso di Dbar si divide la distribuzione per due contributi in massa della cattura, Dbar e pbar.
Si confrontano tre distribuzioni normalizzate: Dbar, pbar (da Dbar) e pbar puri. 
Resistiuisce a schermo la frazione di eventi di ciascuna distribuzione al di sotto di un certo gate temporale (50 ns).
Le informazioni raccolte vengono salvate nel file di output ../sources/build/prompt_gate.root

ATTENZIONE: le distribuzioni mostrate sono inadatte a stabilire se il valore del gate e' adatto o no!!
*/

void ADHD_gate(){

  //output file 
  TFile *f_out = new TFile("../sources/build/prompt_gate.root","recreate");

  //variables name to be plotted
  const int NVar = 3;
  TString var[NVar] = {"Hod3T", "Hod2T", "HGasT"};

  

  //selection for prompt Dbar captures
  TString s_dbar[NVar] = {
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) && (900000. > Hod3T) && (CaptureM>1875) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) && (900000. > Hod2T) && (CaptureM>1875) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (900000. > HGasT) && (CaptureM>1875) )*HGasE"};
    //selection with 50 ns gate
  TString s_dbar_cut[NVar] = { 
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) && (50. > Hod3T) && (CaptureM>1875) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) && (50. > Hod2T) && (CaptureM>1875) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (50. > HGasT) && (CaptureM>1875) )*HGasE" };

  //selection for prompt pbar captures from Dbar 
  // richiedere (900000. > Hod3T) porta a considerare anche eventi rari nelle code che alterano il binnaggio dell'istogramma
  // per questo metto un taglio molto lasco a 200 ns in modo da non altarere la distribuzione e migliorare la visualizzazione
  // teoricamente e' scorretto, tieni a mente questa cosa se sorgono problemi. 
  TString s_dbar_pbar[NVar] = {
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) && (200. > Hod3T) && (CaptureM<1000) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) && (200. > Hod2T) && (CaptureM<1000) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (900000. > HGasT) && (CaptureM<1000) )*HGasE"};
    //selection with 50 ns gate
  TString s_dbar_pbar_cut[NVar] = {
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) && (50. > Hod3T) && (CaptureM<1000) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) && (50. > Hod2T) && (CaptureM<1000) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (50. > HGasT) && (CaptureM<1000) )*HGasE"};

  //selection for prompt pbar
  TString s_pbarOnly[NVar] = {
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) )*HGasE"};
    //selection with 50 ns gate
  TString s_pbarOnly_cut[NVar] = {
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod3T) && (50 > Hod3T) )*Hod3E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > Hod2T) && (50 > Hod2T) )*Hod2E",
    "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (50 > HGasT) )*HGasE"};
  
  

    // histo names
  TString h_dbarNames[NVar] = { "h_dbar_Hod3", "h_dbar_Hod2", "h_dbar_HeCal"};
  TString h_dbarNames_cut[NVar] = { "h_1_Hod3", "h_1_Hod2", "h_1_HeCal"};

  TString h_dbar_pbarNames[NVar] = { "h_pbar_Hod3", "h_pbar_Hod2", "h_pbar_HeCal"};
  TString h_dbar_pbarNames_cut[NVar] = { "h_2_Hod3", "h_2_Hod2", "h_2_HeCal"};

  TString h_pbarOnlyNames[NVar] = { "h_pbarOnly_Hod3", "h_pbarOnly_Hod2", "h_pbarOnly_HeCal"};
  TString h_pbarOnlyNames_cut[NVar] = { "h_3_Hod3", "h_3_Hod2", "h_3_HeCal"};

  TString s_1Titles[NVar]{
    "#bar{D} capture; Hod3T [ns]; Arbitrary units",
    "#bar{D} capture; Hod2T [ns]; Arbitrary units",
    "#bar{D} capture; HGasT [ns]; Arbitrary units"};

   TString s_2Titles[NVar]{
    "#bar{p} capture from #bar{D}; Hod3T [ns]; Arbitrary units",
    "#bar{p} capture from #bar{D}; Hod2T [ns]; Arbitrary units",
    "#bar{p} capture from #bar{D}; HGasT [ns]; Arbitrary units"};
  
   TString s_3Titles[NVar]{
    "#bar{p} capture; Hod3T [ns]; Arbitrary units",
    "#bar{p} capture; Hod2T [ns]; Arbitrary units",
    "#bar{p} capture; HGasT [ns]; Arbitrary units"};

  // selection for prompt SigmaMinus from Dbar (not used)
  // TString s_dbar_HeCal_sigmaMinus = "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (900000. > HGasT) && (CaptureM>1000) && (CaptureM<1200) )*HGasE";
  //   selection with 50 ns gate
  // TString s_dbar_HeCal_sigmaMinus_cut = "( (HeCal_MaxEprompt_CopyNo==HeCal_MaxEDelayed_CopyNo) && (HeCal_MaxEprompt_CopyNo==HGasCopyNo) && (CaptureT > HGasT) && (50. > HGasT) && (CaptureM>1000) && (CaptureM<1200) )*HGasE";


  // ADHD->Draw("HGasT>>h_sigmaMinus",""+s_dbar_HeCal_sigmaMinus+"","hist");
  // TH1F *h_sigmaMinus = (TH1F*)gDirectory->Get("h_sigmaMinus");

  // h_sigmaMinus->SetTitle("#Sigma^{-} capture from #bar{D}");
  // h_sigmaMinus->SetLineWidth(3); h_sigmaMinus->SetLineColor(kViolet);

  // ADHD->Draw("HGasT>>h_3",""+s_dbar_HeCal_sigmaMinus_cut+"","hist");
  // // TH1F *h_3 = (TH1F*)gDirectory->Get("h_3");

  // // printf("Fraction below 50 ns #Sigma^{-} capture from #bar{D}: %f\n", h_3->Integral()/h_sigmaMinus->Integral());

  // // h_sigmaMinus->Scale(1./h_sigmaMinus->Integral());




  // ---- Plotting ----
  TCanvas *c[NVar];
  TString cNames[NVar] = {"c_Hod3", "c_Hod2", "c_HeCal"};
  TString DNames[NVar] = {"HOD3", "HOD2", "HECAL"};
  TCanvas *c_dummy = new TCanvas("c_dummy", "c_dummy", 800, 800);


  // getting the distirbutions
  for(int i=0; i<NVar; i++){

    c_dummy->cd();


    // ---- DBAR projectile ----
    TFile *File = TFile::Open("../sources/build/Dbar_10e7_rec_corrected.root","read");
   TTree* ADHD = (TTree*) File->Get("RHits_bar");


    //getting the histogram and setting the style
    ADHD->Draw(""+var[i]+">>"+h_dbarNames[i]+"",""+s_dbar[i]+"","hist");

    TH1F *h_1 = (TH1F*)gDirectory->Get(h_dbarNames[i].Data());
    h_1->SetTitle(s_1Titles[i]);
    h_1->SetLineWidth(3); h_1->SetLineColor(kRed);

    ADHD->Draw(""+var[i]+">>"+h_dbar_pbarNames[i]+"",""+s_dbar_pbar[i]+"","hist");
  
    TH1F *h_2 = (TH1F*)gDirectory->Get(h_dbar_pbarNames[i].Data());
    h_2->SetTitle(s_2Titles[i]);
    h_2->SetLineWidth(3); h_2->SetLineColor(kBlue);



    ADHD->Draw(""+var[i]+">>"+h_dbarNames_cut[i]+"",""+s_dbar_cut[i]+"","hist");

    TH1F *h_1_cut = (TH1F*)gDirectory->Get(h_dbarNames_cut[i].Data());

    ADHD->Draw(""+var[i]+">>"+h_dbar_pbarNames_cut[i]+"",""+s_dbar_pbar_cut[i]+"","hist");
    TH1F *h_2_cut = (TH1F*)gDirectory->Get(h_dbar_pbarNames_cut[i].Data());


    //fraction of the selected evebnts
    printf("--------- %s ---------\n", DNames[i].Data());
    printf("Fraction below 50 ns #bar{D} capture: %f\n", h_1_cut->Integral()/h_1->Integral());
    printf("Fraction below 50 ns #bar{p} capture from #bar{D}: %f\n", h_2_cut->Integral()/h_2->Integral());



    //normalization
    h_1->Scale(1./h_1->Integral());
    h_2->Scale(1./h_2->Integral());

    //Writing the histograms in the output file
    f_out->cd();

    h_1->Write(h_dbarNames[i].Data());
    h_2->Write(h_dbar_pbarNames[i].Data());
    // h_sigmaMinus->Write("h_sigmaMinus");



    // ---- PBAR projectile ----
    TFile *File2 = TFile::Open("../sources/build/pbar_10e7_rec.root","read");
    TTree* t = (TTree*) File2->Get("RHits_bar");

    t->Draw(""+var[i]+">>"+h_pbarOnlyNames[i]+"",""+s_pbarOnly[i]+"","hist");
    TH1F *h_3 = (TH1F*)gDirectory->Get(h_pbarOnlyNames[i].Data());
    h_3->SetTitle(s_3Titles[i]);
    h_3->SetLineWidth(3); h_3->SetLineColor(kGreen);

      //50 second gate histos
    t->Draw(""+var[i]+">>"+h_pbarOnlyNames_cut[i]+"",""+s_pbarOnly_cut[i]+"","hist");
    TH1F *h_3_cut = (TH1F*)gDirectory->Get(h_pbarOnlyNames_cut[i].Data());

      //fraction of the selected evebnts
    printf("Fraction below 50 ns #bar{p} capture: %f\n", h_3_cut->Integral()/h_3->Integral());


      //normalization
    h_3->Scale(1./h_3->Integral());

    //Writing the histograms in the output file
    f_out->cd();
    h_3->Write(h_pbarOnlyNames[i].Data());



    c[i] = new TCanvas(Form("c%d",i), cNames[i], 800, 800);
    c[i]->cd();
    h_3->DrawCopy("hist");
    h_2->DrawCopy("hist same");
    h_1->DrawCopy("hist same");
    c[i]->BuildLegend();


  }

  f_out->Write();

  return;
}
