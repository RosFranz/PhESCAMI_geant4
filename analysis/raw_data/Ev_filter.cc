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
Programma con percorso al file da selezionare come parametro di ingresso.
Legge il tree, seleziona solo gli eventi con catture in He!
DA USARE per file CON ANTIPARTICELLE.
*/

void Ev_filter(TString fname,TString cuts = "CaptureT>0"){
  
  TFile* F1 = new TFile(Form("%s.root",fname.Data()), "OPEN");
  if(!F1) {
    printf("\n !!!!!!!!!!! \n Do not put the '.root' in the name file!\n !!!!!!!!!!! \n");
    return;
  }
  TTree* tree = (TTree*) F1->Get("Hits");

  cout << "Output file = ../filtred_data/" << fname.Data() << "_filtred.root" << endl;
  cout << "applied cuts = " << cuts.Data() << endl;
  
  TFile* F2 = new TFile(Form("../filtred_data/%s_filtred.root",fname.Data()),"recreate");
  TTree* tree2 = (TTree*) tree->CopyTree(cuts.Data());

  // Saving the new tree
  F2->Write();
  F2->Close();
  tree->Delete();
  F1->Close();

  return;
}
