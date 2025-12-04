#include <Riostream.h>
#include <map>
#include <TH1.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <TRandom3.h>
using namespace std;

void EnergySpectra()
{
    TH1 *h1 = new TH1F("h1", "Uniform Spectrum log_{10}(E_{k}); log_{10}(E_{k}); Entries", 100,1.,3.);
    //TH1 *h1 = new TH1F("h1", "Uniform Spectrum log_{10}(E_{k}); log_{10}(E_{k}); Entries", 100,1.6989700,2.6989700);
    //TH1 *h1 = new TH1F("h1", "Uniform Spectrum log_{10}(E_{k}); log_{10}(E_{k}); Entries", 100, 3.9120230,6.2146081);
    h1->Sumw2();
    TH1 *h2 = new TH1F("h2", "Kinetic energy spectrum; E_{k}; Entries", 100, 10, 1000);
    //TH1 *h2 = new TH1F("h2", "Kinetic energy spectrum; E_{k}; Entries", 100, 50, 500);
    h2->Sumw2();
    TH1 *h3 = new TH1F("h3", "Kinetic energy spectrum (from uniform); E_{k}; Entries", 100, 40, 600);
    h3->Sumw2();
    TH1 *h4 = new TH1F("h4", "Kinetic energy spectrum (from uniform 10 - 1000 ); E_{k}; Entries", 100, 10, 1000);
    h4->Sumw2();
    //enegy spectrum from MIP protons from 0.1 - 10 GeV
    TH1 *h5 = new TH1F("h5", "MIP kinetic energy spectrum (from uniform 10 - 10^{4} MeV); E_{k}; Entries", 1000, 10, 10000);
    h5->Sumw2();

    // generate uniform random numbers
    TRandom3 *r = new TRandom3();
    r->SetSeed(0);
    long int nEvents = 1000000;
    for (int i = 0; i < nEvents; i++) {
        //double rndm = r->Uniform(1.6989700,2.6989700);//log10(50) e log10(500)
        double rndm = r->Uniform(1.,3.);//log10(10) e log10(1000)
        double rndm2 = r->Uniform(0.,1.);
        //double rndm = r->Uniform(3.9120230,6.2146081); //loge(50) e loge(500)
        h1->Fill(rndm);
        //h2->Fill(exp(rndm));
        h2->Fill(pow(10,rndm));
        h3->Fill(pow(10.,rndm2)*50.);
        h4->Fill(pow(10.,2.*rndm2+1));
        h5->Fill(pow(10.,3*rndm2+1));
    }

    cout << "Integral h2: " << h2->Integral()/nEvents << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    h1->DrawNormalized("hist");
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    h2->DrawNormalized("hist");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    h3->DrawNormalized("hist");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    h4->DrawNormalized("hist");

    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    h5->DrawNormalized("hist");
}