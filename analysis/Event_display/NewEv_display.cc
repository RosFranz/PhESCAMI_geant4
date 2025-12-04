#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <Riostream.h>
#include <map>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <TGraph2D.h>



using namespace std;

void NewEv_display(TString inputFile, TString outputDir)
{
    // ---- 1st parameter is the INPUT FILE ---- 2nd parameter the OUTPUT DIR ----- //

    //Programma per avere una ricostruzione 3D degli eventi
    //bisogna specificare l'ide degli eventi che si vogliono ricostruire


    TFile *f_input = new TFile(inputFile,"OPEN");

    gStyle->SetOptStat(0);

    if (!f_input){
        printf("File not found");
        return;
    } 


    // ----------------- GETTING THE TREE ----------------- //
    bool isBar = false;
    TTree *t;
    //Check if the file contains anti-particle
    if(inputFile.Contains("bar")) {
        t = (TTree*)f_input->Get("Hits_bar");
        printf("\nThis is an anti-particle file\n");
        isBar = true;
        }

    else {
        t = (TTree*)f_input->Get("RHits");
        printf("\nThis is an ordinary-particle file\n");
    }
    
    const Long64_t entries = t->GetEntries();
    printf("Number of entries: %lld\n\n", entries);

    // ---------- CONSTANT USED FOR THE EVENT DISPLAY ---------- //
    const int NEvents = 3;//5;
    Long64_t IDs[NEvents] = {48, 49, 77};

    // MC TRUTH

    double MCEnergy, MCMomemtum;
    double MCMomX, MCMomY, MCMomZ;
    double MCPosX, MCPosY, MCPosZ;
    double MCMass, MCCharge;
    int EventID;

        // setting the branches
    t->SetBranchAddress("MCEnergy",&MCEnergy);
    t->SetBranchAddress("MCMomentum",&MCMomemtum);
    t->SetBranchAddress("MCMass",&MCMass);
    t->SetBranchAddress("MCCharge",&MCCharge);
        // momentum components
    t->SetBranchAddress("MCMomX",&MCMomX);
    t->SetBranchAddress("MCMomY",&MCMomY);
    t->SetBranchAddress("MCMomZ",&MCMomZ);
        // position components    const int NEvents = int(MaxEventNumberID-MinEventNumberID)+1;
    t->SetBranchAddress("MCPosX",&MCPosX);
    t->SetBranchAddress("MCPosY",&MCPosY);
    t->SetBranchAddress("MCPosZ",&MCPosZ);
    t->SetBranchAddress("EventID",&EventID);

    // PRIMARY PARTICLE POSITION
        // array of TH2
    TH2D* h_primaryXY[NEvents];
    TH2D* h_primaryYZ[NEvents];
    TH2D* h_primaryXZ[NEvents];

    TGraph2D* g_primaryXYZ[NEvents];


    // POSITION, ENERGY AND TIME OF THE CAPTURE 
        //IN HELIUM
    vector<double> *CaptureX=0, *CaptureY=0, *CaptureZ=0, *CaptureT=0;
    TBranch *b_CaptureX=0, *b_CaptureY=0, *b_CaptureZ=0, *b_CaptureT=0;
    vector<double> *CaptureE=0, *CaptureM=0, *CaptureF=0;
    TBranch *b_CaptureE=0, *b_CaptureM=0, *b_CaptureF=0;

    TH2D* h_captureXY[NEvents];
    TH2D* h_captureYZ[NEvents];
    TH2D* h_captureXZ[NEvents];
    TH2D* h_captureET[NEvents];

    TGraph2D* g_captureXYZ[NEvents];

        //IN THE HODOSCOPES
    vector<double> *HodCaptureX=0, *HodCaptureY=0, *HodCaptureZ=0, *HodCaptureT=0;
    TBranch *b_HodCaptureX=0, *b_HodCaptureY=0, *b_HodCaptureZ=0, *b_HodCaptureT=0;
    vector<double> *HodCaptureE=0, *HodCaptureM=0, *HodCaptureF=0;
    TBranch *b_HodCaptureE=0, *b_HodCaptureM=0, *b_HodCaptureF=0;

    vector<int> *HodCaptureVolume=0;
    TBranch *b_HodCaptureVolume=0;        


    if(isBar){
        t->SetBranchAddress("CaptureX",&CaptureX,&b_CaptureX);
        t->SetBranchAddress("CaptureY",&CaptureY,&b_CaptureY);
        t->SetBranchAddress("CaptureZ",&CaptureZ,&b_CaptureZ);
        t->SetBranchAddress("CaptureT",&CaptureT,&b_CaptureT);

        
        t->SetBranchAddress("CaptureE",&CaptureE,&b_CaptureE);
        t->SetBranchAddress("CaptureM",&CaptureM,&b_CaptureM);
        t->SetBranchAddress("CaptureF",&CaptureF,&b_CaptureF);

        t->SetBranchAddress("HodCaptureX",&HodCaptureX,&b_HodCaptureX);
        t->SetBranchAddress("HodCaptureY",&HodCaptureY,&b_HodCaptureY);
        t->SetBranchAddress("HodCaptureZ",&HodCaptureZ,&b_HodCaptureZ);
        t->SetBranchAddress("HodCaptureT",&HodCaptureT,&b_HodCaptureT);

        
        t->SetBranchAddress("HodCaptureE",&HodCaptureE,&b_HodCaptureE);
        t->SetBranchAddress("HodCaptureM",&HodCaptureM,&b_HodCaptureM);
        t->SetBranchAddress("HodCaptureF",&HodCaptureF,&b_HodCaptureF);

        t->SetBranchAddress("HodCaptureVolume",&HodCaptureVolume,&b_HodCaptureVolume);
    }  


    
    TH2D* h_Hod1captureXY[NEvents];
    TH2D* h_Hod1captureYZ[NEvents];
    TH2D* h_Hod1captureXZ[NEvents];
    TH2D* h_Hod1captureET[NEvents];

    TGraph2D* g_Hod1captureXYZ[NEvents];

    TH2D* h_Hod2captureXY[NEvents];
    TH2D* h_Hod2captureYZ[NEvents];
    TH2D* h_Hod2captureXZ[NEvents];
    TH2D* h_Hod2captureET[NEvents];

    TGraph2D* g_Hod2captureXYZ[NEvents];
   

    // INNER SCINTILLATOR
    vector<double> *Hod2X=0, *Hod2Y=0, *Hod2Z=0, *Hod2T=0, *Hod2E=0;
    TBranch *b_Hod2X=0, *b_Hod2Y=0, *b_Hod2Z=0, *b_Hod2T=0, *b_Hod2E=0;
    t->SetBranchAddress("Hod2X",&Hod2X,&b_Hod2X);
    t->SetBranchAddress("Hod2Y",&Hod2Y,&b_Hod2Y);
    t->SetBranchAddress("Hod2Z",&Hod2Z,&b_Hod2Z);
    t->SetBranchAddress("Hod2T",&Hod2T,&b_Hod2T);
    t->SetBranchAddress("Hod2E",&Hod2E,&b_Hod2E);

        // spatial array of TH2 before the capture
    TH2D* h_InnScintXY[NEvents];
    TH2D* h_InnScintYZ[NEvents];
    TH2D* h_InnScintXZ[NEvents];
    TH2D* h_InnScintET[NEvents];

    TGraph2D* g_InnScintXYZ[NEvents];

        // spatial array of TH2 after the capture
    TH2D* h_InnScintXY_after[NEvents];
    TH2D* h_InnScintYZ_after[NEvents];
    TH2D* h_InnScintXZ_after[NEvents];
    TH2D* h_InnScintET_after[NEvents];

    TGraph2D* g_InnScintXYZ_after[NEvents];


    // EXTERNAL SCINTILLATOR
    vector<double> *Hod3X=0, *Hod3Y=0, *Hod3Z=0, *Hod3T=0, *Hod3E=0;
    TBranch *b_Hod3X=0, *b_Hod3Y=0, *b_Hod3Z=0, *b_Hod3T=0, *b_Hod3E=0;
    t->SetBranchAddress("Hod3X",&Hod3X,&b_Hod3X);
    t->SetBranchAddress("Hod3Y",&Hod3Y,&b_Hod3Y);
    t->SetBranchAddress("Hod3Z",&Hod3Z,&b_Hod3Z);
    t->SetBranchAddress("Hod3T",&Hod3T,&b_Hod3T);
    t->SetBranchAddress("Hod3E",&Hod3E,&b_Hod3E);

        // spatial array of TH2 before the capture
    TH2D* h_ExtScintXY[NEvents];
    TH2D* h_ExtScintYZ[NEvents];
    TH2D* h_ExtScintXZ[NEvents];
    TH2D* h_ExtScintET[NEvents];

    TGraph2D* g_ExtScintXYZ[NEvents];

        // spatial array of TH2 after the capture
    TH2D* h_ExtScintXY_after[NEvents];
    TH2D* h_ExtScintYZ_after[NEvents];
    TH2D* h_ExtScintXZ_after[NEvents];
    TH2D* h_ExtScintET_after[NEvents];

    TGraph2D* g_ExtScintXYZ_after[NEvents];

    

    // SPATIAL AND TIME RESOLUTION
        //spatial dimension
    double maxDim = 1620; // (1500mm) 1.5m max dimension of the detector
        //spatial resolution
    double SCDist = 5; // (5mm) 0.5cm Plastic hit resolution
    int nBinsSCDist = int(2.*maxDim/SCDist);

        //temporal dimension
    double maxTime = 2000; //[micro s] 900 us fixed delay for Z=-1 massive particle capture
        //temporal resolution
    double SCTres = 0.1; // (0.1ns) 100 ps plastic hit resolution
    int nBinsSCTime = int(maxTime/SCTres);
     
    
    // HE SCINTILLATOR
    vector<double> *HGasX=0, *HGasY=0, *HGasZ=0, *HGasT=0, *HGasE=0;
    TBranch *b_HGasX=0, *b_HGasY=0, *b_HGasZ=0, *b_HGasT=0, *b_HGasE=0;
    t->SetBranchAddress("HGasX",&HGasX,&b_HGasX);
    t->SetBranchAddress("HGasY",&HGasY,&b_HGasY);
    t->SetBranchAddress("HGasZ",&HGasZ,&b_HGasZ);
    t->SetBranchAddress("HGasT",&HGasT,&b_HGasT);
    t->SetBranchAddress("HGasE",&HGasE,&b_HGasE);
    
    vector<int> *HGasCopyNo=0;
    TBranch *b_HGasCopyNo=0;
    t->SetBranchAddress("HGasCopyNo",&HGasCopyNo,&b_HGasCopyNo);

        // spatial array of TH2 before the capture
    TH2D* h_HeXY[NEvents];
    TH2D* h_HeYZ[NEvents];
    TH2D* h_HeXZ[NEvents];
    TH2D* h_HeET[NEvents];

    TGraph2D* g_HeXYZ[NEvents];

        // spatial array of TH2 after the capture
    TH2D* h_HeXY_after[NEvents];
    TH2D* h_HeYZ_after[NEvents];
    TH2D* h_HeXZ_after[NEvents];
    TH2D* h_HeET_after[NEvents];

    TGraph2D* g_HeXYZ_after[NEvents];
    

    // SPATIAL AND TIME RESOLUTION                                    
    double HeDist = 20; // (20mm) 2cm He hit resolution
    int nBinsHeDist = int(2.*maxDim/HeDist);

        //temporal resolution
    double HeTres = 0.2; // 200 ps He hit resolution
    int nBinsHeTime = int(maxTime/HeTres);
    //cout << "nBinsHeTime = " << nBinsHeTime << endl;
    
    
    

    // ------------------- READING THE TREE ------------------------ //

    int before = 0; int after = 0;
    int ext = 0; int inn = 0;
    bool flag1[NEvents]; bool flag2[NEvents];
    for(int i = 0; i < NEvents; i++){
        flag1[i] = false;
        flag2[i] = false;
    }


    int counter = 0;    
   
    for(Int_t entry=0; entry<entries; entry++){

        std::stringstream strID;
        TString ID;
        t->GetEntry(entry);

        if(EventID == IDs[counter]){


            printf("EVENT: %d  \t tree entry: %d \n",EventID, entry);

                //setting the histo names
            strID << EventID;
            ID = strID.str();

            // EXTERNAL SCINTILLATOR
                // before the capture
            h_ExtScintXY[counter] = new TH2D("h_ExtScintXY_"+ID,"Ext. scint.(XY) ID ="+ID+";X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintYZ[counter] = new TH2D("h_ExtScintYZ_"+ID,"Ext. scint.(YZ) ID ="+ID+";Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintXZ[counter] = new TH2D("h_ExtScintXZ_"+ID,"Ext. scint.(XZ) ID ="+ID+";X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintET[counter] = new TH2D("h_ExtScintET_"+ID,"Ext. scint.(ET) ID ="+ID+";E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // after the capture
            if(isBar){
                h_ExtScintXY_after[counter] = new TH2D("h_ExtScintXY_after_"+ID,"Ext. scint.(XY) (T>100 #mu s);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_ExtScintYZ_after[counter] = new TH2D("h_ExtScintYZ_after_"+ID,"Ext. scint.(YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_ExtScintXZ_after[counter] = new TH2D("h_ExtScintXZ_after_"+ID,"Ext. scint.(XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_ExtScintET_after[counter] = new TH2D("h_ExtScintET_after_"+ID,"Ext. scint.(ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);


                // printf("Hod3X->size() = %zu\n", Hod3X->size());
                printf("CaptureT->size() = %zu\n", CaptureT->size());
            
                for(int k=0; k<Hod3X->size(); k++) {
                    if(Hod3T->at(k) < CaptureT->at(0)){
                        h_ExtScintXY[counter]->Fill(Hod3X->at(k), Hod3Y->at(k));
                        h_ExtScintYZ[counter]->Fill(Hod3Y->at(k), Hod3Z->at(k));
                        h_ExtScintXZ[counter]->Fill(Hod3X->at(k), Hod3Z->at(k));
                        h_ExtScintET[counter]->Fill(Hod3E->at(k), Hod3T->at(k));
                        before++;
                    }
                    else{
                        h_ExtScintXY_after[counter]->Fill(Hod3X->at(k), Hod3Y->at(k));
                        h_ExtScintYZ_after[counter]->Fill(Hod3Y->at(k), Hod3Z->at(k));
                        h_ExtScintXZ_after[counter]->Fill(Hod3X->at(k), Hod3Z->at(k));
                        h_ExtScintET_after[counter]->Fill(Hod3E->at(k), Hod3T->at(k));
                        after++;
                    }
                }

                double H3X[before+2]; H3X[0]=-maxDim; H3X[before+1]=maxDim;
                double H3Y[before+2]; H3Y[0]=-maxDim; H3Y[before+1]=maxDim;
                double H3Z[before+2]; H3Z[0]=-maxDim; H3Z[before+1]=maxDim;

                double H3X_after[after+2]; H3X_after[0]=-maxDim; H3X_after[after+1]=maxDim;
                double H3Y_after[after+2]; H3Y_after[0]=-maxDim; H3Y_after[after+1]=maxDim;
                double H3Z_after[after+2]; H3Z_after[0]=-maxDim; H3Z_after[after+1]=maxDim;

                before = 0; after = 0;

                for(int k=0; k<Hod3X->size(); k++){
                    if(Hod3T->at(k) < CaptureT->at(0)){
                        H3X[before+1]=Hod3X->at(k);
                        H3Y[before+1]=Hod3Y->at(k);
                        H3Z[before+1]=Hod3Z->at(k);
                        before++;
                    }
                    else{
                        H3X_after[after+1]=Hod3X->at(k);
                        H3Y_after[after+1]=Hod3Y->at(k);
                        H3Z_after[after+1]=Hod3Z->at(k);
                        after++;   
                    }
                }
                
                g_ExtScintXYZ[counter] = new TGraph2D("Ext_"+ID, "Ext_"+ID, before+2, H3X, H3Z, H3Y);
                g_ExtScintXYZ[counter]->GetXaxis()->SetTitle("X [mm]"); g_ExtScintXYZ[counter]->GetYaxis()->SetTitle("Z [mm]"); g_ExtScintXYZ[counter]->GetZaxis()->SetTitle("Y [mm]");
                if(isBar) g_ExtScintXYZ_after[counter] = new TGraph2D("Ext_after_"+ID, "Ext_after", after+2, H3X_after, H3Z_after, H3Y_after);

                before = 0; after = 0;
            }
            
            // ordinary particles
            else{
                
                for(int k=0; k<Hod3X->size(); k++) {
                    h_ExtScintXY[counter]->Fill(Hod3X->at(k), Hod3Y->at(k));
                    h_ExtScintYZ[counter]->Fill(Hod3Y->at(k), Hod3Z->at(k));
                    h_ExtScintXZ[counter]->Fill(Hod3X->at(k), Hod3Z->at(k));
                    h_ExtScintET[counter]->Fill(Hod3E->at(k), Hod3T->at(k));
                    before++;
                }

                double H3X[before+2]; H3X[0]=-maxDim; H3X[before+1]=maxDim;
                double H3Y[before+2]; H3Y[0]=-maxDim; H3Y[before+1]=maxDim;
                double H3Z[before+2]; H3Z[0]=-maxDim; H3Z[before+1]=maxDim;

                before = 0;

                for(int k=0; k<Hod3X->size(); k++){

                    H3X[before+1]=Hod3X->at(k);
                    H3Y[before+1]=Hod3Y->at(k);
                    H3Z[before+1]=Hod3Z->at(k);
                    before++;
                }
                
                g_ExtScintXYZ[counter] = new TGraph2D("Ext_"+ID, "Ext_"+ID, before+2, H3X, H3Z, H3Y);
                g_ExtScintXYZ[counter]->GetXaxis()->SetTitle("X [mm]"); g_ExtScintXYZ[counter]->GetYaxis()->SetTitle("Z [mm]"); g_ExtScintXYZ[counter]->GetZaxis()->SetTitle("Y [mm]");

                before = 0;
            }

            

            // INNER SCINTILLATOR
                // before the capture
            h_InnScintXY[counter] = new TH2D("h_InnScintXY_"+ID,"Inn. scint(XY);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintYZ[counter] = new TH2D("h_InnScintYZ_"+ID,"Inn. scint(YZ);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintXZ[counter] = new TH2D("h_InnScintXZ_"+ID,"Inn. scint(XZ);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintET[counter] = new TH2D("h_InnScintET_"+ID,"Inn. scint(ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // after the capture
            if(isBar){
                h_InnScintXY_after[counter] = new TH2D("h_InnScintXY_after_"+ID,"Inn. scint(XY) (T>100 #mu s);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_InnScintYZ_after[counter] = new TH2D("h_InnScintYZ_after_"+ID,"Inn. scint(YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_InnScintXZ_after[counter] = new TH2D("h_InnScintXZ_after_"+ID,"Inn. scint(XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_InnScintET_after[counter] = new TH2D("h_InnScintET_after_"+ID,"Inn. scint(ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                printf("Hod2X->size() = %zu\n", Hod2X->size());

                    // FILLING
                for(int k=0; k<Hod2X->size(); k++){
                    if(Hod2T->at(k) < CaptureT->at(0)){
                        h_InnScintXY[counter]->Fill(Hod2X->at(k), Hod2Y->at(k));
                        h_InnScintYZ[counter]->Fill(Hod2Y->at(k), Hod2Z->at(k));
                        h_InnScintXZ[counter]->Fill(Hod2X->at(k), Hod2Z->at(k));
                        h_InnScintET[counter]->Fill(Hod2E->at(k), Hod2T->at(k));
                        before++;
                    }
                    else{
                        h_InnScintXY_after[counter]->Fill(Hod2X->at(k), Hod2Y->at(k));
                        h_InnScintYZ_after[counter]->Fill(Hod2Y->at(k), Hod2Z->at(k));
                        h_InnScintXZ_after[counter]->Fill(Hod2X->at(k), Hod2Z->at(k));
                        h_InnScintET_after[counter]->Fill(Hod2E->at(k), Hod2T->at(k));
                        after++;
                    }
                }

                double H2X[before];// H2X[0]=-maxDim; H2X[before+1]=maxDim;
                double H2Y[before];// H2Y[0]=-maxDim; H2Y[before+1]=maxDim;
                double H2Z[before];// H2Z[0]=-maxDim; H2Z[before+1]=maxDim;

                double H2X_after[after];// H2X_after[0]=-maxDim; H2X_after[after+1]=maxDim;
                double H2Y_after[after];// H2Y_after[0]=-maxDim; H2Y_after[after+1]=maxDim;
                double H2Z_after[after];// H2Z_after[0]=-maxDim; H2Z_after[after+1]=maxDim;

                before = 0; after = 0;

                for(int k=0; k<Hod2X->size(); k++){
                    if(Hod2T->at(k) < CaptureT->at(0)){
                        H2X[before]=Hod2X->at(k);
                        H2Y[before]=Hod2Y->at(k);
                        H2Z[before]=Hod2Z->at(k);
                        before++;
                    }
                    else{
                        H2X_after[after]=Hod2X->at(k);
                        H2Y_after[after]=Hod2Y->at(k);
                        H2Z_after[after]=Hod2Z->at(k);
                        after++;
                    }
                }

                g_InnScintXYZ[counter] = new TGraph2D("Inn_"+ID, "Inn", before, H2X, H2Z, H2Y);
                if(isBar) g_InnScintXYZ_after[counter] = new TGraph2D("Inn_after_"+ID, "Inn_after", after, H2X_after, H2Z_after, H2Y_after);

                before = 0; after = 0;
            }

            //ordinary particle
            else{
                for(int k=0; k<Hod2X->size(); k++){
                    h_InnScintXY[counter]->Fill(Hod2X->at(k), Hod2Y->at(k));
                    h_InnScintYZ[counter]->Fill(Hod2Y->at(k), Hod2Z->at(k));
                    h_InnScintXZ[counter]->Fill(Hod2X->at(k), Hod2Z->at(k));
                    h_InnScintET[counter]->Fill(Hod2E->at(k), Hod2T->at(k));
                    before++;
                }

                double H2X[before];// H2X[0]=-maxDim; H2X[before+1]=maxDim;
                double H2Y[before];// H2Y[0]=-maxDim; H2Y[before+1]=maxDim;
                double H2Z[before];// H2Z[0]=-maxDim; H2Z[before+1]=maxDim;

                before = 0;

                for(int k=0; k<Hod2X->size(); k++){
                    H2X[before]=Hod2X->at(k);
                    H2Y[before]=Hod2Y->at(k);
                    H2Z[before]=Hod2Z->at(k);
                    before++;
                }

                g_InnScintXYZ[counter] = new TGraph2D("Inn_"+ID, "Inn", before, H2X, H2Z, H2Y);

                before = 0;

            }
            
            

            // HELIUM SCINTILLATOR
                // before the capture
            h_HeXY[counter] = new TH2D("h_HeXY_"+ID,"He (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeYZ[counter] = new TH2D("h_HeYZ_"+ID,"He (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeXZ[counter] = new TH2D("h_HeXZ_"+ID,"He (XZ);Z [mm];X [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeET[counter] = new TH2D("h_HeET_"+ID,"He (ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // He after the capture
            if(isBar){
                h_HeXY_after[counter] = new TH2D("h_HeXY_after_"+ID,"He (XY) (T>100 #mu s);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_HeYZ_after[counter] = new TH2D("h_HeYZ_after_"+ID,"He (YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_HeXZ_after[counter] = new TH2D("h_HeXZ_after_"+ID,"He (XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_HeET_after[counter] = new TH2D("h_HeET_after_"+ID,"He (ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                printf("HGasX->size() = %zu\n", HGasX->size());

                    // FILLING
                for(int k=0; k<HGasX->size(); k++){
                    if(HGasT->at(k) < CaptureT->at(0)){
                        h_HeXY[counter]->Fill(HGasX->at(k),HGasY->at(k));
                        h_HeYZ[counter]->Fill(HGasY->at(k),HGasZ->at(k));
                        h_HeXZ[counter]->Fill(HGasX->at(k),HGasZ->at(k));
                        h_HeET[counter]->Fill(HGasE->at(k),HGasT->at(k));
                        before++;
                    }
                    else {
                        h_HeXY_after[counter]->Fill(HGasX->at(k),HGasY->at(k));
                        h_HeYZ_after[counter]->Fill(HGasY->at(k),HGasZ->at(k));
                        h_HeXZ_after[counter]->Fill(HGasX->at(k),HGasZ->at(k));
                        h_HeET_after[counter]->Fill(HGasE->at(k),HGasT->at(k));
                        after++;
                    }
                }

                double GasX[before];// GasX[0]=-maxDim; GasX[before+1]=maxDim;
                double GasY[before];// GasY[0]=-maxDim; GasY[before+1]=maxDim;
                double GasZ[before];// GasZ[0]=-maxDim; GasZ[before+1]=maxDim;

                double GasX_after[after];// GasX_after[0]=-maxDim; GasX_after[after+1]=maxDim;
                double GasY_after[after];// GasY_after[0]=-maxDim; GasY_after[after+1]=maxDim;
                double GasZ_after[after];// GasZ_after[0]=-maxDim; GasZ_after[after+1]=maxDim;

                before = 0; after = 0;

                for(int k=0; k<HGasX->size(); k++){
                    if(HGasT->at(k) < CaptureT->at(0)){
                        GasX[before]=HGasX->at(k);
                        GasY[before]=HGasY->at(k);
                        GasZ[before]=HGasZ->at(k);
                        before++;
                    }
                    else {
                        GasX_after[after]=HGasX->at(k);
                        GasY_after[after]=HGasY->at(k);
                        GasZ_after[after]=HGasZ->at(k);
                        after++;
                    }
                }

                g_HeXYZ[counter] = new TGraph2D("He_"+ID, "He", before, GasX, GasZ, GasY);
                if(isBar) g_HeXYZ_after[counter] = new TGraph2D("He_after_"+ID, "He_after", after, GasX_after, GasZ_after, GasY_after);

                before = 0; after = 0;
            }


            //ordinary particles
            else{

                    // FILLING
                for(int k=0; k<HGasX->size(); k++){
                    h_HeXY[counter]->Fill(HGasX->at(k),HGasY->at(k));
                    h_HeYZ[counter]->Fill(HGasY->at(k),HGasZ->at(k));
                    h_HeXZ[counter]->Fill(HGasX->at(k),HGasZ->at(k));
                    h_HeET[counter]->Fill(HGasE->at(k),HGasT->at(k));
                    before++;
                }

                double GasX[before];// GasX[0]=-maxDim; GasX[before+1]=maxDim;
                double GasY[before];// GasY[0]=-maxDim; GasY[before+1]=maxDim;
                double GasZ[before];// GasZ[0]=-maxDim; GasZ[before+1]=maxDim;

                before = 0;

                for(int k=0; k<HGasX->size(); k++){
                        GasX[before]=HGasX->at(k);
                        GasY[before]=HGasY->at(k);
                        GasZ[before]=HGasZ->at(k);
                        before++;
                }

                g_HeXYZ[counter] = new TGraph2D("He_"+ID, "He", before, GasX, GasZ, GasY);

                before = 0;
            }

            


            // HELIUM CAPTURE POSITION, ENERGY AND TIME 
            if(isBar){

                h_captureXY[counter] = new TH2D("h_captureXY_"+ID,"Capture (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_captureYZ[counter] = new TH2D("h_captureYZ_"+ID,"Capture (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_captureXZ[counter] = new TH2D("h_captureXZ_"+ID,"Capture (XZ);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
                h_captureET[counter] = new TH2D("h_captureET_"+ID,"Capture (ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                double HeCX[CaptureX->size()];// HeCX[0]=-maxDim; HeCX[CaptureX->size()+1]=maxDim;
                double HeCY[CaptureX->size()];// HeCY[0]=-maxDim; HeCY[CaptureX->size()+1]=maxDim;
                double HeCZ[CaptureX->size()];// HeCZ[0]=-maxDim; HeCZ[CaptureX->size()+1]=maxDim;

                    // FILLING
                printf("CaptureX->size() = %zu\n", CaptureX->size());
                
                for(int k=0; k<CaptureX->size(); k++){
                    h_captureXY[counter]->Fill(CaptureX->at(k),CaptureY->at(k));
                    h_captureYZ[counter]->Fill(CaptureY->at(k),CaptureZ->at(k));
                    h_captureXZ[counter]->Fill(CaptureX->at(k),CaptureZ->at(k));
                    h_captureET[counter]->Fill(CaptureE->at(k),CaptureT->at(k));
                    HeCX[k]=CaptureX->at(k);
                    HeCY[k]=CaptureY->at(k);
                    HeCZ[k]=CaptureZ->at(k);
                    printf(" Helium capture --> X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",CaptureX->at(k),CaptureY->at(k),CaptureZ->at(k),CaptureE->at(k),CaptureT->at(k),CaptureM->at(k),CaptureF->at(k));
                    if(k==1){
                        printf("!!!! --- SECOND CAPTURE IDENTIFIED --- !!!!\n");
                        printf("X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",CaptureX->at(k),CaptureY->at(k),CaptureZ->at(k),CaptureE->at(k),CaptureT->at(k),CaptureM->at(k),CaptureF->at(k));
                        //printf("Primary particle charge: %f ; mass: %f\n", MCCharge, MCMass);
                    }
                }
                
                

                g_captureXYZ[counter] = new TGraph2D("HeCapture_"+ID, "HeCapture", CaptureX->size(), HeCX, HeCZ, HeCY);



                // HODOSCOPE CAPTURES
                h_Hod1captureXY[counter] = new TH2D("h_InnScintCaptureXY_"+ID,"Inn. scint. capture ;X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod1captureYZ[counter] = new TH2D("h_InnScintCaptureYZ_"+ID,"Inn. scint. capture;Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod1captureXZ[counter] = new TH2D("h_InnScintCaptureXZ_"+ID,"Inn. scint. capture;X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod1captureET[counter] = new TH2D("h_InnScintCaptureET_"+ID,"Inn. scint. capture;E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                h_Hod2captureXY[counter] = new TH2D("h_ExtScintCaptureXY_"+ID,"Ext. scint. capture ;X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod2captureYZ[counter] = new TH2D("h_ExtScintCaptureYZ_"+ID,"Ext. scint. capture ;Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod2captureXZ[counter] = new TH2D("h_ExtScintCaptureXZ_"+ID,"Ext. scint. capture ;X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
                h_Hod2captureET[counter] = new TH2D("h_ExtScintCaptureET_"+ID,"Ext. scint. capture ;E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);


                printf("HodCaptureX->size() = %zu\n", HodCaptureX->size());



                for(int k=0; k<HodCaptureX->size(); k++){
                    if(HodCaptureT->at(k)>0 && HodCaptureVolume->at(k)==2){
                        h_Hod2captureXY[counter]->Fill(HodCaptureX->at(k),HodCaptureY->at(k));
                        h_Hod2captureYZ[counter]->Fill(HodCaptureY->at(k),HodCaptureZ->at(k));
                        h_Hod2captureXZ[counter]->Fill(HodCaptureX->at(k),HodCaptureZ->at(k));
                        h_Hod2captureET[counter]->Fill(HodCaptureE->at(k),HodCaptureT->at(k));
                        printf(" EXTERNAL HOD CAPTURE X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",HodCaptureX->at(k),HodCaptureY->at(k),HodCaptureZ->at(k),HodCaptureE->at(k),HodCaptureT->at(k),HodCaptureM->at(k),HodCaptureF->at(k));
                        ext++;
                    }
                    if(HodCaptureT->at(k)>0 && HodCaptureVolume->at(k)==1){
                        h_Hod1captureXY[counter]->Fill(HodCaptureX->at(k),HodCaptureY->at(k));
                        h_Hod1captureYZ[counter]->Fill(HodCaptureY->at(k),HodCaptureZ->at(k));
                        h_Hod1captureXZ[counter]->Fill(HodCaptureX->at(k),HodCaptureZ->at(k));
                        h_Hod1captureET[counter]->Fill(HodCaptureE->at(k),HodCaptureT->at(k));
                        printf(" INTERNAL HOD CAPTURE X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f ; k: %i\n",HodCaptureX->at(k),HodCaptureY->at(k),HodCaptureZ->at(k),HodCaptureE->at(k),HodCaptureT->at(k),HodCaptureM->at(k),HodCaptureF->at(k), k);
                        inn++;
                    }
                }

                double Hod1CX[inn]; double Hod1CY[inn]; double Hod1CZ[inn];
                double Hod2CX[ext]; double Hod2CY[ext]; double Hod2CZ[ext];

                inn = 0; ext = 0;

                for(int k=0; k<HodCaptureX->size(); k++){
                    if(HodCaptureT->at(k)>0 && HodCaptureVolume->at(k)==2){
                        Hod2CX[ext] = HodCaptureX->at(k);
                        Hod2CY[ext] = HodCaptureY->at(k);
                        Hod2CZ[ext] = HodCaptureZ->at(k);
                        ext++;
                    }
                    if(HodCaptureT->at(k)>0 && HodCaptureVolume->at(k)==1){
                        Hod1CX[inn] = HodCaptureX->at(k);
                        Hod1CY[inn] = HodCaptureY->at(k);
                        Hod1CZ[inn] = HodCaptureZ->at(k);

                        cout << "TEST " << Hod1CX[inn] << " " << Hod1CY[inn] << " " << Hod1CZ[inn] << endl;
                        inn++;
                    }
                }
                if(inn!=0) {
                    g_Hod1captureXYZ[counter] = new TGraph2D("Hod1Capture_"+ID, "Hod1Capture", inn, Hod1CX, Hod1CZ, Hod1CY);
                    flag1[counter] = true;
                }
                else g_Hod1captureXYZ[counter] = new TGraph2D();
                if(ext!=0){
                    g_Hod2captureXYZ[counter] = new TGraph2D("Hod2Capture_"+ID, "Hod2Capture", ext, Hod2CX, Hod2CZ, Hod2CY);
                    flag2[counter] = true;
                } 
                else g_Hod2captureXYZ[counter] = new TGraph2D();

                inn = 0; ext = 0;

            }
            
            



            // PRIMARY PARTICLE POSITION
            h_primaryXY[counter] = new TH2D("h_primaryXY_"+ID,"Primary (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_primaryYZ[counter] = new TH2D("h_primaryYZ_"+ID,"Primary (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_primaryXZ[counter] = new TH2D("h_primaryXZ_"+ID,"Primary (XZ);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);

            // getting the tree entry
            t->GetEntry(entry);

            h_primaryXY[counter]->Fill(MCPosX,MCPosY);
            h_primaryYZ[counter]->Fill(MCPosY,MCPosZ);
            h_primaryXZ[counter]->Fill(MCPosX,MCPosZ);

            g_primaryXYZ[counter] = new TGraph2D("Primary_"+ID, "Primary", 1, &MCPosX, &MCPosZ, &MCPosY);
            double Dx = (MCMomX*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;
            double Dy = (MCMomY*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;
            double Dz = (MCMomZ*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;

            printf(" Primary particle --> X: %f ; Y: %f ; Z: %f E: %f ; M: %f ; Px: %f ; Py: %f ; Pz: %f ; Dx: %f ; Dy: %f ; Dz: %f \n",MCPosX,MCPosY,MCPosZ,MCEnergy,MCMass,MCMomX,MCMomY,MCMomZ,Dx,Dy,Dz);

            counter ++;
        }
    }

    //f_input->Close();
    

    // ------------------- EVENT DISPLAY ------------------- //
    // better select a restricted number of events to visualize

    TFile *f_output = new TFile(outputDir+"/eventDisplay.root","RECREATE");
    f_output->cd();

    TCanvas *canvas[NEvents];
    TCanvas *c2[NEvents];

    for(int i=0;i<NEvents;i++){ // eneter the event range to be displayed

        canvas[i] = new TCanvas(Form("canvas_%lld",IDs[i]), Form("ev_%lld",IDs[i]),800,800);
        canvas[i]->cd();
        canvas[i]->Divide(2,2,0.001,0.001);
        
        // XY PLANE   
        canvas[i]->cd(1);     
            // before capture
        h_ExtScintXY[i]->SetMarkerColor(kRed);
        h_ExtScintXY[i]->SetMarkerSize(1.1);
        h_ExtScintXY[i]->SetMarkerStyle(20);
        h_ExtScintXY[i]->Draw("");
        

        h_InnScintXY[i]->SetMarkerColor(kBlue);
        h_InnScintXY[i]->SetMarkerSize(1.1);
        h_InnScintXY[i]->SetMarkerStyle(20);
        h_InnScintXY[i]->Draw("same");

        h_HeXY[i]->SetMarkerColor(kGreen);
        h_HeXY[i]->SetMarkerSize(1.1);
        h_HeXY[i]->SetMarkerStyle(20);
        h_HeXY[i]->Draw("same");

        h_primaryXY[i]->SetMarkerColor(kViolet);
        h_primaryXY[i]->SetMarkerSize(2.);
        h_primaryXY[i]->SetMarkerStyle(29);
        h_primaryXY[i]->Draw("same");
        
        // helium capture
        if(isBar){
            h_captureXY[i]->SetMarkerColor(kBlack);
            h_captureXY[i]->SetMarkerSize(2.);
            h_captureXY[i]->SetMarkerStyle(29);
            h_captureXY[i]->Draw("same");

            // hodoscopes capture
            if(h_Hod1captureXY[i]->GetEntries()>0){
                h_Hod1captureXY[i]->SetMarkerColor(kOrange+7);
                h_Hod1captureXY[i]->SetMarkerSize(2.);
                h_Hod1captureXY[i]->SetMarkerStyle(29);
                h_Hod1captureXY[i]->Draw("same");
            }
            if(h_Hod2captureXY[i]->GetEntries()>0){
                h_Hod2captureXY[i]->SetMarkerColor(kCyan+1);
                h_Hod2captureXY[i]->SetMarkerSize(2.);
                h_Hod2captureXY[i]->SetMarkerStyle(29);
                h_Hod2captureXY[i]->Draw("same");
            }

                //after capture
            h_ExtScintXY_after[i]->SetMarkerColor(kRed+2);
            h_ExtScintXY_after[i]->SetMarkerSize(1.1);
            h_ExtScintXY_after[i]->SetMarkerStyle(20);
            h_ExtScintXY_after[i]->Draw("same");

            h_InnScintXY_after[i]->SetMarkerColor(kBlue+2);
            h_InnScintXY_after[i]->SetMarkerSize(1.1);
            h_InnScintXY_after[i]->SetMarkerStyle(20);
            h_InnScintXY_after[i]->Draw("same");

            h_HeXY_after[i]->SetMarkerColor(kGreen+2);
            h_HeXY_after[i]->SetMarkerSize(1.1);
            h_HeXY_after[i]->SetMarkerStyle(20);
            h_HeXY_after[i]->Draw("same");
        }
        

        gPad->BuildLegend();



        // YZ PLANE
        canvas[i]->cd(2);
            // before capture
        h_ExtScintYZ[i]->SetMarkerColor(kRed);
        h_ExtScintYZ[i]->SetMarkerSize(1.1);
        h_ExtScintYZ[i]->SetMarkerStyle(20);
        h_ExtScintYZ[i]->Draw("");

        h_InnScintYZ[i]->SetMarkerColor(kBlue);
        h_InnScintYZ[i]->SetMarkerSize(1.1);
        h_InnScintYZ[i]->SetMarkerStyle(20);
        h_InnScintYZ[i]->Draw("same");

        h_HeYZ[i]->SetMarkerColor(kGreen);
        h_HeYZ[i]->SetMarkerSize(1.1);
        h_HeYZ[i]->SetMarkerStyle(20);
        h_HeYZ[i]->Draw("same");

        h_primaryYZ[i]->SetMarkerColor(kViolet);
        h_primaryYZ[i]->SetMarkerSize(2.);
        h_primaryYZ[i]->SetMarkerStyle(29);
        h_primaryYZ[i]->Draw("same");

        if(isBar){
            // helium capture
            h_captureYZ[i]->SetMarkerColor(kBlack);
            h_captureYZ[i]->SetMarkerSize(2.);
            h_captureYZ[i]->SetMarkerStyle(29);
            h_captureYZ[i]->Draw("same");

            // hodoscopes capture
            if(h_Hod1captureYZ[i]->GetEntries()>0){
                h_Hod1captureYZ[i]->SetMarkerColor(kOrange+7);
                h_Hod1captureYZ[i]->SetMarkerSize(2.);
                h_Hod1captureYZ[i]->SetMarkerStyle(29);
                h_Hod1captureYZ[i]->Draw("same");
            }
            if(h_Hod2captureYZ[i]->GetEntries()>0){
                h_Hod2captureYZ[i]->SetMarkerColor(kCyan+1);
                h_Hod2captureYZ[i]->SetMarkerSize(2.);
                h_Hod2captureYZ[i]->SetMarkerStyle(29);
                h_Hod2captureYZ[i]->Draw("same");
            }

                // after capture
            h_ExtScintYZ_after[i]->SetMarkerColor(kRed+2);
            h_ExtScintYZ_after[i]->SetMarkerSize(1.1);
            h_ExtScintYZ_after[i]->SetMarkerStyle(20);
            h_ExtScintYZ_after[i]->Draw("same");

            h_InnScintYZ_after[i]->SetMarkerColor(kBlue+2);
            h_InnScintYZ_after[i]->SetMarkerSize(1.1);
            h_InnScintYZ_after[i]->SetMarkerStyle(20);
            h_InnScintYZ_after[i]->Draw("same");

            h_HeYZ_after[i]->SetMarkerColor(kGreen+2);
            h_HeYZ_after[i]->SetMarkerSize(1.1);
            h_HeYZ_after[i]->SetMarkerStyle(20);
            h_HeYZ_after[i]->Draw("same");

        }


        // XZ PLANE
        canvas[i]->cd(3);
            // before capture
        h_ExtScintXZ[i]->SetMarkerColor(kRed);
        h_ExtScintXZ[i]->SetMarkerSize(1.1);
        h_ExtScintXZ[i]->SetMarkerStyle(20);
        h_ExtScintXZ[i]->Draw("");

        h_InnScintXZ[i]->SetMarkerColor(kBlue);
        h_InnScintXZ[i]->SetMarkerSize(1.1);
        h_InnScintXZ[i]->SetMarkerStyle(20);
        h_InnScintXZ[i]->Draw("same");

        h_HeXZ[i]->SetMarkerColor(kGreen);
        h_HeXZ[i]->SetMarkerSize(1.1);
        h_HeXZ[i]->SetMarkerStyle(20);
        h_HeXZ[i]->Draw("same");

        h_primaryXZ[i]->SetMarkerColor(kViolet);
        h_primaryXZ[i]->SetMarkerSize(2.);
        h_primaryXZ[i]->SetMarkerStyle(29);
        h_primaryXZ[i]->Draw("same");

        if(isBar){
            // helium capture
            h_captureXZ[i]->SetMarkerColor(kBlack);
            h_captureXZ[i]->SetMarkerSize(2.);
            h_captureXZ[i]->SetMarkerStyle(29);
            h_captureXZ[i]->Draw("same");
            // hodoscopes capture
            if(h_Hod1captureXZ[i]->GetEntries()>0){
                h_Hod1captureXZ[i]->SetMarkerColor(kOrange+7);
                h_Hod1captureXZ[i]->SetMarkerSize(2.);
                h_Hod1captureXZ[i]->SetMarkerStyle(29);
                h_Hod1captureXZ[i]->Draw("same");
            }
            if(h_Hod2captureXZ[i]->GetEntries()>0){
                h_Hod2captureXZ[i]->SetMarkerColor(kCyan+1);
                h_Hod2captureXZ[i]->SetMarkerSize(2.);
                h_Hod2captureXZ[i]->SetMarkerStyle(29);
                h_Hod2captureXZ[i]->Draw("same");
            }

                // after capture
            h_ExtScintXZ_after[i]->SetMarkerColor(kRed+2);
            h_ExtScintXZ_after[i]->SetMarkerSize(1.1);
            h_ExtScintXZ_after[i]->SetMarkerStyle(20);
            h_ExtScintXZ_after[i]->Draw("same");

            h_InnScintXZ_after[i]->SetMarkerColor(kBlue+2);
            h_InnScintXZ_after[i]->SetMarkerSize(1.1);
            h_InnScintXZ_after[i]->SetMarkerStyle(20);
            h_InnScintXZ_after[i]->Draw("same");

            h_HeXZ_after[i]->SetMarkerColor(kGreen+2);
            h_HeXZ_after[i]->SetMarkerSize(1.1);
            h_HeXZ_after[i]->SetMarkerStyle(20);
            h_HeXZ_after[i]->Draw("same");
        }


        // ET PLANE
        canvas[i]->cd(4);
        
        h_ExtScintET[i]->SetMarkerColor(kRed);
        h_ExtScintET[i]->SetMarkerSize(1.1);
        h_ExtScintET[i]->SetMarkerStyle(20);
        h_ExtScintET[i]->Draw("");

        h_InnScintET[i]->SetMarkerColor(kBlue);
        h_InnScintET[i]->SetMarkerSize(1.1);
        h_InnScintET[i]->SetMarkerStyle(20);
        h_InnScintET[i]->Draw("same");

        h_HeET[i]->SetMarkerColor(kGreen);
        h_HeET[i]->SetMarkerSize(1.1);
        h_HeET[i]->SetMarkerStyle(20);
        h_HeET[i]->Draw("same");


        if(isBar){
            // helium capture
            h_captureET[i]->SetMarkerColor(kBlack);
            h_captureET[i]->SetMarkerSize(2.);
            h_captureET[i]->SetMarkerStyle(29);
            h_captureET[i]->Draw("same");

            // hodoscopes capture
            if(h_Hod1captureET[i]->GetEntries()>0){
                h_Hod1captureET[i]->SetMarkerColor(kOrange+7);
                h_Hod1captureET[i]->SetMarkerSize(2.);
                h_Hod1captureET[i]->SetMarkerStyle(29);
                h_Hod1captureET[i]->Draw("same");
            }
            if(h_Hod2captureET[i]->GetEntries()>0){
                h_Hod2captureET[i]->SetMarkerColor(kCyan+1);
                h_Hod2captureET[i]->SetMarkerSize(2.);
                h_Hod2captureET[i]->SetMarkerStyle(29);
                h_Hod2captureET[i]->Draw("same");
            }

                // after capture
            h_ExtScintET_after[i]->SetMarkerColor(kRed+2);
            h_ExtScintET_after[i]->SetMarkerSize(1.1);
            h_ExtScintET_after[i]->SetMarkerStyle(20);
            h_ExtScintET_after[i]->Draw("same");

            h_InnScintET_after[i]->SetMarkerColor(kBlue+2);
            h_InnScintET_after[i]->SetMarkerSize(1.1);
            h_InnScintET_after[i]->SetMarkerStyle(20);
            h_InnScintET_after[i]->Draw("same");

            h_HeET_after[i]->SetMarkerColor(kGreen+2);
            h_HeET_after[i]->SetMarkerSize(1.1);
            h_HeET_after[i]->SetMarkerStyle(20);
            h_HeET_after[i]->Draw("same");
        }
        

            //saving
        //canvasET[i].SaveAs(Form("images/ETplane%i.png",i));

        canvas[i]->Write();

        c2[i] = new TCanvas(Form("c2_%lld",IDs[i]), Form("c2_%lld",IDs[i]),800,800);
        c2[i]->cd();

        g_ExtScintXYZ[i]->SetMarkerColor(kRed);
        g_ExtScintXYZ[i]->SetMarkerSize(1.4);
        g_ExtScintXYZ[i]->SetMarkerStyle(20);
        g_ExtScintXYZ[i]->Draw("P");

        g_InnScintXYZ[i]->SetMarkerColor(kBlue);
        g_InnScintXYZ[i]->SetMarkerSize(1.4);
        g_InnScintXYZ[i]->SetMarkerStyle(20);
        g_InnScintXYZ[i]->Draw("P, same");

        g_HeXYZ[i]->SetMarkerColor(kGreen);
        g_HeXYZ[i]->SetMarkerSize(1.1);
        g_HeXYZ[i]->SetMarkerStyle(20);
        g_HeXYZ[i]->Draw("P, same");

        g_primaryXYZ[i]->SetMarkerColor(kViolet);
        g_primaryXYZ[i]->SetMarkerSize(2.);
        g_primaryXYZ[i]->SetMarkerStyle(29);
        g_primaryXYZ[i]->Draw("P, same");

        if(isBar){
            g_ExtScintXYZ_after[i]->SetMarkerColor(kRed+2);
            g_ExtScintXYZ_after[i]->SetMarkerSize(1.4);
            g_ExtScintXYZ_after[i]->SetMarkerStyle(20);
            g_ExtScintXYZ_after[i]->Draw("P, same");

            g_InnScintXYZ_after[i]->SetMarkerColor(kBlue+2);
            g_InnScintXYZ_after[i]->SetMarkerSize(1.4);
            g_InnScintXYZ_after[i]->SetMarkerStyle(20);
            g_InnScintXYZ_after[i]->Draw("P, same");

            g_HeXYZ_after[i]->SetMarkerColor(kGreen+2);
            g_HeXYZ_after[i]->SetMarkerSize(1.1);
            g_HeXYZ_after[i]->SetMarkerStyle(20);
            g_HeXYZ_after[i]->Draw("P, same");

            g_captureXYZ[i]->SetMarkerColor(kBlack);
            g_captureXYZ[i]->SetMarkerSize(2.);
            g_captureXYZ[i]->SetMarkerStyle(29);
            g_captureXYZ[i]->Draw("P, same");

            if(flag1[i]){
                g_Hod1captureXYZ[i]->SetMarkerColor(kOrange+7);
                g_Hod1captureXYZ[i]->SetMarkerSize(2.);
                g_Hod1captureXYZ[i]->SetMarkerStyle(29);
                g_Hod1captureXYZ[i]->Draw("P, same");
            }
            if(flag2[i]){
                g_Hod2captureXYZ[i]->SetMarkerColor(kCyan+1);
                g_Hod2captureXYZ[i]->SetMarkerSize(2.);
                g_Hod2captureXYZ[i]->SetMarkerStyle(29);
                g_Hod2captureXYZ[i]->Draw("P, same");
            }
        }

        c2[i]->BuildLegend();
        c2[i]->Write();
    }
    f_output->Close();
    return;
}