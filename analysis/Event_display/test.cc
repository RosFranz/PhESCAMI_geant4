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

void test(TString inputFile, TString outputDir)
{

    // ---- 1st parameter is the INPUT FILE ---- 2nd parameter the OUTPUT DIR ----- //

    //ANTI-PROTON   
    //ANTI-Deuteron Newdbar.root (capture as vectors)

    TFile *f_input = new TFile(inputFile,"OPEN");

    gStyle->SetOptStat(0);

    if (!f_input){
        printf("File not found");
        printf("Suggested file: \n");
        printf("\t Anti-proton: \n");
        printf("\t Anti-deuteron: Newdbar.root/\n");
        return;
    } 


    // ----------------- GETTING THE TREE ----------------- //
    TTree *t = (TTree*)f_input->Get("Hits");
    const Long64_t entries = t->GetEntries();
    cout << "Number of entries: " << entries << endl;

    // ---------- CONSTANT USED FOR THE EVENT DISPLAY ---------- //
    const int NEvents = 1;//5;
    Long64_t IDs[NEvents] = {0};

                            //3269, 5454, 22625, 46942,
    /* 
                            74753, 78994, 96772, 100166, 101691,
                            109183, 114848, 150089, 153127, 155020, 166475, 177413, 197132, 209315, 212156,
                            217625, 224344, 230889, 252254, 288616, 291680}; //array of event IDs
    */
    
    TTreeReader reader(t);

    // MC TRUTH

    double MCEnergy, MCMomentum;
    double MCMomX, MCMomY, MCMomZ;
    double MCPosX, MCPosY, MCPosZ;
    double MCMass, MCCharge;

        // setting the branches
    t->SetBranchAddress("MCEnergy",&MCEnergy);
    t->SetBranchAddress("MCMomentum",&MCMomentum);
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


    // PRIMARY PARTICLE POSITION
        // array of TH2
    TH2D* h_primaryXY[NEvents];
    TH2D* h_primaryYZ[NEvents];
    TH2D* h_primaryXZ[NEvents];

    TGraph2D* g_primaryXYZ[NEvents];


    // POSITION, ENERGY AND TIME OF THE CAPTURE 
        //IN HELIUM
    TTreeReaderArray<double> CaptureX(reader, "CaptureX");
    TTreeReaderArray<double> CaptureY(reader, "CaptureY");
    TTreeReaderArray<double> CaptureZ(reader, "CaptureZ");
    TTreeReaderArray<double> CaptureE(reader, "CaptureE");
    TTreeReaderArray<double> CaptureT(reader, "CaptureT");
    TTreeReaderArray<double> CaptureM(reader, "CaptureM");
    TTreeReaderArray<double> CaptureF(reader, "CaptureF");

    TH2D* h_captureXY[NEvents];
    TH2D* h_captureYZ[NEvents];
    TH2D* h_captureXZ[NEvents];
    TH2D* h_captureET[NEvents];

    TGraph2D* g_captureXYZ[NEvents];

        //IN THE HODOSCOPES
    TTreeReaderArray<double> HodCaptureX(reader, "HodCaptureX");
    TTreeReaderArray<double> HodCaptureY(reader, "HodCaptureY");
    TTreeReaderArray<double> HodCaptureZ(reader, "HodCaptureZ");
    TTreeReaderArray<double> HodCaptureE(reader, "HodCaptureE");
    TTreeReaderArray<double> HodCaptureT(reader, "HodCaptureT");
    TTreeReaderArray<double> HodCaptureM(reader, "HodCaptureM");
    TTreeReaderArray<double> HodCaptureF(reader, "HodCaptureF");
    TTreeReaderArray<int> HodCaptureVolume(reader, "HodCaptureVolume");

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
    TTreeReaderArray<double> Hod2E(reader, "Hod2E");
    TTreeReaderArray<double> Hod2T(reader, "Hod2T");
    TTreeReaderArray<double> Hod2X(reader, "Hod2X");
    TTreeReaderArray<double> Hod2Y(reader, "Hod2Y");
    TTreeReaderArray<double> Hod2Z(reader, "Hod2Z");

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
    TTreeReaderArray<double> Hod3E(reader, "Hod3E");
    TTreeReaderArray<double> Hod3T(reader, "Hod3T");
    TTreeReaderArray<double> Hod3X(reader, "Hod3X");
    TTreeReaderArray<double> Hod3Y(reader, "Hod3Y");
    TTreeReaderArray<double> Hod3Z(reader, "Hod3Z");

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
    TTreeReaderArray<double> HGasE(reader,"HGasE");
    TTreeReaderArray<double> HGasX(reader,"HGasX");
    TTreeReaderArray<double> HGasY(reader,"HGasY");
    TTreeReaderArray<double> HGasZ(reader,"HGasZ");
    TTreeReaderArray<double> HGasT(reader,"HGasT");
    TTreeReaderArray<int> HGasCopyNo(reader,"HGasCopyNo");

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

    int j = 0;
    int before = 0; int after = 0;
    int ext = 0; int inn = 0;
    bool flag1[NEvents]; bool flag2[NEvents];
    for(int i = 0; i < NEvents; i++){
        flag1[i] = false;
        flag2[i] = false;
    }

    Long64_t entry = 0;
    int counter = 0;
    
   
    while(reader.Next()){

        std::stringstream strID;
        TString ID;

        if(entry == IDs[counter]){

            printf("EVENT: %lli \n",entry);

                //setting the histo names
            strID << entry;
            ID = strID.str();

            // EXTERNAL SCINTILLATOR
                // before the capture
            h_ExtScintXY[counter] = new TH2D("h_ExtScintXY_"+ID,"Ext. scint.(XY) ID ="+ID+";X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintYZ[counter] = new TH2D("h_ExtScintYZ_"+ID,"Ext. scint.(YZ) ID ="+ID+";Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintXZ[counter] = new TH2D("h_ExtScintXZ_"+ID,"Ext. scint.(XZ) ID ="+ID+";X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintET[counter] = new TH2D("h_ExtScintET_"+ID,"Ext. scint.(ET) ID ="+ID+";E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // after the capture
            h_ExtScintXY_after[counter] = new TH2D("h_ExtScintXY_after_"+ID,"Ext. scint.(XY) (T>100 #mu s);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintYZ_after[counter] = new TH2D("h_ExtScintYZ_after_"+ID,"Ext. scint.(YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintXZ_after[counter] = new TH2D("h_ExtScintXZ_after_"+ID,"Ext. scint.(XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_ExtScintET_after[counter] = new TH2D("h_ExtScintET_after_"+ID,"Ext. scint.(ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                // FILLING
            for(auto i:Hod3X) {
                if(Hod3T.At(j)/pow(10.,3) < 100.){
                    h_ExtScintXY[counter]->Fill(i,Hod3Y.At(j));
                    h_ExtScintYZ[counter]->Fill(Hod3Y.At(j),Hod3Z.At(j));
                    h_ExtScintXZ[counter]->Fill(i,Hod3Z.At(j));
                    h_ExtScintET[counter]->Fill(Hod3E.At(j),Hod3T.At(j)/pow(10.,3));
                    before++;
                }
                else{
                    h_ExtScintXY_after[counter]->Fill(i,Hod3Y.At(j));
                    h_ExtScintYZ_after[counter]->Fill(Hod3Y.At(j),Hod3Z.At(j));
                    h_ExtScintXZ_after[counter]->Fill(i,Hod3Z.At(j));
                    h_ExtScintET_after[counter]->Fill(Hod3E.At(j),Hod3T.At(j)/pow(10.,3));
                    after++;
                }
                j++;
            }
            double H3X[before+2]; H3X[0]=-maxDim; H3X[before+1]=maxDim;
            double H3Y[before+2]; H3Y[0]=-maxDim; H3Y[before+1]=maxDim;
            double H3Z[before+2]; H3Z[0]=-maxDim; H3Z[before+1]=maxDim;

            double H3X_after[after+2]; H3X_after[0]=-maxDim; H3X_after[after+1]=maxDim;
            double H3Y_after[after+2]; H3Y_after[0]=-maxDim; H3Y_after[after+1]=maxDim;
            double H3Z_after[after+2]; H3Z_after[0]=-maxDim; H3Z_after[after+1]=maxDim;

            j = 0; before = 0; after = 0;

            for(auto i:Hod3X){
                if(Hod3T.At(j)/pow(10.,3) < 100.){
                    H3X[before+1]=i;
                    H3Y[before+1]=Hod3Y.At(j);
                    H3Z[before+1]=Hod3Z.At(j);
                    before++;
                }
                else{
                    H3X_after[after+1]=i;
                    H3Y_after[after+1]=Hod3Y.At(j);
                    H3Z_after[after+1]=Hod3Z.At(j);
                    after++;   
                }
                j++;
            }
            
            g_ExtScintXYZ[counter] = new TGraph2D("Ext_"+ID, "Ext_"+ID, before+2, H3X, H3Y, H3Z);
            g_ExtScintXYZ[counter]->GetXaxis()->SetTitle("X [mm]"); g_ExtScintXYZ[counter]->GetYaxis()->SetTitle("Y [mm]"); g_ExtScintXYZ[counter]->GetZaxis()->SetTitle("Z [mm]");
            g_ExtScintXYZ_after[counter] = new TGraph2D("Ext_after_"+ID, "Ext_after", after+2, H3X_after, H3Y_after, H3Z_after);

            j = 0; before = 0; after = 0;

            // INNER SCINTILLATOR
                // before the capture
            h_InnScintXY[counter] = new TH2D("h_InnScintXY_"+ID,"Inn. scint(XY);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintYZ[counter] = new TH2D("h_InnScintYZ_"+ID,"Inn. scint(YZ);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintXZ[counter] = new TH2D("h_InnScintXZ_"+ID,"Inn. scint(XZ);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintET[counter] = new TH2D("h_InnScintET_"+ID,"Inn. scint(ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // after the capture
            h_InnScintXY_after[counter] = new TH2D("h_InnScintXY_after_"+ID,"Inn. scint(XY) (T>100 #mu s);X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintYZ_after[counter] = new TH2D("h_InnScintYZ_after_"+ID,"Inn. scint(YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintXZ_after[counter] = new TH2D("h_InnScintXZ_after_"+ID,"Inn. scint(XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_InnScintET_after[counter] = new TH2D("h_InnScintET_after_"+ID,"Inn. scint(ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                // FILLING
            for(auto i:Hod2X){
                if(Hod2T.At(j)/pow(10.,3) < 100.){
                    h_InnScintXY[counter]->Fill(i,Hod2Y.At(j));
                    h_InnScintYZ[counter]->Fill(Hod2Y.At(j),Hod2Z.At(j));
                    h_InnScintXZ[counter]->Fill(i,Hod2Z.At(j));
                    h_InnScintET[counter]->Fill(Hod2E.At(j),Hod2T.At(j)/pow(10.,3));
                    before++;
                }
                else{
                    h_InnScintXY_after[counter]->Fill(i,Hod2Y.At(j));
                    h_InnScintYZ_after[counter]->Fill(Hod2Y.At(j),Hod2Z.At(j));
                    h_InnScintXZ_after[counter]->Fill(i,Hod2Z.At(j));
                    h_InnScintET_after[counter]->Fill(Hod2E.At(j),Hod2T.At(j)/pow(10.,3));
                    after++;
                }
                j++;
            } 

            double H2X[before];// H2X[0]=-maxDim; H2X[before+1]=maxDim;
            double H2Y[before];// H2Y[0]=-maxDim; H2Y[before+1]=maxDim;
            double H2Z[before];// H2Z[0]=-maxDim; H2Z[before+1]=maxDim;

            double H2X_after[after];// H2X_after[0]=-maxDim; H2X_after[after+1]=maxDim;
            double H2Y_after[after];// H2Y_after[0]=-maxDim; H2Y_after[after+1]=maxDim;
            double H2Z_after[after];// H2Z_after[0]=-maxDim; H2Z_after[after+1]=maxDim;

            j = 0; before = 0; after = 0;

            for(auto i:Hod2X){
                if(Hod2T.At(j)/pow(10.,3) < 100.){
                    H2X[before]=i;
                    H2Y[before]=Hod2Y.At(j);
                    H2Z[before]=Hod2Z.At(j);
                    before++;
                }
                else{
                    H2X_after[after]=i;
                    H2Y_after[after]=Hod2Y.At(j);
                    H2Z_after[after]=Hod2Z.At(j);
                    after++;
                }
                j++;
            }

            g_InnScintXYZ[counter] = new TGraph2D("Inn_"+ID, "Inn", before, H2X, H2Y, H2Z);
            g_InnScintXYZ_after[counter] = new TGraph2D("Inn_after_"+ID, "Inn_after", after, H2X_after, H2Y_after, H2Z_after);

            j = 0; before = 0; after = 0;

            // HELIUM SCINTILLATOR
                // before the capture
            h_HeXY[counter] = new TH2D("h_HeXY_"+ID,"He (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeYZ[counter] = new TH2D("h_HeYZ_"+ID,"He (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeXZ[counter] = new TH2D("h_HeXZ_"+ID,"He (XZ);Z [mm];X [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeET[counter] = new TH2D("h_HeET_"+ID,"He (ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);
                // He after the capture
            h_HeXY_after[counter] = new TH2D("h_HeXY_after_"+ID,"He (XY) (T>100 #mu s);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeYZ_after[counter] = new TH2D("h_HeYZ_after_"+ID,"He (YZ) (T>100 #mu s);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeXZ_after[counter] = new TH2D("h_HeXZ_after_"+ID,"He (XZ) (T>100 #mu s);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_HeET_after[counter] = new TH2D("h_HeET_after_"+ID,"He (ET) (T>100 #mu s);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

                // FILLING
            for(auto i:HGasX){
                if(HGasT.At(j)/pow(10.,3) < 100){
                    h_HeXY[counter]->Fill(i,HGasY.At(j));
                    h_HeYZ[counter]->Fill(HGasY.At(j),HGasZ.At(j));
                    h_HeXZ[counter]->Fill(i,HGasZ.At(j));
                    h_HeET[counter]->Fill(HGasE.At(j),HGasT.At(j)/pow(10.,3));
                    before++;
                }
                else {
                    h_HeXY_after[counter]->Fill(i,HGasY.At(j));
                    h_HeYZ_after[counter]->Fill(HGasY.At(j),HGasZ.At(j));
                    h_HeXZ_after[counter]->Fill(i,HGasZ.At(j));
                    h_HeET_after[counter]->Fill(HGasE.At(j),HGasT.At(j)/pow(10.,3));
                    after++;
                }
                j++;
            }

            double GasX[before];// GasX[0]=-maxDim; GasX[before+1]=maxDim;
            double GasY[before];// GasY[0]=-maxDim; GasY[before+1]=maxDim;
            double GasZ[before];// GasZ[0]=-maxDim; GasZ[before+1]=maxDim;

            double GasX_after[after];// GasX_after[0]=-maxDim; GasX_after[after+1]=maxDim;
            double GasY_after[after];// GasY_after[0]=-maxDim; GasY_after[after+1]=maxDim;
            double GasZ_after[after];// GasZ_after[0]=-maxDim; GasZ_after[after+1]=maxDim;

            j = 0; before = 0; after = 0;

            for(auto i:HGasX){
                if(HGasT.At(j)/pow(10.,3) < 100){
                    GasX[before]=i;
                    GasY[before]=HGasY.At(j);
                    GasZ[before]=HGasZ.At(j);
                    before++;
                }
                else {
                    GasX_after[after]=i;
                    GasY_after[after]=HGasY.At(j);
                    GasZ_after[after]=HGasZ.At(j);
                    after++;
                }
                j++;
            }

            g_HeXYZ[counter] = new TGraph2D("He_"+ID, "He", before, GasX, GasY, GasZ);
            g_HeXYZ_after[counter] = new TGraph2D("He_after_"+ID, "He_after", after, GasX_after, GasY_after, GasZ_after);

            j = 0; before = 0; after = 0;


            // HELIUM CAPTURE POSITION, ENERGY AND TIME 
            h_captureXY[counter] = new TH2D("h_captureXY_"+ID,"Capture (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_captureYZ[counter] = new TH2D("h_captureYZ_"+ID,"Capture (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_captureXZ[counter] = new TH2D("h_captureXZ_"+ID,"Capture (XZ);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_captureET[counter] = new TH2D("h_captureET_"+ID,"Capture (ET);E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

            double HeCX[CaptureX.GetSize()];// HeCX[0]=-maxDim; HeCX[CaptureX.GetSize()+1]=maxDim;
            double HeCY[CaptureX.GetSize()];// HeCY[0]=-maxDim; HeCY[CaptureX.GetSize()+1]=maxDim;
            double HeCZ[CaptureX.GetSize()];// HeCZ[0]=-maxDim; HeCZ[CaptureX.GetSize()+1]=maxDim;

                // FILLING
            for(auto i:CaptureX){
                h_captureXY[counter]->Fill(i,CaptureY.At(j));
                h_captureYZ[counter]->Fill(CaptureY.At(j),CaptureZ.At(j));
                h_captureXZ[counter]->Fill(i,CaptureZ.At(j));
                h_captureET[counter]->Fill(CaptureE.At(j),CaptureT.At(j)/pow(10.,3));
                HeCX[j]=i;
                HeCY[j]=CaptureY.At(j);
                HeCZ[j]=CaptureZ.At(j);
                printf(" Helium capture --> X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",CaptureX.At(j),CaptureY.At(j),CaptureZ.At(j),CaptureE.At(j),CaptureT.At(j),CaptureM.At(j),CaptureF.At(j));
                if(j==1){
                    printf("!!!! --- SECOND CAPTURE IDENTIFIED --- !!!!\n");
                    printf("X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",CaptureX.At(j),CaptureY.At(j),CaptureZ.At(j),CaptureE.At(j),CaptureT.At(j),CaptureM.At(j),CaptureF.At(j));
                    //printf("Primary particle charge: %f ; mass: %f\n", MCCharge, MCMass);
                }
                j++;
            }

            g_captureXYZ[counter] = new TGraph2D("HeCapture_"+ID, "HeCapture", CaptureX.GetSize(), HeCX, HeCY, HeCZ);

            j=0;

            // HODOSCOPE CAPTURES
            h_Hod1captureXY[counter] = new TH2D("h_InnScintCaptureXY_"+ID,"Inn. scint. capture ;X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod1captureYZ[counter] = new TH2D("h_InnScintCaptureYZ_"+ID,"Inn. scint. capture;Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod1captureXZ[counter] = new TH2D("h_InnScintCaptureXZ_"+ID,"Inn. scint. capture;X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod1captureET[counter] = new TH2D("h_InnScintCaptureET_"+ID,"Inn. scint. capture;E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

            h_Hod2captureXY[counter] = new TH2D("h_ExtScintCaptureXY_"+ID,"Ext. scint. capture ;X [mm];Y [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod2captureYZ[counter] = new TH2D("h_ExtScintCaptureYZ_"+ID,"Ext. scint. capture ;Y [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod2captureXZ[counter] = new TH2D("h_ExtScintCaptureXZ_"+ID,"Ext. scint. capture ;X [mm];Z [mm]",nBinsSCDist,-maxDim,maxDim, nBinsSCDist,-maxDim,maxDim);
            h_Hod2captureET[counter] = new TH2D("h_ExtScintCaptureET_"+ID,"Ext. scint. capture ;E [MeV];T [#mu s]",300,0.,30., 10000,0., maxTime/2.);

            for(auto i:HodCaptureX){
                if(HodCaptureT.At(j)>0 && HodCaptureVolume.At(j)==2){
                    h_Hod2captureXY[counter]->Fill(i,HodCaptureY.At(j));
                    h_Hod2captureYZ[counter]->Fill(HodCaptureY.At(j),HodCaptureZ.At(j));
                    h_Hod2captureXZ[counter]->Fill(i,HodCaptureZ.At(j));
                    h_Hod2captureET[counter]->Fill(HodCaptureE.At(j),HodCaptureT.At(j)/pow(10.,3));
                    printf(" EXTERNAL HOD CAPTURE X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f \n",HodCaptureX.At(j),HodCaptureY.At(j),HodCaptureZ.At(j),HodCaptureE.At(j),HodCaptureT.At(j),HodCaptureM.At(j),HodCaptureF.At(j));
                    ext++;
                }
                if(HodCaptureT.At(j)>0 && HodCaptureVolume.At(j)==1){
                    h_Hod1captureXY[counter]->Fill(i,HodCaptureY.At(j));
                    h_Hod1captureYZ[counter]->Fill(HodCaptureY.At(j),HodCaptureZ.At(j));
                    h_Hod1captureXZ[counter]->Fill(i,HodCaptureZ.At(j));
                    h_Hod1captureET[counter]->Fill(HodCaptureE.At(j),HodCaptureT.At(j)/pow(10.,3));
                    printf(" INTERNAL HOD CAPTURE X: %f ; Y: %f ; Z: %f E: %f ; T: %f ; M: %f ; F: %f ; j: %i\n",HodCaptureX.At(j),HodCaptureY.At(j),HodCaptureZ.At(j),HodCaptureE.At(j),HodCaptureT.At(j),HodCaptureM.At(j),HodCaptureF.At(j), j);
                    inn++;
                }
                j++;
            }

            double Hod1CX[inn]; double Hod1CY[inn]; double Hod1CZ[inn];
            double Hod2CX[ext]; double Hod2CY[ext]; double Hod2CZ[ext];

            j = 0; inn = 0; ext = 0;

            for(auto i:HodCaptureX){
                if(HodCaptureT.At(j)>0 && HodCaptureVolume.At(j)==2){
                    Hod2CX[ext] = i;
                    Hod2CY[ext] = HodCaptureY.At(j);
                    Hod2CZ[ext] = HodCaptureZ.At(j);
                    ext++;
                }
                if(HodCaptureT.At(j)>0 && HodCaptureVolume.At(j)==1){
                    Hod1CX[inn] = i;
                    Hod1CY[inn] = HodCaptureY.At(j);
                    Hod1CZ[inn] = HodCaptureZ.At(j);

                    cout << "TEST " << Hod1CX[inn] << " " << Hod1CY[inn] << " " << Hod1CZ[inn] << endl;
                    inn++;
                }
                j++;
            }

            if(inn!=0) {
                g_Hod1captureXYZ[counter] = new TGraph2D("Hod1Capture_"+ID, "Hod1Capture", inn, Hod1CX, Hod1CY, Hod1CZ);
                flag1[counter] = true;
            }
            else g_Hod1captureXYZ[counter] = new TGraph2D();
            if(ext!=0){
                g_Hod2captureXYZ[counter] = new TGraph2D("Hod2Capture_"+ID, "Hod2Capture", ext, Hod2CX, Hod2CY, Hod2CZ);
                flag2[counter] = true;
            } 
            else g_Hod2captureXYZ[counter] = new TGraph2D();

            j = 0; inn = 0; ext = 0;

            // PRIMARY PARTICLE POSITION
            h_primaryXY[counter] = new TH2D("h_primaryXY_"+ID,"Primary (XY);X [mm];Y [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_primaryYZ[counter] = new TH2D("h_primaryYZ_"+ID,"Primary (YZ);Y [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);
            h_primaryXZ[counter] = new TH2D("h_primaryXZ_"+ID,"Primary (XZ);X [mm];Z [mm]",nBinsHeDist,-maxDim,maxDim, nBinsHeDist,-maxDim,maxDim);

            // getting the tree entry
            t->GetEntry(entry);

            h_primaryXY[counter]->Fill(MCPosX,MCPosY);
            h_primaryYZ[counter]->Fill(MCPosY,MCPosZ);
            h_primaryXZ[counter]->Fill(MCPosX,MCPosZ);

            g_primaryXYZ[counter] = new TGraph2D("Primary_"+ID, "Primary", 1, &MCPosX, &MCPosY, &MCPosZ);
            double Dx = (MCMomX*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;
            double Dy = (MCMomY*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;
            double Dz = (MCMomZ*2.)/(MCEnergy*2.+MCMass)*7.151369*299.792458;

            printf(" Primary particle --> X: %f ; Y: %f ; Z: %f E: %f ; M: %f ; Px: %f ; Py: %f ; Pz: %f ; Dx: %f ; Dy: %f ; Dz: %f \n",MCPosX,MCPosY,MCPosZ,MCEnergy,MCMass,MCMomX,MCMomY,MCMomZ,Dx,Dy,Dz);

            counter ++;
        }

        entry++;
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

            //saving
        //canvasET[i].SaveAs(Form("images/ETplane%i.png",i));

        canvas[i]->Write();

        c2[i] = new TCanvas(Form("c2_%lld",IDs[i]), Form("c2_%lld",IDs[i]),800,800);
        c2[i]->cd();

        g_ExtScintXYZ[i]->SetMarkerColor(kRed);
        g_ExtScintXYZ[i]->SetMarkerSize(1.4);
        g_ExtScintXYZ[i]->SetMarkerStyle(20);
        g_ExtScintXYZ[i]->Draw("P");

        g_ExtScintXYZ_after[i]->SetMarkerColor(kRed+2);
        g_ExtScintXYZ_after[i]->SetMarkerSize(1.4);
        g_ExtScintXYZ_after[i]->SetMarkerStyle(20);
        g_ExtScintXYZ_after[i]->Draw("P, same");

        g_InnScintXYZ[i]->SetMarkerColor(kBlue);
        g_InnScintXYZ[i]->SetMarkerSize(1.4);
        g_InnScintXYZ[i]->SetMarkerStyle(20);
        g_InnScintXYZ[i]->Draw("P, same");

        g_InnScintXYZ_after[i]->SetMarkerColor(kBlue+2);
        g_InnScintXYZ_after[i]->SetMarkerSize(1.4);
        g_InnScintXYZ_after[i]->SetMarkerStyle(20);
        g_InnScintXYZ_after[i]->Draw("P, same");

        g_HeXYZ[i]->SetMarkerColor(kGreen);
        g_HeXYZ[i]->SetMarkerSize(1.1);
        g_HeXYZ[i]->SetMarkerStyle(20);
        g_HeXYZ[i]->Draw("P, same");

        g_HeXYZ_after[i]->SetMarkerColor(kGreen+2);
        g_HeXYZ_after[i]->SetMarkerSize(1.1);
        g_HeXYZ_after[i]->SetMarkerStyle(20);
        g_HeXYZ_after[i]->Draw("P, same");

        g_captureXYZ[i]->SetMarkerColor(kBlack);
        g_captureXYZ[i]->SetMarkerSize(2.);
        g_captureXYZ[i]->SetMarkerStyle(29);
        g_captureXYZ[i]->Draw("P, same");

        g_primaryXYZ[i]->SetMarkerColor(kViolet);
        g_primaryXYZ[i]->SetMarkerSize(2.);
        g_primaryXYZ[i]->SetMarkerStyle(29);
        g_primaryXYZ[i]->Draw("P, same");

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

        c2[i]->BuildLegend();


    }
    f_output->Close();
    return;
}