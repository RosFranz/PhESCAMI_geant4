#include <vector>
#include "globals.hh"

class MyDet
{
public:

  static MyDet* gethead()
  {
    if (!_head) _head = new MyDet();
    return _head;
  };

  ~MyDet(){_head=0;};
  MyDet(){SetDefault();};

  ///constructor
  G4int MCEvnum;
  G4int TRGEvnum;
  G4double GenSurface;

  //CAPTURE IN HELIUM
  std::vector<G4double> MyCaptureT;
  std::vector<G4double> MyCaptureX;
  std::vector<G4double> MyCaptureY;
  std::vector<G4double> MyCaptureZ;
  std::vector<G4double> MyCaptureE;
  std::vector<G4double> MyCaptureF;
  std::vector<G4double> MyCaptureM;
  std::vector<G4int> MyCaptureCopyNo;

    //daughters
  std::vector<G4double> MyCaptureDaughterM;
  std::vector<G4double> MyCaptureDaughterX;
  std::vector<G4double> MyCaptureDaughterY;
  std::vector<G4double> MyCaptureDaughterZ;
  std::vector<G4long> MyCaptureDaughterPDGcode;

  //CAPTURE IN HODOSCOPES
  std::vector<G4double> HodCaptureT;
  std::vector<G4double> HodCaptureX;
  std::vector<G4double> HodCaptureY;
  std::vector<G4double> HodCaptureZ;
  std::vector<G4double> HodCaptureE;
  std::vector<G4double> HodCaptureF;
  std::vector<G4double> HodCaptureM;
  std::vector<G4int> HodCaptureVolume;

    //daughters
  std::vector<G4double> HodCaptureDaughterM;
  std::vector<G4double> HodCaptureDaughterX;
  std::vector<G4double> HodCaptureDaughterY;
  std::vector<G4double> HodCaptureDaughterZ;
  std::vector<G4long> HodCaptureDaughterPDGcode;

  G4double MyDelay;
  G4double MyMass;
  
  G4double MyVesselThick;
  G4double MyVesselRadius;
  G4double MyVesselLength;
  G4double MyVesselTubsLength;
  G4double MyVesselPlugLength;
  G4double MyVesselPlugRadius;
  
  G4double MyHeDensity;
  G4double MyHeTres;
  G4double MyHeDist;

  G4double MySCTres;
  G4double MySCDist;
  G4double MyDistance;
  G4double MySthick;
  G4double MyRS2;
  G4int MyNBars2;
  G4double MyRS3;
  G4int MyNBars3;

  G4int Visible, OnlyHod;
  G4int antiparticle;


private:

  static  MyDet* _head;

  void SetDefault()
  {
    MCEvnum=0;
    TRGEvnum=0;
    GenSurface = 0;
    //Helium capture
    MyCaptureT.resize(1,0);
    MyCaptureX.resize(1,0);
    MyCaptureY.resize(1,0);
    MyCaptureZ.resize(1,0);
    MyCaptureE.resize(1,0);
    MyCaptureF.resize(1,0);
    MyCaptureM.resize(1,0);
    MyCaptureCopyNo.resize(1,0);
    //daughters
    MyCaptureDaughterM.resize(1,0);
    MyCaptureDaughterX.resize(1,0);
    MyCaptureDaughterY.resize(1,0);
    MyCaptureDaughterZ.resize(1,0);
    MyCaptureDaughterPDGcode.resize(1,0);

    //Hodoscope capture
    HodCaptureT.resize(1,0);
    HodCaptureX.resize(1,0);
    HodCaptureY.resize(1,0);
    HodCaptureZ.resize(1,0);
    HodCaptureE.resize(1,0);
    HodCaptureF.resize(1,0);
    HodCaptureM.resize(1,0);
    HodCaptureVolume.resize(1,0);
    //daughters
    HodCaptureDaughterM.resize(1,0);
    HodCaptureDaughterX.resize(1,0);
    HodCaptureDaughterY.resize(1,0);
    HodCaptureDaughterZ.resize(1,0);
    HodCaptureDaughterPDGcode.resize(1,0);
    

    MyDelay = 900000.; //900 us fixed delay for Z=-1 massive particle capture
    MyMass = 900.; //900MeV mass of Z=-1 massive particle capture
    

    //vessel (He) dimensions
    MyVesselThick = 11.; //  vessel thickness
    MyVesselRadius = 216.;//432./2.; // 216mm vessel radius
    MyVesselLength = 660.; // 660mm vessel length
    MyVesselTubsLength = (MyVesselLength - (MyVesselRadius*2.))/2.; // 148mm vessel tubs half length
    MyVesselPlugLength = 31.; // 62/2 mm vessel plug half legth
    MyVesselPlugRadius = 31.; // 31 mm vessel plug radius


    //He resolutions and specifications
    MyHeDensity = (0.114*40./68.)*(310./400.); // g/cm3 He @ pressure 310 bar
    MyHeTres = 0.2; // 200 ps He hit resolution                                                         
    MyHeDist = 20.; // (1= 1mm) 2cm He hit resolution 

    //Scintillators resolutions
    MySCDist = 5.; // 0.5cm Plastic hit resolution  
    MySCTres = 0.1; // 100 ps plastic hit resolution
    //dimensions
    MySthick = 2.; // 4 mm Scint half thickness
    MyDistance = 200.; // 200 mm distance between the scintillator layers
    //outer layer
    MyRS3 = 1510.; // mm External ScG4int Half Lenght  
    MyNBars3 = 64; // 64x inner bars
    //inner layer
    MyRS2 = MyRS3 - MyDistance; // mm internal ScG4int Half Lenght
    MyNBars2 = 64; // 64x inner bars
    
    Visible = 0; // open one side if = 1
    OnlyHod = 0; // flag to display only one side of hodoscopes
    antiparticle = 1; // flag to choose wich ntuple produce in run action (1 for antiparticles)
    return;
  };

};

