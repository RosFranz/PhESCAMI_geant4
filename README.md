# PhESCAMI_geant4
PHeSCAMI 2024 Monte Carlo simulation

## Name
PHeSCAMI 2024 Monte Carlo simulation

### Detector description 

Geant4 Monte Carlo simulation to study the sensitivity of a large-acceptance detector based on the detection priciple of the [PHeSCAMI](https://www.tifpa.infn.it/projects/prin2022-phescami/) project.
The overall volume is 3x3x3 m^3, the PHeSCAMI detector is made by two subdetectors: the Time of Flight (TOF) and the Helium Calorimeter (HeCal).

### Detector description (TOF)
The TOF is composed by two segmented concentric cubes. Each cube has 64 slabs of plastic scintillators with a thickness of 4 [mm].
The distance between the two layers is 20 [cm]; the external TOF layer is 3.02 [m] long, while the inner layer is 2.62 [m] long.
The spatial resolution is 5 [mm], while the temporal is 0.2 [ns].
All this information are contained in the sources/src/include/MyDet.hh file. Be aware that the default GEANT 4 unit for length is [mm] and [ns] for time. 
Further info in:
- sources/src/include/ADHDDetectorConstruction.hh
- sources/src/include/ADHDHodoscopeHit.hh
- sources/src/include/ADHDHodoscopeSD.hh 

### Detector description (HeCal)
The HeCal is composed by 75 ellipsoid Helium tanks, containg He at 310 [bar]. All the dimensions, and the Helium density are specified here: sources/src/include/MyDet.hh
Each tank has two alluminum plugs, the spatial resolution is 2 [cm] and 0.2 [ns] for the temporal resolution. 
Further info in:
- sources/src/include/ADHDDetectorConstruction.hh
- sources/src/include/ADHDHeGasHit.hh
- sources/src/include/ADHDHeGasSD.hh 

## Installation
GEANT4 and root are required to run this project. 
The project runs with GEANT4 v11.1.2 and root v6.24. Any other combination has never been tested.
- How to install geant4: https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/index.html. \
    - Video tutorial: https://www.youtube.com/watch?v=Lxb4WZyKeCE&list=PLLybgCU6QCGWgzNYOV0SKen9vqg4KXeVL

- For root see: https://root.cern/install/


## Repository structure
The repo is divided in two folders:
- analysis
- sources

Within sources it recommended to create the build directory and a folder for the simulation output (i.e. ROOT_files). On the path sources/src/ there are all the relevant files for the Monte Carlo simulation. Alle the *.cc files are in sources/src/src, while *.hh files are in sources/src/include. 
In sources/src there are many *.mac (macro) files, these files can be used to control the simulation. For example it is possible to choose the particle type, the energy spectra, beam position, output file name and so on. Most of the macro present by deafult have a very large number of simulated events. 

Some reconstruction and analysis macro are stored in the analysis folder. The most relevant one are inside analysis/macros and relies on root. 

## Usage
### First time 
First usage:
- Clone the repository
- Inside the repository folder: 'mkdir phescami-geant4/sources/build'
- 'mkdir phescami-geant4/sources/ROOT_files'
- 'cd phescami-geant4/sources/build'
- 'cmake ../src'
- 'make -j4'
- './exampleADHD' 
- The gui should appear on your terminal. To use a macro do: ./exampleADHD MyMacro.mac

### Event reconstruction
For event reconstruction use phescami-geant4/analysis/macros/ADHD_gateRange.cc for realistic reconstuction (prompt and delayed event are evaluated depending from the gate value). 
For event reconstruction based on MC truth use phescami-geant4/analysis/macros/ADHD_gateRange_MCTruth.cc.
The macros has some input-output file name hard coded that must be changed. As a suggestion, input files can be stored in sources/ROOT_files/ (you have to create the directory).


### Antimatter perfomances
If you have already simulated some pbar and dbar sample you can see some pbar-dbar discrimination perfomance. After running ADHD_gateRange.cc, type:

- 'root -l phescami-geant4/analysis/macros/IFAE_preparation.cc'
- 'root -l phescami-geant4/analysis/macros/IFAE_plot.cc'

Remember to always check file paths. 

### Ordinary matter background

To see matter studies and expected rates you must simulate protons, He, C, electrons. 
Than after reconstruction, run:

- 'root -l phescami-geant4/analysis/macros/trigger_plotter.cc'


## Support
francesco.rossi-9@unitn.it

## Contributing
Any contribution is well accepted. 

## Authors and acknowledgment
Francesco Rossi
Francesco Nozzoli

## License
For open source projects, say how it is licensed.

