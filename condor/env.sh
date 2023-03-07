source /cvmfs/geant4.cern.ch/geant4/10.7/x86_64-centos7-gcc8-optdeb/bin/geant4.sh       # Geant 4 environment
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8/x86_64-centos7/setup.sh                     # gcc
source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/qt5/5.12.4/x86_64-centos7-gcc9-opt/Qt5-env.sh
#source /cvmfs/sft.cern.ch/lcg/releases/qt5/5.12.1-05305/x86_64-centos7-gcc8-opt/Qt5-env.sh 
#source /cvmfs/geant4.cern.ch/geant4/10.7/x86_64-centos7-gcc8-optdeb-MT/CMake-setup.sh   # CMake
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.18.3/Linux-x86_64/bin/:${PATH}       # CMake to be used
export QT_QPA_PLATFORM=offscreen
