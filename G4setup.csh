setenv G4WORKDIR ${HOME}/Geant4
source ${G4WORKDIR}/env.csh
setenv G4ROOTLIBS 1
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4INSTALL}/lib/${G4SYSTEM}:${G4WORKDIR}/lib/${G4SYSTEM}:${CLHEP_LIB_DIR}

