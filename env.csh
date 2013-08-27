######################################
#
# Clean all G4 envs
unsetenv  CLHEP_BASE_DIR
unsetenv  CLHEP_INCLUDE_DIR
unsetenv  CLHEP_LIB
unsetenv  CLHEP_LIB_DIR

unsetenv  G4DEBUG
unsetenv  G4INCLUDE
unsetenv  G4INSTALL

unsetenv  G4LEDATA
unsetenv  G4LEVELGAMMADATA
unsetenv  G4NEUTRONHPDATA
unsetenv  G4RADIOACTIVEDATA
unsetenv  G4ABLADATA
unsetenv  G4REALSURFACEDATA
unsetenv  G4NEUTRONXSDATA
unsetenv  G4PIIDATA

unsetenv  G4LIB
unsetenv  G4LIB_BUILD_G3TOG4
unsetenv  G4LIB_BUILD_SHARED
unsetenv  G4LIB_BUILD_STATIC
unsetenv  G4LIB_USE_DLL

unsetenv  G4LIB_BUILD_ZLIB
unsetenv  G4LIB_BUILD_GDML
unsetenv  G4LIB_USE_G3TOG4
unsetenv  G4LIB_USE_GRANULAR
unsetenv  G4LIB_USE_ZLIB

unsetenv  G4SYSTEM

unsetenv  G4UI_NONE
unsetenv  G4UI_BUILD_WIN32_SESSION
unsetenv  G4UI_BUILD_XAW_SESSION
unsetenv  G4UI_BUILD_XM_SESSION
unsetenv  G4UI_USE_TCSH
unsetenv  G4UI_USE_WIN32
unsetenv  G4UI_USE_XAW
unsetenv  G4UI_USE_XM
unsetenv  G4UI_USE_QT

unsetenv  G4VIS_NONE
unsetenv  G4VIS_BUILD_DAWN_DRIVER
unsetenv  G4VIS_BUILD_OIWIN32_DRIVER
unsetenv  G4VIS_BUILD_OIX_DRIVER
unsetenv  G4VIS_BUILD_OPENGLWIN32_DRIVER
unsetenv  G4VIS_BUILD_OPENGLXM_DRIVER
unsetenv  G4VIS_BUILD_OPENGLX_DRIVER
unsetenv  G4VIS_BUILD_RAYTRACERX_DRIVER
unsetenv  G4VIS_BUILD_VRML_DRIVER
unsetenv  G4VIS_BUILD_OPENGLQT_DRIVER

unsetenv  G4VIS_USE_DAWN
unsetenv  G4VIS_USE_OIWIN32
unsetenv  G4VIS_USE_OIX
unsetenv  G4VIS_USE_OPENGLWIN32
unsetenv  G4VIS_USE_OPENGLX
unsetenv  G4VIS_USE_OPENGLXM
unsetenv  G4VIS_USE_RAYTRACERX
unsetenv  G4VIS_USE_VRML
unsetenv  G4VIS_USE_OPENGLQT

######################################
#
# g4system.U
#
#+
setenv G4SYSTEM "Linux-g++"

#
# g4dirs.U
#
#+
if ( X/cern/geant4.9.4.p04 != X/cern/geant4 ) then
    setenv G4INSTALL "/cern/geant4/src/geant4"
else
    setenv G4INSTALL "/cern/geant4.9.4.p04"
endif

#+
if ( Xn != Xn ) then
    if ( X/cern/geant4.9.4.p04 != X/cern/geant4 ) then
        setenv G4INCLUDE "/cern/geant4/include/geant4"
    else
        setenv G4INCLUDE "/cern/geant4/include"
    endif

endif

#+
if ( X/cern/geant4.9.4.p04/lib != X ) then
    if ( X/cern/geant4.9.4.p04 != X/cern/geant4 ) then
        setenv G4LIB "/cern/geant4/lib/geant4"
    else
        setenv G4LIB "/cern/geant4/lib"
    endif

endif

#+
if ( X/cern/geant4.9.4.p04/data/PhotonEvaporation2.1 != X ) then
    setenv G4LEVELGAMMADATA "/cern/geant4.9.4.p04/data/PhotonEvaporation2.1"
endif

#+
if ( X/cern/geant4.9.4.p04/data/RadioactiveDecay3.3 != X ) then
    setenv G4RADIOACTIVEDATA "/cern/geant4.9.4.p04/data/RadioactiveDecay3.3"
endif

#+
if ( X/cern/geant4.9.4.p04/data/G4EMLOW6.19 != X ) then
    setenv G4LEDATA "/cern/geant4.9.4.p04/data/G4EMLOW6.19"
endif

#+
if ( X/cern/geant4.9.4.p04/data/G4NDL3.14 != X ) then
    setenv G4NEUTRONHPDATA "/cern/geant4.9.4.p04/data/G4NDL3.14"
endif

#+
if ( X/cern/geant4.9.4.p04/data/G4ABLA3.0 != X ) then
    setenv G4ABLADATA "/cern/geant4.9.4.p04/data/G4ABLA3.0"
endif

#+
if ( X/cern/geant4.9.4.p04/data/RealSurface1.0 != X ) then
    setenv G4REALSURFACEDATA "/cern/geant4.9.4.p04/data/RealSurface1.0"
endif

#+
if ( X/cern/geant4.9.4.p04/data/G4NEUTRONXS1.0 != X ) then
    setenv G4NEUTRONXSDATA "/cern/geant4.9.4.p04/data/G4NEUTRONXS1.0"
endif

#+
if ( X/cern/geant4.9.4.p04/data/G4PII1.2 != X ) then
    setenv G4PIIDATA "/cern/geant4.9.4.p04/data/G4PII1.2"
endif



#
# g4clhep.U
#
if ( X/cern/clhep != X ) then
    setenv CLHEP_BASE_DIR "/cern/clhep"
endif

#+
if ( X/cern/clhep/include != X ) then
    setenv CLHEP_INCLUDE_DIR "/cern/clhep/include"
endif

#+
if ( X/cern/clhep/lib != X ) then
    setenv CLHEP_LIB_DIR "/cern/clhep/lib"
endif

#+
if ( XCLHEP != X ) then
    setenv CLHEP_LIB "CLHEP"
endif

#+
#
# g4debug
#
if ( Xy == Xy ) then
    setenv G4DEBUG 1
endif


#
# g4ui
#
#+
if ( Xn == Xy ) then
    setenv G4UI_NONE 1
endif

# Check for Windows!
if ( "X$G4SYSTEM" != "XWIN32-VC" && "X$G4SYSTEM" != "XWIN32-VC7" ) then
    if ( Xn != Xy ) then
        setenv G4UI_USE_TCSH 1
    endif
endif

#+
if ( Xy == Xy ) then
    setenv G4UI_BUILD_XAW_SESSION 1
endif

#+
if ( Xy == Xy ) then
    setenv G4UI_USE_XAW 1
endif

#+
if ( Xy == Xy ) then
    setenv G4UI_BUILD_XM_SESSION 1
endif

#+
if ( Xy == Xy ) then
    setenv G4UI_USE_XM 1
endif

#+
if ( Xn == Xy ) then
    setenv G4UI_BUILD_WIN32_SESSION 1
endif

#+
if ( Xn == Xy ) then
    setenv G4UI_USE_WIN32 1
endif

#+
if ( Xy == Xy ) then
    setenv G4UI_BUILD_QT_SESSION 1
endif 

#+
if ( Xy == Xy ) then
    setenv G4UI_USE_QT 1
endif 



#
# g4vis
#
#+
if ( Xn == Xy ) then
    setenv G4VIS_NONE 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_DAWN_DRIVER 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_OPENGLX_DRIVER 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_BUILD_OPENGLWIN32_DRIVER 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_BUILD_OIX_DRIVER 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_RAYTRACERX_DRIVER 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_BUILD_OIWIN32_DRIVER 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_VRML_DRIVER 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_BUILD_OPENGLQT_DRIVER 1
endif 


#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_DAWN 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_OPENGLX 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_OPENGLXM 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_USE_OPENGLWIN32 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_USE_OIX 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_RAYTRACERX 1
endif

#+
if ( Xn == Xy ) then
    setenv G4VIS_USE_OIWIN32 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_VRML 1
endif

#+
if ( Xy == Xy ) then
    setenv G4VIS_USE_OPENGLQT 1
endif

#+
if ( X != X )  then
    setenv OGLHOME ""
endif 

#+
if ( X != X )  then
    setenv OIVHOME ""
endif 


#+
if ( Xy != X )  then
    setenv XMFLAGS ""
endif 

#+
if ( Xy != X )  then
    setenv XMLIBS " -lXm -lXpm "
endif 

#+
if ( Xy != X )  then
    setenv XMFLAGS ""
endif 

#+
if ( Xy != X )  then
    setenv XMLIBS " -lXm -lXpm "
endif 

#+
if ( Xy != X )  then
    setenv XAWFLAGS ""
endif 

#+
if ( Xy != X )  then
    setenv XAWLIBS " -lXaw "
endif 


#
# Qt Flags and Libs, messy, but needed for backward compatibility
#+
if ( "Xy" == "Xy" || "Xy" == "Xy" )  then
    setenv QTFLAGS "-I/usr/include -I/usr/include/QtCore -I/usr/include/QtGui"
    setenv QTLIBS "-L/usr/lib64 -lQtCore_debug -lQtGui"
    setenv QTMOC "/usr/bin/moc-qt4"
endif

if ( "Xy" == "Xy" || "Xy" == "Xy" )  then
    if ( "X$QTFLAGS" == "X" )  then
        setenv QTFLAGS "-I/usr/include -I/usr/include/QtCore -I/usr/include/QtGui  -I/usr/include/QtOpenGL"
    else
        setenv QTFLAGS "$QTFLAGS  -I/usr/include/QtOpenGL"
    endif

    if ( "X$QTMOC" == "X" )  then
        setenv QTMOC "/usr/bin/moc-qt4"
    endif

    setenv GLQTLIBS "-L/usr/lib64 -lQtCore_debug -lQtGui -lQtOpenGL"
endif


#
# Use GDML module
#
#+
if ( Xn == Xy ) then
    setenv G4LIB_BUILD_GDML 1
endif 

if ( Xn == Xy ) then
    setenv XERCESCROOT ""
endif

if ( Xn == Xy )  then
    setenv G4LIB_BUILD_G3TOG4 1
endif 

if ( Xn == Xy )  then
    setenv G4LIB_USE_G3TOG4 1
endif 


if ( Xn == Xy )  then
    setenv G4LIB_BUILD_ZLIB 1
endif 

if ( X == Xy )  then
    setenv G4LIB_USE_ZLIB 1
endif 

#+
#
# g4shared
#
if ( Xy == Xy ) then
    setenv G4LIB_BUILD_SHARED 1
endif

if ( Xn == Xy ) then
    setenv G4LIB_BUILD_STATIC 1
endif


if ( Xn == Xy ) then
    setenv G4LIB_USE_DLL 1
endif

#+
#
# g4granular
#
if ( Xn == Xy ) then
    setenv G4LIB_USE_GRANULAR 1
endif

#####################################################################



#+
#
# G4WORKDIR
#
if ( ${?G4WORKDIR} ) then
    echo " "
else
    # Check for Windows!
    if ( "X$G4SYSTEM" == "XWIN32-VC" || "X$G4SYSTEM" == "XWIN32-VC7" ) then
        setenv G4WORKDIR "c:/geant4"
    else # if Unix
        setenv G4WORKDIR $HOME/geant4
    endif
endif

#
# *NIX Shared Libraries
# If we built Geant4 with shared libraries, we need to add the Gean4
# library directory to (DY)LD_LIBRARY_PATH.
# In all cases, external shared library directories should be added to
# (DY)LD_LIBRARY_PATH

if ( "X$G4SYSTEM" != "XDarwin-g++" ) then
    if ( ${?LD_LIBRARY_PATH} ) then
        if ( ${?G4LIB_BUILD_SHARED} ) then
            setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}	
        endif

        setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CLHEP_LIB_DIR}

        if ( ${?G4LIB_BUILD_GDML} ) then
            setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${XERCESCROOT}/lib
        endif

    else
        if ( ${?G4LIB_BUILD_SHARED} ) then
            setenv LD_LIBRARY_PATH ${G4LIB}/${G4SYSTEM}
            setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CLHEP_LIB_DIR}
        else
            setenv LD_LIBRARY_PATH ${CLHEP_LIB_DIR}
        endif

        if ( ${?G4LIB_BUILD_GDML} ) then
            setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${XERCESCROOT}/lib
        endif
    endif
else
    #
    # Darwin Shared Libraries
    # we repeat the above logic, but for DYLD_LIBRARY_PATH
    #   
    if ( ${?DYLD_LIBRARY_PATH} ) then
        if ( ${?G4LIB_BUILD_SHARED} ) then
            setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}	
        endif

        setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${CLHEP_LIB_DIR}

        if ( ${?G4LIB_BUILD_GDML} ) then
            setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${XERCESCROOT}/lib
        endif
    else
        if ( ${?G4LIB_BUILD_SHARED} ) then
            setenv DYLD_LIBRARY_PATH ${G4LIB}/${G4SYSTEM}
            setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${CLHEP_LIB_DIR}
        else
            setenv DYLD_LIBRARY_PATH ${CLHEP_LIB_DIR}
        endif

        if ( ${?G4LIB_BUILD_GDML} ) then
            setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${XERCESCROOT}/lib
        endif
    endif
endif


setenv PATH ${PATH}:${G4WORKDIR}/bin/${G4SYSTEM}

