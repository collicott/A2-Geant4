//--Author	D Glazier   6th May 2008
//
// TOF apparatus .
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TOF
//
// User coded version of ToF apparatus.  Includes identification
// of charged particles through veto layer and further particle id
// from deltaE-E and tof analysis. Calculates the energy of nucleons
// from tof

#ifndef __TA2TOF_h__
#define __TA2TOF_h__

#include "TA2Apparatus.h"            // Acqu-Root histogrammer
#include "TA2LongScint.h"            // Acqu-Root histogrammer
#include <iostream>


class TA2TOF : public TA2Apparatus {
protected:
  TA2LongScint* fToFWall;                     // Glasgow photon tagger
  TA2LongScint* fToFVeto;
 
  Double_t fVetoThetaMax; //ToF hit and veto hit must be less than this angle apart
  Double_t fLayerThetaMax; //two tof bar hits must be less than this angle apart
                           //to be considered to be from the same particle

  Int_t fNLayers;  //number of layers of TOF bars, not inc. vetos 
  Double_t *fCoincVetoE;  //Energy of veto in coinc with tof (deltaE)
  Double_t *fCoincToFE;   //Energy of ToF hit in coinc with veto (E)
  Double_t *fCoincToFT;   //Time of ToF hit
  Double_t *fToFE;   //Energy of ToF hit (E)
  Double_t *fToFT;   //Time of ToF hit
  Double_t *fVetoTheta;  //look for correlated veto hits
  Double_t *fLayerTheta; //look for correlated layer hits
  Int_t     fNcharged;         // Number of tof clusters identified as charged

public:
  TA2TOF( const char*, TA2System* );  // pass ptr to analyser
  virtual ~TA2TOF();                  // destructor
  virtual void PostInit();                    // some setup after parms read in
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void LoadVariable(  );              // display setup
  virtual void Cleanup();                     // reset at end of event
  virtual void SetConfig( Char_t*, Int_t );   // setup decode in implement
  // switch to stardard or MC
  virtual void Reconstruct();
  // in your analyis you have to cast the "Taps-Apparatus" to get this functions

  ClassDef(TA2TOF,1)
};

//-----------------------------------------------------------------------------
inline void TA2TOF::Cleanup( )
{
  // Clear any arrays holding variables
 TA2DataManager::Cleanup();
 fCoincVetoE[0]=EBufferEnd;
 fCoincToFE[0]=EBufferEnd;
 fCoincToFT[0]=EBufferEnd;
 fToFE[0]=EBufferEnd;
 fToFT[0]=EBufferEnd;
 fVetoTheta[0]=EBufferEnd;
 fLayerTheta[0]=EBufferEnd;
 fNcharged=0;

 if( !fNparticle ) return;
 Int_t i;
 for( i=0; i<fNparticle; i++ ){
   fPDG_ID[i] = 0;
 }
 fP4tot.SetXYZT(0.,0.,0.,0.);
 fNparticle = 0;
}

//-----------------------------------------------------------------------------

#endif
