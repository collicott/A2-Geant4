//--Author	JRM Annand   30th Sep 2003
//--Rev 	JRM Annand...15th Oct 2003 ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004 3v8 compatible
//--Rev 	JRM Annand...19th Feb 2004 User code
//--Rev 	JRM Annand...16th May 2005 ReadDecoded (bug G3 output)
//--Update	JRM Annand... 2nd Jun 2005 ReadRecoded bug fix (D.Glazier)
//--Update	D.Glazier...24th Aug 2007 Add time for ReadDecoded 
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2PlasticPID
//
// Internal Particle Identification Detector for the Crystal Ball
// Cylindrical array of plastic scintillators
//

#ifndef __TA2PlasticPID_h__
#define __TA2PlasticPID_h__

#include "MCBranchID.h"
#include "TA2Detector.h"

class TA2PlasticPID : public TA2Detector {
 private:
 public:
  TA2PlasticPID( const char*, TA2System* );// Normal use
  virtual ~TA2PlasticPID();
  //  virtual void SetConfig( char*, int );// decode and load setup info
  virtual void LoadVariable( );            // display/cut setup
  //  virtual void PostInit( );            // initialise using setup info
  virtual void Decode( );              // hits -> energy procedure
  virtual void SaveDecoded( );             // save local analysis
  virtual void ReadDecoded( ){             // read back previous analysis
    if(fName.Contains("PID")) PIDReadDecoded();
    if(fName.Contains("BaF2")) TAPSReadDecoded();
  }
  virtual void PIDReadDecoded( );
  virtual void TAPSReadDecoded( );

  ClassDef(TA2PlasticPID,1)
};

//---------------------------------------------------------------------------
inline void TA2PlasticPID::TAPSReadDecoded( )
{
  // Read back...
  //   either previously analysed data from Root Tree file
  //   or MC simulation results, assumed to have the same data structure
  // Bug fix remove fNhits-- line, (DG 27/5/05, implemented JRMA 2/6/05)

  UInt_t i,j;
  Double_t total = 0.0;
  fNhits = *(Int_t*)(fEvent[EI_nvtaps]);
  Float_t* energy = (Float_t*)(fEvent[EI_evtaps]);
  Int_t* index = (Int_t*)(fEvent[EI_ivtaps]);

  for( i=0; i<fNhits; i++ ){
    j = *index++;
    if( !j ){                     // is it s real hit
      //      fNhits--;
      i--;
      energy++;
      continue;
    }
    j--;
    fHits[i] = j;
    fEnergy[j] = (*energy++) * 1000.;
    fEnergyOR[i] = fEnergy[j];
    total += fEnergy[j];
  }
  fHits[i] = EBufferEnd;
  fEnergyOR[i] = EBufferEnd;
  fTotalEnergy = total;
}
inline void TA2PlasticPID::PIDReadDecoded( )
{
  // Read back...
  //   either previously analysed data from Root Tree file
  //   or MC simulation results, assumed to have the same data structure
  // Bug fix remove fNhits-- line, (DG 27/5/05, implemented JRMA 2/6/05)

  UInt_t i,j;
  Double_t total = 0.0;
  fNhits = *(Int_t*)(fEvent[EI_vhits]);
  Float_t* energy = (Float_t*)(fEvent[EI_eveto]);
  Float_t* time;
  if(fIsTime) time= (Float_t*)(fEvent[EI_tveto]);
  Int_t* index = (Int_t*)(fEvent[EI_iveto]);

  for( i=0; i<fNhits; i++ ){
    j = *index++;
    if( !j ){                     // is it s real hit
      //      fNhits--;
      i--;
      energy++;
      if(fIsTime)time++;
      continue;
    }
    j--;
    fHits[i] = j;
    fEnergy[j] = (*energy++) * 1000.;
    fEnergyOR[i] = fEnergy[j];
    if(fIsTime){
      fTime[j] = (*time++);
      fTimeOR[i] = fTime[j];
    }     
    total += fEnergy[j];
  }
  fHits[i] = EBufferEnd;
  fEnergyOR[i] = EBufferEnd;
  if(fIsTime)fTimeOR[i] = EBufferEnd;
  fTotalEnergy = total;
}

//---------------------------------------------------------------------------
inline void TA2PlasticPID::Decode( )
{
  // Run basic TA2Detector decode, then anything else

  DecodeBasic();
  // Anything else here           
}

#endif
