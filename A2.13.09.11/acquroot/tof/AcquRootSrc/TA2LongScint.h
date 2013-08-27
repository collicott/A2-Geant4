//--Author	JRM Annand   27th Apr 2003
//--Rev 	JRM Annand...30th Sep 2003  New Detector bases
//--Rev 	JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004  3v8 compatible
//--Rev 	JRM Annand...25th Jul 2005  SetConfig hierarchy
//--Update	JRM Annand...18th Jan 2007  Inherit from TA2Detector
//--Update	D Glazier... 6th May 2008   Implement ReadDecoded
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2LongScint
//
// Decoding and calibration methods for a long bar of plastic scintillator
// Assumes one PMT at each end of bar -> position sensitivity, mean time etc.
//

#ifndef __TA2LongScint_h__
#define __TA2LongScint_h__

#include "MCBranchID.h"
#include "TA2Detector.h"
#include "LongBar_t.h"
#include <iostream.h>

class TA2LongScint : public TA2Detector {
 private:
 protected:
  LongBar_t** fBar;                    // Array of scintillator bar info
  Double_t* fMeanEnergy;               // Mean pulse heights
  Double_t* fMeanEnergyOR;             // Mean pulse-height OR
  Double_t* fMeanTime;                 // Mean times
  Double_t* fMeanTimeOR;               // Mean-time OR
  Double_t* fTimeDiff;                 // Time difference
  Double_t* fTimeDiffOR;               // Time difference OR
  TVector3* fBarPos;                   // Hit position on Bar
  Int_t* fBarHits;                     // Scint-bar hits
  Int_t* fBarHitCount;                 // # elements fired for each bar
  Int_t* fElement2Bar;                 // element-bar correspondence
  Int_t fNBarHits;                     // # scint-bar hits
  Int_t fNBar;                         // # scintillator bars
  Int_t fNbar;                         // # scintillator bars running counter
  Int_t fNVetos;      //#veto detector associated with wall (=-1 if it is a veto detector for ReadDecoded dglazier
 public:
  TA2LongScint( const char*, TA2System* );// Normal use
  virtual ~TA2LongScint();
  virtual void PostInit( );            // initialise using setup info
  virtual void ParseDisplay(char*);    // display setup
  virtual void SetConfig(Char_t*, Int_t);
  virtual void LoadVariable();
  virtual void Decode( );              // hits -> energy procedure
  virtual void DecodeBar();            // combine hits to bar
  virtual void Cleanup( );             // end-of-event cleanup
  virtual void SaveDecoded( );         // save local analysis
  virtual void ReadDecoded( );         // read back previous analysis
  virtual void TimeCorrect( Double_t );// apply start time correction
  //
  LongBar_t** GetBar(){ return fBar; }
  LongBar_t* GetBar(Int_t i){ return fBar[i]; };
  Double_t* GetMeanEnergy(){ return fMeanEnergy; }
  Double_t GetMeanEnergy(Int_t i){ return fMeanEnergy[i]; }
  Double_t* GetMeanEnergyOR(){ return fMeanEnergyOR; }
  Double_t* GetMeanTime(){ return fMeanTime; }
  Double_t GetMeanTime(Int_t i){ return fMeanTime[i]; }
  Double_t* GetMeanTimeOR(){ return fMeanTimeOR; }
  Double_t* GetTimeDiff(){ return fTimeDiff; }
  Double_t GetTimeDiff(Int_t i){ return fTimeDiff[i]; }
  Double_t* GetTimeDiffOR(){ return fTimeDiffOR; }
  TVector3* GetBarPos(){ return fBarPos; }
  TVector3 GetBarPos(Int_t i){ return fBarPos[i]; }
  Int_t* GetBarHits(){ return fBarHits; }
  Int_t GetNBarHits(){ return fNBarHits; }
  Int_t GetNBar(){ return fNBar; }
  Int_t GetNbar(){ return fNbar; }

  ClassDef(TA2LongScint,1)
};

//---------------------------------------------------------------------------
inline void TA2LongScint::Decode( )
{
  // Decode the NaI ADC and TDC information
  // Decode raw TDC and Scaler information into
  // Hit pattern, "Energy" pattern, aligned OR etc.

  TA2Detector::Decode();
  DecodeBar();
}

//---------------------------------------------------------------------------
inline void TA2LongScint::DecodeBar()
{
  // Step through PMT hits array
  // back reference PMTs to bars and require 2 PMT hits for a bar hit
  // Calculate bar quantities....mean time, position etc.
  Int_t j;
  fNBarHits = 0;
  for( UInt_t i=0; i<fNhits; i++ ){   // step through decoded hits
    j = fElement2Bar[ fHits[i] ];     // ref bar to PMT
    fBarHitCount[j]++;
    if( fBarHitCount[j] == 2 ){       // require both PMTs fire
      fBarHits[fNBarHits] = j;        // record bar hit
      fBar[j]->TMean();               // mean time
      fBar[j]->GMean();               // mean energy
      fBar[j]->TDiff();               // time difference
      fBar[j]->Pos();                 // hit position
      // Save OR quanties for diagnostic spectra
      fMeanEnergyOR[fNBarHits] = fMeanEnergy[j];
      fMeanTimeOR[fNBarHits] = fMeanTime[j];
      fTimeDiffOR[fNBarHits] = fTimeDiff[j];
      fNBarHits++;                    // hit counter
    }
  }
  // mark end of variable length arrays
  fBarHits[fNBarHits] = EBufferEnd;
  fMeanEnergyOR[fNBarHits] = EBufferEnd;
  fMeanTimeOR[fNBarHits] = EBufferEnd;
  fTimeDiffOR[fNBarHits] = EBufferEnd;
}

//---------------------------------------------------------------------------
inline void TA2LongScint::ReadDecoded( )
{
  Double_t total = 0.0;
  Int_t Nhits = *(Int_t*)(fEvent[EI_ntof]);
  Float_t* energy = (Float_t*)(fEvent[EI_tofe]);
  Float_t* time = (Float_t*)(fEvent[EI_toft]);
  Float_t* xpos = (Float_t*)(fEvent[EI_tofx]);
  Float_t* ypos = (Float_t*)(fEvent[EI_tofy]);
  Int_t* index = (Int_t*)(fEvent[EI_tofi]);

  Int_t i,j;
  TVector3 hpos;
  TVector3 align;
  fNBarHits=0;
  for( i=0; i<Nhits; i++ ){
    j = *index++;
    //Veto detectors will take the first fNVetos indexes from G4
    if(fNVetos==-1){//it is veto detector
      if(j>=fNBar) {*energy++;*time++;continue;}//this bar is main tof, not veto
    }
    else if(j<fNVetos){*energy++;*time++;continue;}//this bar from G4 is veto, but this is not veto detector
    else j-=fNVetos;//this is a main tof det and this is a tof bar from sim

    if(j>=fNBar){*energy++;*time++;continue;}
    fBarHits[fNBarHits] = j;
    fMeanEnergy[j] = (*energy++) * 1000.;
    fMeanEnergyOR[fNBarHits] = fMeanEnergy[j];
    total += fMeanEnergy[j];
    fMeanTime[j] = (*time++) ;
    fMeanTimeOR[fNBarHits] = fMeanTime[j];
    hpos.SetXYZ(*xpos++,*ypos++,0);
    hpos-=fBar[j]->GetMeanPos();
    align=fBar[j]->GetAlignment();
    hpos.SetXYZ(hpos.x()*align.x(),hpos.y()*align.y(),hpos.z()*align.z()); //position of hit along bar
    fTimeDiff[j]=(1./fBar[j]->GetCeff()*hpos).X()+(1./fBar[j]->GetCeff()*hpos).Y();//convert pos to time
    fTimeDiffOR[fNBarHits]=fTimeDiff[j];
    fBar[j]->Pos();//Calculate position from TDiff 
    fNBarHits++;
   }
  fBarHits[fNBarHits] = EBufferEnd;
  fMeanEnergyOR[fNBarHits] = EBufferEnd;
  fMeanTimeOR[fNBarHits] = EBufferEnd;
  fTimeDiffOR[fNBarHits] = EBufferEnd;
  fTotalEnergy = total;
}

//---------------------------------------------------------------------------
inline void TA2LongScint::TimeCorrect( Double_t toff )
{
  // Apply correction to time arrays...
  // e.g. correction to common TDC start/stop time
  for(Int_t j=0; j<fNBarHits; j++) fMeanTime[ fBarHits[j] ] -= toff;
  for(UInt_t j=0; j<fNhits; j++) fTime[ fHits[j] ] -= toff;
}

//---------------------------------------------------------------------------
inline void TA2LongScint::Cleanup( )
{
  // end-of-event cleanup
  TA2Detector::Cleanup();
  Int_t j;
  for(Int_t i=0; i<fNBarHits; i++){
    j = fBarHits[i];
    fBarHitCount[j] = 0;
    fBar[j]->Cleanup();
  }
  fNBarHits=0;
}


#endif
