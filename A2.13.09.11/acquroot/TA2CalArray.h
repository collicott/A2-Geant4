//--Author	JRM Annand   27th Apr 2003
//--Rev 	JRM Annand...30th Sep 2003  New Detector bases
//--Rev 	JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004  3v8 compatible
//--Rev 	JRM Annand... 7th Jun 2005  ReadDecoded...use fEnergyScale
//--Update	JRM Annand...25th Oct 2005  ReadDecoded...energy thresh
//--Update	D.Glazier ...24th Aug 2007  ReadDecoded...include energy smearing
//--Update	D.Glazier ...24th Aug 2007  ReadDecoded...include detector time
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2CalArray
//
// Decoding and calibration methods for the Crystall Ball array of NaI(Tl)
// Configured as forward wall to work in conjunction with the CB
// This can use standard TA2ClusterDetector methods or ones redifined here
//

#ifndef __TA2CalArray_h__
#define __TA2CalArray_h__

#include "MCBranchID.h"
#include "TA2ClusterDetector.h"
#include "iostream.h"

class TA2CalArray : public TA2ClusterDetector {
 private:
 public:
  TA2CalArray( const char*, TA2System* );// Normal use
  virtual ~TA2CalArray();
  virtual void PostInit( );            // initialise using setup info
  virtual void ParseDisplay(char*);    // display setup
  virtual void Decode( );              // hits -> energy procedure
  virtual void Cleanup( );             // end-of-event cleanup
  virtual void SaveDecoded( );         // save local analysis
  virtual void ReadDecoded( );         // read back previous analysis

  ClassDef(TA2CalArray,1)
};

//---------------------------------------------------------------------------
inline void TA2CalArray::Decode( )
{
  // Decode the NaI ADC and TDC information
  // Decode raw TDC and Scaler information into
  // Hit pattern, "Energy" pattern, aligned OR etc.

  TA2ClusterDetector::Decode();
}

//---------------------------------------------------------------------------
inline void TA2CalArray::ReadDecoded( )
{
  // Read Crystal Ball energies from  GEANT-3 simulation output
  // See MCBranchID.h for GEANT-3 output details.
  // Add energy thresholds 25/10/05

  UInt_t i,k;
  Int_t j;

  fNhits = *(Int_t*)(fEvent[EI_nhits]);                // # crystals fired
  fTotalEnergy = *(Float_t*)(fEvent[EI_eNaI]) * 1000;  // total deposited
  Float_t* energy = (Float_t*)(fEvent[EI_ecryst]);     // energies in crystals
  Float_t* time;
  if(fIsTime)time = (Float_t*)(fEvent[EI_tcryst]);     // time in crystals
  Int_t* index = (Int_t*)(fEvent[EI_icryst]);          // crystal indices
  Double_t E,T;                                          // energy
  Double_t EthrLow,EthrHigh;                           // energy thresholds
  Double_t TthrLow,TthrHigh;                           // time thresholds
  Double_t Emc;
  Double_t res1=0.04;
  Double_t res2=0.00;
  for( i=0,k=0; i<fNhits; i++ ){                       // loop over hits
    // slightly less ugly decoding of icryst
    j = index[i] % 10000;                              // AcquRoot index
    if( j == -1 ) continue;                            // check spurious index
    EthrLow = fElement[j]->GetADCcut()->GetLowThr();   // low thresh
    EthrHigh = fElement[j]->GetADCcut()->GetHighThr(); // high thresh
    Emc=1000 * (*energy++);
    E =  fEnergyScale*(Emc+res1*fRandom->Gaus()*TMath::Power(Emc,0.75)+res2*fRandom->Gaus()*Emc);             // G3 output in GeV
    if( (E < EthrLow) || (E > EthrHigh) ) {*time++;continue;}    // energy inside thresh?
    fEnergy[j] = E;                                    // save energy
    fEnergyOR[k] = E;
    if(fIsTime){
      T=*time++;
      TthrLow = fElement[j]->GetTDCcut()->GetLowThr();   // low thresh
      TthrHigh = fElement[j]->GetTDCcut()->GetHighThr(); // high thresh
      if( (T < TthrLow) || (T > TthrHigh) ) continue;    // time inside thresh?
      fTime[j] = T;                                    // save time
      fTimeOR[k] = T;
    }
    fHits[k] = j;                                      // store hit
    k++;                                               // update # good hits
  }
  fHits[k] = EBufferEnd;                               // end markers
  fEnergyOR[k] = EBufferEnd;
  if(fIsTime)fTimeOR[k] = EBufferEnd;
  fNhits = k;                                          // # hits inside thresh
}

//---------------------------------------------------------------------------
inline void TA2CalArray::Cleanup( )
{
  // end-of-event cleanup
  TA2ClusterDetector::Cleanup();
}
#endif
