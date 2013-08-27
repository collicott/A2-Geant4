//--Author      V Lisin      28th Jun 2004  original DAPHNE fortran -> C
//--Rev         JRM Annand... 8th Jul 2004  incorporate AcquRoot C++ class
//--Update      JRM Annand... 7th Sep 2004  alpha-version for online
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2CylMWPC
//
// Decoding and calibration methods for cylindrical DAPHNE pattern
// multi-wire chamber charged-particle trackers for the Crystall Ball.
// 2 or 3 Cylindrical MWPC (aligned along the photon beam axis)
// with helically-wound cathode strip readout.
// Wires are read into TDCs or pattern units, strips into ADCs
//

#ifndef __TA2CylMWPC_h__
#define __TA2CylMWPC_h__

#include "TRandom3.h"

#include "MCBranchID.h"
#include <iostream>
#include "TA2WireChamber.h"
#include "TA2CylStripSven.h"
#include "TA2CylWireSven.h"
#include "TA2WCLayerSven.h"
#include "TA2Apparatus.h"
#include "TA2Calorimeter.h"
#include "TH2.h"
#include "TF1.h"
using namespace std;
/*#define TWOPI   6.283185307
#define PI      3.141592654
#define RADDEG 57.295779515*/

#define PI      TMath::Pi()
#define TWOPI   TMath::TwoPi()
#define RADDEG  TMath::RadToDeg()
#define DEGRAD  TMath::DegToRad()

class TA2CylMWPC : public TA2WireChamber {
protected:
  Double_t* fR;                   // chamber wire radii
  Double_t* fRtE;                 // chamber external strip radii
  Double_t* fRtI;                 // chamber inner strip radii
  Double_t* fC1;
  Double_t* fC2;
  Double_t* fD1;
  Double_t* fD2;
  Double_t* fNextPhiInt;
  Double_t fDPhiMax;              // max difference in wire-strip phi
  Double_t fDPhi12Max;            // max difference in phi chambers 1,2
  Double_t fMaxVertRadius;
  Double_t fMaxVertZ;
  Double_t fMinVertZ;
  Double_t** fZIntersect;
  Double_t** fPhiSWDiff;
  Double_t** fPhiWire;
  Double_t* fPhiDiff;
  Int_t* fNIntersect;
  Double_t* fVertX;
  Double_t* fVertY;
  Double_t* fVertZ;
  Double_t* fVertR;
  TVector3** fPseudoVertex;
  Double_t* fPseudoVertX;
  Double_t* fPseudoVertY;
  Double_t* fPseudoVertZ;
  Double_t* fPseudoVertR;
  Double_t* fDPhi12;
  Double_t fSTesting;     //Spare single-hit spectrum for testing purposes

  Double_t fClustEnI[2][256];
  Double_t fClustEnE[2][256];
  Int_t DeadStrPatt[4][128];
  Double_t TargZShift;

  Bool_t Simulate;
  char FileSigmaTheta[256];
  char FileSigmaPhi[256];
  char FileEfficiency[256];
  Double_t SigmaTheta[36];
  Double_t SigmaPhi[36];
  Double_t Efficiency[36];
  TRandom3* pRandoms;
  TH2* fhZtheta;
  TF1* fZres;
  //----------------------------------------------------------------------------
  //This is Sven's Debug Code. Don't use it!
  /*Int_t fA1;
  Int_t fA2;
  Int_t fB1;
  Int_t fB2;*/
  //----------------------------------------------------------------------------
  TA2Apparatus* fCB;
  TA2Detector* fPID;
  Double_t fPhiPlastic[24];
  Int_t StripWire(Int_t, Double_t, Double_t, Double_t[], Double_t*);

  //Addition dglazier for MC redecoded
  TVector3 fLayerv[100][4]; //put the cathode hit vectors 4=number of layers
 
public:
  TA2CylMWPC( const char*, TA2System* );// Normal use
  virtual ~TA2CylMWPC();
  virtual void DeleteArrays();         // flush local new store
  virtual void SaveDecoded( );         // save local analysis
  virtual void ReadDecoded( );         // read back previous analysis
  virtual void PostInit();             // some pre-event setup
  virtual void LoadVariable();
  virtual void SetChamberParm(Int_t, Double_t*);// save chamber param
  virtual void Decode();               // decode event info
  virtual void IntersectLayers( Int_t );
  virtual void MakeTracks();
  virtual void SetConfig( char*, Int_t );        // decode and load setup info
  virtual void ParseMisc(char* line);
  virtual void ReadParamFiles();
  virtual Double_t AbsDiffPhiSven(Double_t, Double_t);
  Double_t GetDPhiMax(){ return fDPhiMax; }
  Double_t** GetZIntersect(){ return fZIntersect; }
  Double_t* GetZIntersect(Int_t i){ return fZIntersect[i]; }
  Double_t** GetPhiSWDiff(){ return fPhiSWDiff; }
  Double_t* GetPhiSWDiff(Int_t i){ return fPhiSWDiff[i]; }
  Double_t** GetPhiWire(){ return fPhiWire; }
  Double_t* GetPhiWire(Int_t i){ return fPhiWire[i]; }
  Double_t* GetPhiDiff(){ return fPhiDiff; }
  Double_t* GetVertX(){ return fVertX; }
  Double_t* GetVertY(){ return fVertY; }
  Double_t* GetVertZ(){ return fVertZ; }
  Double_t* GetPseudoVertX(){ return fPseudoVertX; }
  Double_t* GetPseudoVertY(){ return fPseudoVertY; }
  Double_t* GetPseudoVertZ(){ return fPseudoVertZ; }
  Double_t GetClustEnI(Int_t nChamber, Int_t nCluster){return fClustEnI[nChamber][nCluster];}
  Double_t GetClustEnE(Int_t nChamber, Int_t nCluster){return fClustEnE[nChamber][nCluster];}
  Double_t* GetClustEnIPtr(Int_t nChamber){return fClustEnI[nChamber];}
  Double_t* GetClustEnEPtr(Int_t nChamber){return fClustEnE[nChamber];}
  Int_t GetNTrack(){ return fNTrack; }
  Int_t* GetNIntersect(){ return fNIntersect; }
  Int_t GetNIntersect(Int_t i){ return fNIntersect[i]; }

  ClassDef(TA2CylMWPC,1)
};

//---------------------------------------------------------------------------
inline void TA2CylMWPC::ReadDecoded(){
 
  if(!fhZtheta){ gROOT->cd();fhZtheta=new TH2F("DZth","DZth",100,0,180,100,-15,15);}
  if(!fZres){ gROOT->cd();fZres=new TF1("Zres","-0.6+1/sin(x)");}
  fNhits = *(Int_t*)(fEvent[EI_nmwpc]);  //total number of hits in all layers
  Int_t* Index=(Int_t*)(fEvent[EI_imwpc]); //index of layer hit
  Float_t *posx=(Float_t*)(fEvent[EI_mposx]); 
  Float_t *posy=(Float_t*)(fEvent[EI_mposy]); 
  Float_t *posz=(Float_t*)(fEvent[EI_mposz]); 

  Int_t LayNhits[]={0,0,0,0};
  //first read in all hit vectors
  for(Int_t i=0;i<fNhits;i++){
    //    std::cout<<" TA2CylMWPC::ReadDecoded() "<<std::endl<<Index[i]<<" "<<LayNhits[Index[i]]<<" "<<posx[i]<<" "<<posy[i]<<" "<<posz[i]<<std::endl;
    fLayerv[LayNhits[Index[i]]][Index[i]].SetXYZ(posx[i],posy[i],posz[i]);
    //    std::cout<<Index[i]<<" "<<fLayerv[LayNhits[Index[i]]][Index[i]].X()<<" "<<fLayerv[LayNhits[Index[i]]][Index[i]].Y()<<" "<<fLayerv[LayNhits[Index[i]]][Index[i]].Z()<<std::endl;
    LayNhits[Index[i]]++;
  }
  //Now look for intersections between cathode planes
  fNIntersect[0]=fNIntersect[1]=0;
  for(Int_t i=0;i<LayNhits[0];i++)//loop over chamber1 inner 
    for(Int_t j=0;j<LayNhits[1];j++)//loop over chamber1 outer
      if(fLayerv[i][0].Angle(fLayerv[j][1])<fDPhiMax){//use the real data defined dPhi (not quite the same thing though!!) 
	//we have an intersection
	fZIntersect[0][fNIntersect[0]]=(fLayerv[i][0].Z()+fLayerv[j][1].Z())/2+fZres->Eval((fLayerv[i][1]-fLayerv[i][0]).Theta())/2.35*pRandoms->Gaus();
	fhZtheta->Fill((fLayerv[i][1]-fLayerv[i][0]).Theta()*TMath::RadToDeg(),fZIntersect[0][fNIntersect[0]]-(fLayerv[i][0].Z()+fLayerv[j][1].Z())/2);
	fPhiWire[0][fNIntersect[0]]=(fLayerv[i][0]+fLayerv[j][1]).Phi()+PI/fWCLayer[fChamberLayers[0][2]]->GetNElement()*(2*pRandoms->Uniform()-1);
	//	std::cout<<fZIntersect[0][fNIntersect[0]]<<" "<<fPhiWire[0][fNIntersect[0]]*TMath::RadToDeg()<<" "<<fDPhiMax<<" "<<fLayerv[i][0].Perp()<<" "<<fLayerv[j][1].Perp()<<std::endl;
	fNIntersect[0]++;
      }
  for(Int_t i=0;i<LayNhits[2];i++)//loop over chamber2 inner 
    for(Int_t j=0;j<LayNhits[3];j++)//loop over chamber2 outer
      if(fLayerv[i][2].Angle(fLayerv[j][3])<fDPhiMax){//use the real data defined dPhi (not quite the same thing though!!) 
	//we have an intersection
	fZIntersect[1][fNIntersect[1]]=(fLayerv[i][2].Z()+fLayerv[j][3].Z())/2+fZres->Eval((fLayerv[i][1]-fLayerv[i][0]).Theta())/2.35*pRandoms->Gaus();
	fPhiWire[1][fNIntersect[1]]=(fLayerv[i][2]+fLayerv[j][3]).Phi()+PI/fWCLayer[fChamberLayers[1][2]]->GetNElement()*(2*pRandoms->Uniform()-1);
	//      	std::cout<<fZIntersect[1][fNIntersect[1]]<<" "<<fPhiWire[1][fNIntersect[1]]*TMath::RadToDeg()<<" "<<fDPhiMax<<" "<<fLayerv[i][2].Perp()<<" "<<fLayerv[j][3].Perp()<<std::endl;
	fNIntersect[1]++;
      }

  MakeTracks();
}
/* inline void TA2CylMWPC::ReadDecoded() */
/* { */
/*   // Read back... */
/*   //   either previously analysed data from Root Tree file */
/*   //   or MC simulation results, assumed to have the same data structure */
/*   // No MC simulation of MWPC yet */

/*   Float_t* dircos = (Float_t*)(fEvent[EI_dircos]); */
/*   Float_t* vertex = (Float_t*)(fEvent[EI_vertex]); */
/*   Int_t npart = *(Int_t*)(fEvent[EI_npart]); */
/*   Int_t* idpart = (Int_t*)(fEvent[EI_idpart]); */
/*   Int_t proton; */

/*   //Find proton index in simulated data */
/*   for(proton=0; proton<npart; proton++) */
/*     if(idpart[proton]==14) break; //14 = GEANT proton */

/*   //Move pointer to correct positions for proton data */
/*   for(Int_t t=0; t<proton; t++) //For all particles before proton: */
/*     for(Int_t i=0;i<3; i++) dircos++; //Jump over x,y,z components */

/*   //Read x,y,z components for proton */
/*   Double_t x = *dircos++; */
/*   Double_t y = *dircos++; */
/*   Double_t z = *dircos; */

/*   //Calculate theta and phi angles out of cartesian components */
/*   Double_t theta = acos(z); */
/*   Double_t phi = atan2(y, x); */

/*   //Get theta bin */
/*   Double_t fTheta = (theta*RADDEG)/5.0; */
/*   Int_t nTheta = (Int_t)fTheta;; */

/*   //Produce theta and phi errors for smearing */
/*   Double_t dtheta = pRandoms->Gaus(0.0, SigmaTheta[nTheta]); */
/*   Double_t dphi = pRandoms->Gaus(0.0, SigmaPhi[nTheta]); */

/*   //Calculate trigonometrics of smeared angles */
/*   Double_t costheta = TMath::Cos(theta + dtheta); */
/*   Double_t sintheta = TMath::Sin(theta + dtheta); */
/*   Double_t cosphi = TMath::Cos(phi + dphi); */
/*   Double_t sinphi = TMath::Sin(phi + dphi); */

/*   //Do only build tracks within given efficiency */
/*   if((pRandoms->Rndm() < Efficiency[nTheta]) && Simulate) */
/*   { */
/*     //Build track with angles smeared out within given resolutions */
/*     fNTrack = 1; */
/*     fTrackTheta[0] = theta + dtheta; */
/*     fTrackPhi[0] = phi + dphi; */
/*     fTrack[0]->SetTrack(0.0, 0.0, 0.0, cosphi*sintheta, sinphi*sintheta, costheta); */

/*     //Pseudo-vertex is still exact. Should be smeared out also... */
/*     fPseudoVertX[0] = vertex[0] * 10.0; //cm to mm */
/*     fPseudoVertY[0] = vertex[1] * 10.0; //cm to mm */
/*     fPseudoVertZ[0] = vertex[2] * 10.0; //cm to mm */
/*     fPseudoVertR[0] = TMath::Sqrt(fPseudoVertX[0]*fPseudoVertX[0] + fPseudoVertY[0]*fPseudoVertY[0] + fPseudoVertZ[0]*fPseudoVertZ[0]); */
/*   } */
/*   else */
/*     fNTrack = 0; */

/*   //Terminate */
/*   fTrackTheta[fNTrack] = EBufferEnd; */
/*   fTrackPhi[fNTrack] = EBufferEnd; */
/*   fPseudoVertX[fNTrack] = EBufferEnd; */
/*   fPseudoVertY[fNTrack] = EBufferEnd; */
/*   fPseudoVertZ[fNTrack] = EBufferEnd; */
/*   fPseudoVertR[fNTrack] = EBufferEnd; */

/*   //Kill some unused stuff... */
/*   for(Int_t i=0; i<fNLayer; i++) */
/*     (fWCLayer[i]->GetLayerHits())[0] = EBufferEnd; */
/*   for(Int_t n=0; n<fNchamber; n++) */
/*   { */
/*    fPhiSWDiff[n][0] = EBufferEnd; */
/*    fZIntersect[n][0] = EBufferEnd; */
/*    fPhiSWDiff[n][0] = EBufferEnd; */
/*    fPhiWire[n][0] = EBufferEnd; */
/*   } */
/*   fVertX[0] = EBufferEnd; */
/*   fVertY[0] = EBufferEnd; */
/*   fVertZ[0] = EBufferEnd; */
/*   fVertR[0] = EBufferEnd; */
/*   fPhiDiff[0] = EBufferEnd; */
/*   fDPhi12[0] = EBufferEnd; */
/* } */

//---------------------------------------------------------------------------

inline Double_t TA2CylMWPC::AbsDiffPhiSven(Double_t phi1, Double_t phi2)
{
  // Utility to provide absolute difference in phi angle
  // always returns +ve number

  Double_t diff = fabs(phi1-phi2);
  if(diff>PI) diff = TWOPI - diff;

  return fabs(diff);
}

//---------------------------------------------------------------------------

inline void TA2CylMWPC::Decode()
{
  // Do the basic decoding (e.g. ADC -> Energy)
  // Find chamber coordinates from clusters of hits
  // Join coordinates to make tracks

 //----------------------------------------------------------------------------
  //This is Sven's Debug Code. Don't use it!
  /*static int Count = 0;
  static int Cnt0 = 0;
  static int Cnt1 = 0;*/
  //----------------------------------------------------------------------------

  TA2WireChamber::Decode();

  for(Int_t i=0; i<fNchamber; i++)
    IntersectLayers(i);

  //----------------------------------------------------------------------------
  //This is Sven's Debug Code. Don't use it!
  /*fA1 = fA2 = fB1 = fB2 = EBufferEnd;

  TA2CylStripSven* sI0 = (TA2CylStripSven*)fWCLayer[fChamberLayers[0][1]]; //Inner
  TA2CylWireSven* sW0 = (TA2CylWireSven*)fWCLayer[fChamberLayers[0][2]];   //Inner
  TA2CylStripSven* sE0 = (TA2CylStripSven*)fWCLayer[fChamberLayers[0][3]]; //Inner
  TA2CylStripSven* sI1 = (TA2CylStripSven*)fWCLayer[fChamberLayers[1][1]]; //Outer
  TA2CylWireSven* sW1 = (TA2CylWireSven*)fWCLayer[fChamberLayers[1][2]];   //Outer
  TA2CylStripSven* sE1 = (TA2CylStripSven*)fWCLayer[fChamberLayers[1][3]]; //Outer

  //Always require "reasonable" phi range
  Double_t DPhiA = fabs(fPhiWire[0][0]-fPhiWire[0][1]);
  if((DPhiA>1.0)&&(DPhiA<5.0))
  //Case A: require 2 wire clusters in MWPC0 (inner) AND 2 points in MWPC0 (inner) ...
  if((sW0->GetNClust()==2)&&(fNIntersect[0]==2))
  {
    fA1 = 1;
    Cnt0++;
    //... AND 2 or more points in MWPC1 (outer)
    if(fNIntersect[1]>=2)
    {
      fA2 = 1;
      Cnt1++;
    }
    printf("%d : %d %d : %f %f : %d %d %d  : %d - %d %d %d : %d\n",
           Count, Cnt0, Cnt1, fPhiWire[0][0], fPhiWire[0][1],
           sW1->GetNClust(), sI1->GetNClust(), sE1->GetNClust(), fNIntersect[1],
           sW0->GetNClust(), sI0->GetNClust(), sE0->GetNClust(), fNIntersect[0]);
    if(Count==387){
    printf("i0: ");
    for(Int_t x=0; x<sI0->GetNHits(); x++)
      printf("%d ",(sI0->GetLayerHits())[x]);
    printf("\n");
    printf("w0: ");
    for(Int_t x=0; x<sW0->GetNHits(); x++)
      printf("%d ",(sW0->GetLayerHits())[x]);
    printf("\n");
    printf("e0: ");
    for(Int_t x=0; x<sE0->GetNHits(); x++)
      printf("%d ",(sE0->GetLayerHits())[x]);
    printf("\n");
    printf("i1: ");
    for(Int_t x=0; x<sI1->GetNHits(); x++)
      printf("%d ",(sI1->GetLayerHits())[x]);
    printf("\n");
    printf("w1: ");
    for(Int_t x=0; x<sW1->GetNHits(); x++)
      printf("%d ",(sW1->GetLayerHits())[x]);
    printf("\n");
    printf("e1: ");
    for(Int_t x=0; x<sE1->GetNHits(); x++)
      printf("%d ",(sE1->GetLayerHits())[x]);
    printf("\n");}
  }
  Count++;

  //Always require "reasonable" phi range
  Double_t DPhiB = fabs(fPhiWire[1][0]-fPhiWire[1][1]);
  if((DPhiB>1.0)&&(DPhiB<5.0))
  //Case B: require 2 wire clusters in MWPC1 (outer) AND 2 points in MWPC1 (outer)
  //        AND Z coordinates of MWPC1 (outer) within  MWPC 0 (inner) ...
  if((sW1->GetNClust()==2)&&(fNIntersect[1]==2))
  if((sE0->IsInside(fZIntersect[1][0]))&&(sE0->IsInside(fZIntersect[1][1])))
  {
    fB1 = 1;
    //... AND 2 or more points in MWPC0 (inner)
    if(fNIntersect[0]>=2)
      fB2 = 1;
  }*/
  //----------------------------------------------------------------------------

  MakeTracks();
}

#endif
