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

#include "TA2CylMWPC.h"
#include "TA2CrystalBall.h"

//*****************************************************************************************

enum {
  ENWCLayer = 200, ENWCChamber, ENWCLayersInChamber, EWCChamberParm,
  EWCTypePlaneWire, EWCTypeCylWire, EWCTypePlaneDrift,
  EWCTypePlaneStrip, EWCTypeCylStrip,
};

// Command-line key words which determine what to read in
static const Map_t kWCKeys[] = {
  {"Number-Layers:",     ENWCLayer},
  {"Number-Chambers:",   ENWCChamber},
  {"Layers-In-Chamber:", ENWCLayersInChamber},
  {"Chamber-Parameter:", EWCChamberParm},
  {"Plane-Wire:",        EWCTypePlaneWire},
  {"Cyl-Wire:",          EWCTypeCylWire},
  {"Plane-Drift:",       EWCTypePlaneDrift},
  {"Plane-Strip:",       EWCTypePlaneStrip},
  {"Cyl-Strip:",         EWCTypeCylStrip},
  {NULL,          -1}
};

ClassImp(TA2CylMWPC)

//*****************************************************************************************

TA2CylMWPC::TA2CylMWPC( const char* name, TA2System* apparatus )
  :TA2WireChamber(name, apparatus)
{
  // Initial setting of local variables
  // mostly set undefined
  fDPhiMax = 0.08;
  fR = fRtE = fRtI = fNextPhiInt = fC1 = fC2 = NULL;
  fZIntersect = fPhiSWDiff = fPhiWire = NULL;
  fVertX = fVertY = fVertZ = fVertR = fPhiDiff =
    fPseudoVertX = fPseudoVertY = fPseudoVertZ = fPseudoVertR = NULL;
  fNIntersect = NULL;

  for(Int_t t=0; t<4; t++)
    for(Int_t s=0; s<128; s++)
      DeadStrPatt[t][s] = 0;
  TargZShift = 0.0;
  Simulate = false;
}

//*****************************************************************************************

TA2CylMWPC::~TA2CylMWPC()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  DeleteArrays();
}

//*****************************************************************************************

void TA2CylMWPC::DeleteArrays()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  TA2WireChamber::DeleteArrays();
  if( fR ) delete fR;
  if( fRtE ) delete fRtE;
  if( fRtI ) delete fRtI;
  if( fC1 ) delete fC1;
  if( fC2 ) delete fC2;
  if( fD1 ) delete fD1;
  if( fD2 ) delete fD2;
  if( fNextPhiInt ) delete fNextPhiInt;
  for(Int_t i=0; i<fNchamber; i++){
    delete fZIntersect[i];
    delete fPhiSWDiff[i];
    delete fPhiWire[i];
  }
  if( fZIntersect ) delete fZIntersect;
  if( fPhiSWDiff ) delete fPhiSWDiff;
  if( fPhiWire ) delete fPhiWire;
  if( fPhiDiff ) delete fPhiDiff;
  for(Int_t i=0; i<fMaxTrack; i++){
    delete fTrack[i];
    delete fPseudoVertex[i];
  }
  if( fTrack ) delete fTrack;
  if( fPseudoVertex ) delete fPseudoVertex;
  if( fVertexLimits ) delete fVertexLimits;
  if( fVertX ) delete fVertX;
  if( fVertY ) delete fVertY;
  if( fVertZ ) delete fVertZ;
  if( fVertR ) delete fVertR;
  if( fPseudoVertX ) delete fPseudoVertX;
  if( fPseudoVertY ) delete fPseudoVertY;
  if( fPseudoVertZ ) delete fPseudoVertZ;
  if( fPseudoVertR ) delete fPseudoVertR;
  if( fDPhi12 ) delete fDPhi12;
}

//*****************************************************************************************

void TA2CylMWPC::SaveDecoded( )
{
  // Save decoded info to Root Tree file
}

//*****************************************************************************************

void TA2CylMWPC::PostInit( )
{
  // Do some further setup before any display specification

  Int_t* lch;
  TA2CylStripSven* sI;
  TA2CylStripSven* sE;
  TA2CylWireSven* wire;
  TVector3** pPlastic;

  fIsEnergy = fIsTime = fIsPos = ETrue;  // ensure flags set properly

  // Quantites required for each chamber
  fNIntersect = new Int_t[fNchamber];      // # intersects
  fR = new Double_t[fNchamber];            // chamber radii
  fRtI = new Double_t[fNchamber];          // for intersection calc
  fRtE = new Double_t[fNchamber];          // ditto
  fC1 = new Double_t[fNchamber];          // ditto
  fC2 = new Double_t[fNchamber];          // ditto
  fD1 = new Double_t[fNchamber];          // ditto
  fD2 = new Double_t[fNchamber];          // ditto
  fNextPhiInt = new Double_t[fNchamber];        // ditto
  fZIntersect = new Double_t*[fNchamber];  // strip intersections
  fPhiSWDiff = new Double_t*[fNchamber];
  fPhiWire = new Double_t*[fNchamber];     // wire phi

  for(Int_t ich=0; ich<fNchamber; ich++){  // for N chambers;
    lch = fChamberLayers[ich];             // 0th element is # layers
    sI = (TA2CylStripSven*)fWCLayer[ lch[1] ];
    wire = (TA2CylWireSven*)fWCLayer[ lch[2] ];
    sE = (TA2CylStripSven*)fWCLayer[ lch[3] ];
    fR[ich] = wire->GetRadius();
    fRtI[ich] = sI->GetRadius()/sI->GetTgWC();
    fRtE[ich] = sE->GetRadius()/sE->GetTgWC();
    fC1[ich] = fPi*sI->GetTgWC()/wire->GetRadius();
    fC2[ich] = fPi*sI->GetRadius()*sE->GetRadius()/
      (wire->GetRadius()*sI->GetTgWC());
    fD1[ich] = (sE->GetZCor(0)-sI->GetZCor(0))/(fRtI[ich]+ fRtE[ich]);
    fD2[ich] = (sE->GetZCor(0)*fRtI[ich] + sI->GetZCor(0)*fRtE[ich])/(fRtI[ich]+ fRtE[ich]);
    // fRtI[ich] = sI->GetTgWC()/sI->GetRadius();
    // fRtE[ich] = sE->GetTgWC()/sE->GetRadius();
    fNextPhiInt[ich] = f2Pi/(1 - fRtI[ich]/fRtE[ich]);
    fZIntersect[ich] = new Double_t[fMaxIntersect+1];
    fPhiSWDiff[ich] = new Double_t[fMaxIntersect+1];
    fPhiWire[ich] = new Double_t[fMaxIntersect+1];
  }

  // For holding track info
  fTrack = new TA2Track*[fMaxTrack+1];
  fPseudoVertex = new TVector3*[fMaxTrack+1];
  for(Int_t i=0; i<fMaxTrack; i++) fTrack[i] = new TA2Track("CylTracks");
  for(Int_t i=0; i<fMaxTrack; i++)
    fPseudoVertex[i] = new TVector3(0,0,0);
  fPseudoVertX = new Double_t[fMaxTrack+1];
  fPseudoVertY = new Double_t[fMaxTrack+1];
  fPseudoVertZ = new Double_t[fMaxTrack+1];
  fPseudoVertR = new Double_t[fMaxTrack+1];
  fPhiDiff     = new Double_t[fMaxTrack+1];
  fTrackTheta  = new Double_t[fMaxTrack+1];
  fTrackPhi    = new Double_t[fMaxTrack+1];
  fDPhi12      = new Double_t[fMaxTrack+1];

  // Container for track vertex information
  fVertex = new TVector3*[fMaxVertex+1];
  for( Int_t j=0; j<fMaxVertex; j++ ) fVertex[j] = new TVector3(0,0,0);
  fVertX = new Double_t[fMaxVertex+1];
  fVertY = new Double_t[fMaxVertex+1];
  fVertZ = new Double_t[fMaxVertex+1];
  fVertR = new Double_t[fMaxVertex+1];

  //  fCB = (TA2Calorimeter*)((TA2Analysis*)fParent);
  //  fPID = ((TA2Calorimeter*)fCB)->GetPID();
  fCB = (TA2Calorimeter*)fParent;
  fPID = ((TA2CrystalBall*)fCB)->GetPID();
  pPlastic = fPID->GetPosition();
  for(Int_t t=0; t<24; t++)
    fPhiPlastic[t]= pPlastic[t]->Z();

  pRandoms = new TRandom3;
  pRandoms->SetSeed(0); //'Randomize Timer'

  TA2Detector::PostInit();                // standard detector stuff
}

//*****************************************************************************************

void TA2CylMWPC::SetChamberParm( Int_t n, Double_t* parm )
{
  // Store some further chamber parameters
  // There should be 7 provided
  if( n < 8 ) {
    PrintError("","<Invalid DAPHNE-chamber parameters supplied>");
    return;
  }
  fMaxIntersect = (Int_t)parm[0];       // max # strip intersect points
  fMaxTrack = (Int_t)parm[1];           // max # particle tracks
  fMaxVertex = (Int_t)parm[2];          // max # vertices
  fDPhiMax = parm[3];                   // max diff phi wire and phi intersect
  // Vertex limits i=0,1,2  Rmax, +Zmax -Zmax;
  // Confine vertices to target area
  fVertexLimits = new Double_t[3];
  for( Int_t i=0; i<3; i++ ) fVertexLimits[i] = parm[i+4];
  fMaxVertRadius = fVertexLimits[0];
  fMaxVertZ      = fVertexLimits[1];
  fMinVertZ      = fVertexLimits[2];
  fDPhi12Max = parm[7];
}

//*****************************************************************************************

void TA2CylMWPC::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  // Specialised chamber diagnostic stuff
  Char_t name[256];
  Char_t index[8];
  Char_t* varname;

  // standard wire chamber stuff
  TA2WireChamber::LoadVariable();

  // cylindrical wire chamber stuff... intersection coordinates
  // loop round each chamber
  for( Int_t n=0; n<fNChamber; n++ ){
    sprintf(index,"%d",n);
    strcpy(name,"ZIntersect");
    strcat(name,index);
    varname = new Char_t[strlen(name)+1];
    strcpy(varname,name);
    TA2DataManager::LoadVariable(varname, fZIntersect[n], EDMultiX);
    strcpy(name,"PhiSWDiff");
    strcat(name,index);
    varname = new Char_t[strlen(name)+1];
    strcpy(varname,name);
    TA2DataManager::LoadVariable(varname, fPhiSWDiff[n], EDMultiX);
    strcpy(name,"PhiWire");
    strcat(name,index);
    varname = new Char_t[strlen(name)+1];
    strcpy(varname,name);
    TA2DataManager::LoadVariable(varname, fPhiWire[n], EDMultiX);
  }
  // Coordinates of vertices
  TA2DataManager::LoadVariable("Xvertex", fVertX, EDMultiX);
  TA2DataManager::LoadVariable("Yvertex", fVertY, EDMultiX);
  TA2DataManager::LoadVariable("Zvertex", fVertZ, EDMultiX);
  TA2DataManager::LoadVariable("Rvertex", fVertR, EDMultiX);
  // Coordinates of pseudo vertices (single track closest approach to Z axis)
  TA2DataManager::LoadVariable("PsXvertex", fPseudoVertX, EDMultiX);
  TA2DataManager::LoadVariable("PsYvertex", fPseudoVertY, EDMultiX);
  TA2DataManager::LoadVariable("PsZvertex", fPseudoVertZ, EDMultiX);
  TA2DataManager::LoadVariable("PsRvertex", fPseudoVertR, EDMultiX);
  TA2DataManager::LoadVariable("PhiDiff",   fPhiDiff,     EDMultiX);
  TA2DataManager::LoadVariable("DPhi12",    fDPhi12,     EDMultiX);
  //Spare testing stuff
  TA2DataManager::LoadVariable("STesting", &fSTesting, EDSingleX);
  //----------------------------------------------------------------------------
  //This is Sven's Debug Code. Don't use it!
  /*TA2DataManager::LoadVariable("A1", &fA1, EISingleX);
  TA2DataManager::LoadVariable("A2", &fA2, EISingleX);
  TA2DataManager::LoadVariable("B1", &fB1, EISingleX);
  TA2DataManager::LoadVariable("B2", &fB2, EISingleX);
  TA2DataManager::LoadVariable("N1", &fNIntersect[0], EISingleX);
  TA2DataManager::LoadVariable("N2", &fNIntersect[1], EISingleX);*/
  //----------------------------------------------------------------------------
}

//*****************************************************************************************

void TA2CylMWPC::IntersectLayers(Int_t nch)
{
  // For chamber nch, inner hit cluster clI, external hit cluster clE
  // Find intersection points defined by hits on 2 helically-wound strips
  // of a cylindrical cathode "plane". Chamber layers 1,3 are inner,
  // external strip planes; layer 2 is the wire plane.
  // Intersection points are separated by Pi (in phi) and ambiguities are
  // resolved using wire or PID info.

  TA2CylStripSven* sI = (TA2CylStripSven*)fWCLayer[fChamberLayers[nch][1]]; //Inner strip plane
  TA2CylWireSven* sW = (TA2CylWireSven*)fWCLayer[fChamberLayers[nch][2]]; //Wire plane between strip planes
  TA2CylStripSven* sE = (TA2CylStripSven*)fWCLayer[fChamberLayers[nch][3]]; //External strip plane

  Double_t phiW,phiI,phiE,z0,phi0;
  Double_t z[2];
  Double_t phi[2];
  Int_t j,k;
  Int_t nint = 0; //# intersects
  Int_t nsol;

  Int_t nIClust;
  Int_t nWClust;
  Int_t nEClust;

  Int_t MaxClust;
  Int_t nTemp;
  Double_t OpStr;
  Int_t* pPIDElem;

//------------------------------------------------------------------------------
  for(Int_t t=0;t<sI->GetNClust();t++)
    fClustEnI[nch][t]=sI->GetClustEn(t);
  for(Int_t t=0;t<sE->GetNClust();t++)
    fClustEnE[nch][t]=sE->GetClustEn(t);

  fClustEnI[nch][sI->GetNClust()]=EBufferEnd;
  fClustEnE[nch][sE->GetNClust()]=EBufferEnd;

  fZIntersect[nch][0] = EBufferEnd;
  fPhiSWDiff[nch][0] = EBufferEnd;
  fPhiWire[nch][0] = EBufferEnd;
  fNIntersect[nch] = 0;
//------------------------------------------------------------------------------

  MaxClust = ((TA2WCLayerSven*)fWCLayer[fChamberLayers[1][2]])->GetNElement();
  Int_t Wused[MaxClust]; //Assistant local variables:
  Int_t Iused[MaxClust]; //The same dimension as cluster's arrays dimensions; the meaning:
  Int_t Eused[MaxClust]; //how many times the cluster was used for the point production;

  for(Int_t i=0; i<sW->GetNClust(); i++) Wused[i]=0;
  for(Int_t i=0; i<sI->GetNClust(); i++) Iused[i]=0;
  for(Int_t i=0; i<sE->GetNClust(); i++) Eused[i]=0;

  nIClust = sI->GetNClust();
  nWClust = sW->GetNClust();
  nEClust = sE->GetNClust();

//------------------------------------------------------------------------------
  //Step 1: we have i/e strips and wire information:
  //Loop over the decoded hit clusters in inner and outer strip layers
  for(Int_t clI=0; clI<sI->GetNClust(); clI++)
    for(Int_t clE=0; clE<sE->GetNClust(); clE++)
    {
      phiI = sI->GetCGClust(clI);
      phiE = sE->GetCGClust(clE);
      nsol = 0;
      for(j=-2; j<=0; j++)
        for(k=-1; k<=1; k++)
        {
          z0 = fC2[nch]*(phiI+phiE+j+k) + fD2[nch];
          if(sE->IsInside(z0))
          {
            phi0 = fC1[nch]*(fRtE[nch]*(phiE+k) - fRtI[nch]*(phiI+j)) + fD1[nch];
            if(fabs(phi0)<0.0001) phi0=0.0;
            if(fabs(phi0-f2Pi)<0.0001) phi0=f2Pi;
            if((phi0>=0.0) && (phi0<=f2Pi))
            {
              z[nsol] = z0;
              phi[nsol] = phi0;
              nsol++;
              if(nsol==2) goto GotStripIntersects;
            }
          }
        }
GotStripIntersects:
      // two intersects possible....which if any valid
      // resolve using wire phi info
      for(k=0; k<nsol; k++)
        for(j=0; j<sW->GetNClust(); j++)
        {
          if(Wused[j]) continue;
          if(Iused[clI]) continue;
          if(Eused[clE]) continue;
          phiW = sW->GetCGClust(j);
          //printf("\n %d %d %d  : |%f - %f| = %f",
          //      sI->GetNClust(), sW->GetNClust(), sE->GetNClust(),
          //      phiW, phi[k], AbsDiffPhiSven(phiW,phi[k]));
          if(AbsDiffPhiSven(phiW,phi[k])>=fDPhiMax) continue;
          fZIntersect[nch][nint] = z[k];
          fPhiWire[nch][nint] = phiW;
          fPhiSWDiff[nch][nint] = phi[k] - phiW;
          //Marks the used clusters
          Wused[j]++;
          Iused[clI]++;
          Eused[clE]++;
          nWClust--;
          nIClust--;
          nEClust--;
          nint++;
          //printf(" -> Accepted!");
          if(nint>=fMaxIntersect) goto BuffersFull;
          break;
        }
    }
    //printf("\n");
//------------------------------------------------------------------------------
  nTemp = nint;

  //1st iteration: Step 2; for each pair of unused wire and strip clusters finds wire-
  //               strip solution, and accepts the point if the corresponding strip in
  //               the opposite layer is dead
  //2nd iteration: Step 3; the same as before, but without dead elements checking
  for(Int_t r=0; r<2; r++)
  {
      if((nWClust==nEClust)&&(nWClust==nIClust+1))
      for(Int_t w=0; w<sW->GetNClust(); w++)
        for(Int_t e=0; e<sE->GetNClust(); e++)
        {
          if(Wused[w]) continue;
          if(Eused[e]) continue;
          phiW = sW->GetCGClust(w);
          if(StripWire(nch*2+1, sE->GetCGClust(e), phiW, z, &OpStr)==0) continue;
          if(!r) //In which iteration (=step) are we?
            if(DeadStrPatt[nch*2][(Int_t)rint(fabs(OpStr))]==0) continue;
          //accept the point:
          fZIntersect[nch][nint] = z[0];
          fPhiWire[nch][nint] = phiW;
          fPhiSWDiff[nch][nint] = 0.19; //This is a tag how the point has been reconstructed
          //Marks the used clusters
          Wused[w]++;
          Eused[e]++;
          nWClust--;
          nEClust--;
          nint++;
          if(nint>=fMaxIntersect) goto BuffersFull;
        }
      if((nWClust==nIClust)&&(nWClust==nEClust+1))
      for(Int_t w=0; w<sW->GetNClust(); w++)
        for(Int_t i=0; i<sI->GetNClust(); i++)
        {
          if(Wused[w]) continue;
          if(Iused[i]) continue;
          phiW = sW->GetCGClust(w);
          if(StripWire(nch*2, sI->GetCGClust(i), phiW, z, &OpStr)==0) continue;
          if(!r) //In which iteration (=step) are we?
            if(DeadStrPatt[nch*2+1][(Int_t)rint(fabs(OpStr))]==0) continue;
          //accept the point:
          fZIntersect[nch][nint] = z[0];
          fPhiWire[nch][nint] = phiW;
          fPhiSWDiff[nch][nint] = 0.15; //This is a tag how the point has been reconstructed
          //Marks the used clusters
          Wused[w]++;
          Iused[i]++;
          nIClust--;
          nWClust--;
          nint++;
          if(nint>=fMaxIntersect) goto BuffersFull;
        }
    nTemp = nint;
  }
//------------------------------------------------------------------------------
  //Step 4: Missing wire information, instead use PID elements as big 'wires'
  if((nEClust==nIClust)&&(nEClust==nWClust+1))
  {
    // Loop over the decoded hit clusters in inner and outer strip layers
    for(Int_t clI=0; clI<sI->GetNClust(); clI++)
      for(Int_t clE=0; clE<sE->GetNClust(); clE++)
      {
        phiI = sI->GetCGClust(clI);
        phiE = sE->GetCGClust(clE);
        nsol = 0;
        for(j=-2; j<=0; j++)
        {
          for(k=-1; k<=1; k++)
          {
            z0 = fC2[nch]*(phiI+phiE+j+k) + fD2[nch];
            if(sI->IsInside(z0))
            {
              phi0 = fC1[nch]*(fRtE[nch]*(phiE+k) - fRtI[nch]*(phiI+j)) + fD1[nch];
              if((phi0>=0) && (phi0<f2Pi))
              {
                z[nsol] = z0;
                phi[nsol] = phi0;
                nsol++;
                if(nsol==2) goto GotStripIntersects2;
              }
            }
          }
        }
GotStripIntersects2:
        for(j=0; j<(Int_t)fPID->GetNhits(); j++)
        {
          //Use PID element as 'wire':
          if(Iused[clI]) continue;
          if(Eused[clE]) continue;
          pPIDElem = fPID->GetHits();
          phiW = fPhiPlastic[pPIDElem[j]] / RADDEG; //Get phi angle of hit PID element;
          phiW+=PI;                                 //convert to rad and correct coordinate
          if(phiW>TWOPI) phiW-=TWOPI;               //system (-pi..pi -> 0..2pi)
          for(k=0; k<nsol; k++)
          {
            if(AbsDiffPhiSven(phiW,phi[k])>=(15.0/RADDEG)) continue;
            fZIntersect[nch][nint] = z[k];
            fPhiWire[nch][nint] = phi[k];
            fPhiSWDiff[nch][nint] = 0.17; //This is a tag how the point has been reconstructed
            //Marks the used clusters
            Iused[clI]++;
            Eused[clE]++;
            nIClust--;
            nEClust--;
            nint++;
            if(nint>=fMaxIntersect) goto BuffersFull;
            break;
          }
        }
      }
  }
//------------------------------------------------------------------------------
BuffersFull:
  //Write end markers to intersection-point arrays
  fZIntersect[nch][nint] = EBufferEnd;
  fPhiSWDiff[nch][nint] = EBufferEnd;
  fPhiWire[nch][nint] = EBufferEnd;
  fNIntersect[nch] = nint;
//------------------------------------------------------------------------------
  //Correct different phi coordinate systems between CB and WC;
  //perform global shift along z axis.
  for(Int_t t=0; t<nint; t++)
  {
    fPhiWire[nch][t]-=PI;
    fZIntersect[nch][t]+=TargZShift;
  }
//------------------------------------------------------------------------------
  return;
}

//*****************************************************************************************

void TA2CylMWPC::MakeTracks()
{
  // Join r,phi,z coordinates from chambers to produce a "track"
  // 2 3-vectors, origin (inner chamber x,y,z) and direction cosines

  Int_t i,j;
  Double_t z,z0,phi,phi0,r;
  Int_t n = 0;
  TVector3* vertex;

  // Find & store all possible tracks
  // loop over intersections in inside chamber
  for(i=0; i<fNIntersect[0]; i++)
  {
    z0 = fZIntersect[0][i];                     // z, phi inner chamber
    phi0 = fPhiWire[0][i];
    // loop over outer chamber intersections
    for(j=0; j<fNIntersect[1]; j++)
    {
      z = fZIntersect[1][j];                    // z, phi outer chamber
      phi = fPhiWire[1][j];
      // Demand inner and outer phi within specified limit
      fDPhi12[n] = AbsDiffPhiSven(phi0,phi) * RADDEG;
      if(AbsDiffPhiSven(phi0,phi) > fDPhi12Max) continue;
      fTrack[n]->SetTrackCyl(fR[0],phi0,z0,fR[1],phi,z);
      // Look for closest approach to Z axis
      // track must intersect cylinder radius fMaxVertRadius about Z=0
      vertex = fTrack[n]->PseudoZVertex(fMaxVertRadius);
      if(vertex==NULL) continue;
      // check if pseudo-vertex within z limits
      z = vertex->Z();
      if( z < fMinVertZ ) continue;     // -z limit
      if( z > fMaxVertZ ) continue;     // +z limit
      fPhiDiff[n] = phi0 - phi;
      *fPseudoVertex[n] = *vertex;
      fPseudoVertX[n] = vertex->X();
      fPseudoVertY[n] = vertex->Y();
      fPseudoVertZ[n] = vertex->Z();
      fPseudoVertR[n] =
        TMath::Sqrt(vertex->X()*vertex->X() + vertex->Y()*vertex->Y());
      fTrackTheta[n] = fTrack[n]->GetTheta();
      fTrackPhi[n] = fTrack[n]->GetPhi();
      n++;                                     // track accepted
      if(n>=fMaxTrack) goto FinishedTracks;
    }
  }
FinishedTracks:
  // end markers to the pseudo vertex coordinate buffers
  fPseudoVertX[n] = fPseudoVertY[n] = fPseudoVertZ[n] = fPseudoVertR[n] =
    fTrackTheta[n] = fTrackPhi[n] = fPhiDiff[n] = fDPhi12[n] = EBufferEnd;
  fNTrack = n;
  // Now if there are >1 tracks, attempt to reconstruct vertices
  n = 0;
//   for(i=1; i<fNTrack; i++){
//     for(j=0; j<fNTrack; j++){
//       if(j==i) continue;                        // same track
//       vertex = fTrack[i]->TrackVertex( fTrack[j] );
//       z = vertex->Z();
//       r = TMath::Sqrt( vertex->X() * vertex->X() + vertex->Y() * vertex->Y() );
//       if( r > fMaxVertRadius ) continue;            // r limit
//       if( z < fMinVertZ ) continue;                 // -z limit
//       if( z > fMaxVertZ ) continue;                 // +z limit
//       *fVertex[n] = *vertex;
//       fVertX[n] = vertex->X();
//       fVertY[n] = vertex->Y();
//       fVertZ[n] = vertex->Z();
//       fVertR[n] = r;
//       n++;                                   // vertex accepted
//       if( n >= fMaxVertex ) goto FinishedVertices; // no more storage
//     }
//   }
 FinishedVertices:
  fVertX[n] = fVertY[n] = fVertZ[n] = fVertR[n] = EBufferEnd;
  fNVertex = n;
}

//*****************************************************************************************

void TA2CylMWPC::SetConfig( char* line, int key )
{
  // Load config parameters from file or command line
  // Keywords which specify a type of command can be found in
  // the kWCKeys array at the top of the source .cc file


  Int_t n, nelem, maxcl, maxclsize;
  Int_t layer[8];
  Int_t DefStrip;
  Int_t DefLayer;
  Double_t dparm[16];
  Char_t name[64];

  // Cluster specific configuration
  switch( key ){
  case ENWCLayer:
    // no. of active layers (or planes) in chamber
    if( sscanf( line, "%d", &fNLayer ) < 1 ) goto error;
    fWCLayer = new TA2WCLayer*[fNLayer];
    break;
  case ENWCChamber:
    // no. of chambers in detector
    if( sscanf( line, "%d", &fNChamber ) < 1 ) goto error;
    fChamberLayers = new Int_t*[fNChamber];
    break;
  case ENWCLayersInChamber:
    // layers in particular chamber
    if( fNchamber >= fNChamber ){
      PrintError( line, "<Too many WC layers-in-chamber input>");
      return;
    }
    if( ( n = sscanf( line, "%d%d%d%d%d%d%d%d",
                      layer,layer+1,layer+2,layer+3,
                      layer+4,layer+5,layer+6,layer+7 ) ) < 1 ) goto error;
    fChamberLayers[fNchamber]  = new Int_t[n+1];
    fChamberLayers[fNchamber][0]  = n;
    for( Int_t i=1; i<=n; i++ ) fChamberLayers[fNchamber][i] = layer[i-1];
    fNchamber++;
    break;
  case EWCChamberParm:
    // wild-card chamber setup....depends on specific SetChamberParm
    n = sscanf( line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                dparm,    dparm+1,  dparm+2,  dparm+3,
                dparm+4,  dparm+5,  dparm+6,  dparm+7,
                dparm+8,  dparm+9,  dparm+10, dparm+11,
                dparm+12, dparm+13, dparm+14, dparm+15 );
    SetChamberParm( n, dparm );
    break;
  case EWCTypePlaneWire:
    //I'm abusing ths entry for parsing different stuff for MWPC simulation. Introducing a new
    //command-line keyword does not work due to the stupid idea of somebody to define the map
    //and keys also in class TA2WireChamber, so that introducing new words leads to segfaults!
    TA2CylMWPC::ParseMisc(line);
    break;
  case EWCTypeCylWire:
    // Wires arranged in cylinder
    // Has hit TDC, but no pulse height
    // layer parameters.....
    // radius, length,
    // quadratic z correction  coeff 0,1,2
    if( fNlayer >= fNLayer ){
      PrintError( line, "<Too many WC layers input>");
      return;
    }
    fIsEnergy = EFalse;
    fIsTime = fIsPos = ETrue;
    if( sscanf( line, "%s%d%d%d%lf%lf%lf%lf%lf",
                name, &nelem, &maxcl, &maxclsize, dparm, dparm+1,
                dparm+2, dparm+3, dparm+4 ) < 9 ) goto error;
    fWCLayer[fNlayer] = new TA2CylWireSven( name, nelem, maxcl, maxclsize,
                                        this, dparm );
    fNlayer++;
    break;
  case EWCTypePlaneStrip:
    //I'm abusing ths entry for defective cylindrical strips. Introducing a new command-line keyword
    //does not work due to the stupid idea of somebody to define the map and keys also in class
    //TA2WireChamber, so that introducing new words leads to segfaults!
    if(sscanf(line, "%d%d", &DefLayer, &DefStrip) < 2) goto error;
    DeadStrPatt[DefLayer][DefStrip] = 1;
    break;
  case EWCTypePlaneDrift:
    //I'm abusing ths entry for a general shift along z-coordinate. Introducing a new command-line
    //keyword does not work due to the stupid idea of somebody to define the map and keys also in
    //class TA2WireChamber, so that introducing new words leads to segfaults!
    if(sscanf(line, "%lf", &TargZShift) < 1) goto error;
    break;
  case EWCTypeCylStrip:
    // Cathode strips helically "wound" on cylinder
    // Has pulse height but no time
    // layer parameters.....
    // radius, length, TgWC, Z0, pitch,
    // quadratic z correction  coeff 0,1,2
    // quadratic phi correction coeff 0,1,2
    if( fNlayer >= fNLayer ){
      PrintError( line, "<Too many WC layers input>");
      return;
    }
    fIsEnergy = fIsPos = ETrue;
    fIsTime  = EFalse;
    if( sscanf( line, "%s%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                name, &nelem, &maxcl, &maxclsize,
                dparm, dparm+1, dparm+2, dparm+3, dparm+4,
                dparm+5, dparm+6, dparm+7,
                dparm+8, dparm+9, dparm+10 ) < 15 ) goto error;
    fWCLayer[fNlayer] = new TA2CylStripSven( name, nelem, maxcl, maxclsize, //DESWEGEN!!!
                                             this, dparm );
    fNlayer++;
    break;
  default:
    //fIsConfigPass = ETrue;
    // Default detector configuration
    // If this flags it has done nothing carry on with local setup
    TA2Detector::SetConfig( line, key );
    //if( !fIsConfigPass ) return;
    break;;
  }
  return;
 error: PrintError( line, "<WC configuration>" );
  return;
}

//*****************************************************************************************

Int_t TA2CylMWPC::StripWire(Int_t lay, Double_t str, Double_t phiW, Double_t ZetSol[], Double_t* Strip)
{
  Int_t nsol;
  Int_t ik,inner,lay1;
  Double_t f2,z,zet,strip;
  Double_t fZ0Abs;
  TA2CylStripSven* pStrips;

  //Get pointers to different WC layers
  TA2CylStripSven* sI = (TA2CylStripSven*)fWCLayer[ fChamberLayers[lay/2][1] ]; //Inner strip plane
  TA2CylStripSven* sE = (TA2CylStripSven*)fWCLayer[ fChamberLayers[lay/2][3] ]; //External strip plane

  //Get number of hit clusters in different WC layers
  Int_t nsi = sI->GetNElement();
  Int_t nse = sE->GetNElement();

  //Select layer for getting parameters
  if(lay%2)
    pStrips = sE;
  else
    pStrips = sI;

  fZ0Abs = 0.0; //I don't know from where to take it

  nsol = 0;
  ik = lay/2; //chamber's number 0,1;
  if(lay%2==0) //inner layer
  {
    inner=1;
    lay1=lay+1;
    f2=TWOPI*fRtI[ik];
    zet= fRtI[ik]*phiW + f2*str + pStrips->GetZ0();
  }
  else
  {
    inner=0;
    lay1=lay-1;
    f2=TWOPI*fRtE[ik];
    zet=-fRtE[ik]*phiW + f2*str + pStrips->GetZ0();
  }

  if(zet>pStrips->GetEffLength())
  {
    zet = zet - f2;
    if(zet>pStrips->GetEffLength())
      zet = zet - f2;
  }
  if(zet<-pStrips->GetEffLength())
    zet = zet + f2;

  z = zet + zet*zet*pStrips->GetZCor(2) + zet*pStrips->GetZCor(1) + pStrips->GetZCor(0);
  if((z<pStrips->GetEffLength()) && (z>-pStrips->GetEffLength()))
    nsol=1;

  ZetSol[0] = z + fZ0Abs;

  if(inner)
  {
    strip = (zet+fRtE[ik]*phiW)/(TWOPI*fRtE[ik])*nse;
    if(strip<0)
      strip=strip+nse;
    if(strip>=nse)
      strip=strip-nse;
  }
  else
  {
    strip = (zet-fRtI[ik]*phiW)/(TWOPI*fRtI[ik])*nsi;
    if(strip<0)
      strip = strip+nsi;
    if(strip<0)
      strip = strip+nsi;  //sic!
    if(strip>=nsi)
      strip = strip-nsi;
  }
  *Strip=strip;

  return nsol;
}

//*****************************************************************************************

void TA2CylMWPC::ParseMisc(char* line)
{
  char sWord[256];

  //Get keyword
  if(sscanf(line, "%s", sWord)!=1) return;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if(!strcmp("SimulateProton", sWord))
  {
    Simulate = true;
    printf("MWPC proton simulation enabled\n");
    return;
  }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if(!strcmp("EfficiencyFile", sWord) && Simulate)
  {
    sscanf(line, "%*s %s", FileEfficiency);
    printf("MWPC efficiencies in simulation from:\n %s\n", FileEfficiency);
    return;
  }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if(!strcmp("SigmaPhiFile", sWord) && Simulate)
  {
    sscanf(line, "%*s %s", FileSigmaPhi);
    printf("MWPC phi resolutions in simulation from:\n %s\n", FileSigmaPhi);
    return;
  }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if(!strcmp("SigmaThetaFile", sWord) && Simulate)
  {
    sscanf(line, "%*s %s", FileSigmaTheta);
    printf("MWPC phi resolutions in simulation from:\n %s\n", FileSigmaTheta);
    return;
  }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if(!strcmp("ReadParamFiles()", sWord) && Simulate)
  {
    ReadParamFiles();
    return;
  }
}

//*****************************************************************************************

void TA2CylMWPC::ReadParamFiles()
{
  Double_t Param;
  Int_t nTheta;

  FILE* EffFile = fopen(FileEfficiency, "r");
  for(Int_t t=0; t<36; t++)
  {
    fscanf(EffFile, "%d %lf", &nTheta, &Param);
    Efficiency[nTheta] = Param;
    if(Efficiency[nTheta]!=0.0)
      printf("Reading efficiency for MWPC in simulation (Theta = %d...%d)\n", nTheta*5, nTheta*5+5);
  }
  fclose(EffFile);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FILE* ThetaFile = fopen(FileSigmaTheta, "r");
  for(Int_t t=0; t<36; t++)
  {
    fscanf(ThetaFile, "%d %lf", &nTheta, &Param);
    SigmaTheta[nTheta] = Param * DEGRAD;
    if(SigmaTheta[nTheta]!=0.0)
      printf("Reading Theta resolution for MWPC in simulation (Theta = %d...%d)\n", nTheta*5, nTheta*5+5);
  }
  fclose(ThetaFile);

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FILE* PhiFile = fopen(FileSigmaPhi, "r");
  for(Int_t t=0; t<36; t++)
  {
    fscanf(PhiFile, "%d %lf", &nTheta, &Param);
    SigmaPhi[nTheta] = Param * DEGRAD;
    if(SigmaPhi[nTheta]!=0.0)
      printf("Reading Phi resolution for MWPC in simulation (Theta = %d...%d)\n", nTheta*5, nTheta*5+5);
  }
  fclose(PhiFile);
}

//*****************************************************************************************
