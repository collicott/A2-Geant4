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

#include "TA2TOF.h"
#include "TA2KensTagger.h"
#include "TA2KensLadder.h"
#include "TA2Analysis.h"

// Default list of detector classes that the TA2TAPS 
// apparatus may contain
enum { ETOFWall, ETOFVeto};
enum{EToFNLayers=100,EToFCoincTheta};

static Map_t kValidDetectors[] = {
  {"TA2LongScint",     ETOFWall},
  {"TA2LongScint",     ETOFVeto},
  {NULL, 		-1}
};

static const Map_t kTOFKeys[] = {
  {"NToFLayers:",            EToFNLayers},
  {"CoincTheta:",            EToFCoincTheta},
  {NULL,            -1}
};

ClassImp(TA2TOF)

//-----------------------------------------------------------------------------
TA2TOF::TA2TOF( const char* name, TA2System* analysis  )
  :TA2Apparatus( name, analysis, kValidDetectors )
{
  // Zero pointers and counters, add local SetConfig command list
  // to list of other apparatus commands
  fToFWall = NULL;
  fToFVeto = NULL;
  fCoincVetoE= NULL;
  fCoincToFE= NULL;
  fCoincToFT= NULL;
  fToFE= NULL;
  fToFT= NULL;
  fVetoTheta= NULL;
  fLayerTheta= NULL;
  fNcharged=0;
  fNLayers=0;
  fVetoThetaMax=0;
  fLayerThetaMax=0;
  AddCmdList( kTOFKeys );                  // for SetConfig()
}


//-----------------------------------------------------------------------------
TA2TOF::~TA2TOF()
{
  // Free up allocated memory
 
}

//-----------------------------------------------------------------------------
TA2DataManager*  TA2TOF::CreateChild(const char* name, Int_t dclass)
{
  // Create a TA2Detector class for use with the TA2TOF
  // Valid detector types are...
  // 1. TA2TOF_BaF2
  
  if( !name ) name = Map2String( dclass );
  switch (dclass){
  case ETOFWall:
    if(TString(name).Contains("Veto")){
      fToFVeto = new TA2LongScint( name, this );
      return fToFVeto;
    }
    else{
      fToFWall = new TA2LongScint( name, this );
      return fToFWall;
    }
  default:
    return NULL;
  }
}

//---------------------------------------------------------------------------
void TA2TOF::SetConfig( char* line, int key )
{
  // Load TOF parameters from file or command line
  // TOF specific configuration
  switch( key ){
  case EToFNLayers:
    if( sscanf( line, "%d", &fNLayers ) < 1 ) goto error;
   break;
  case EToFCoincTheta:
    if( sscanf( line, "%lf %lf", &fVetoThetaMax,&fLayerThetaMax ) < 1 ) goto error;
   break;
  default:
    // Command not found...possible pass to next config
    TA2Apparatus::SetConfig( line, key );
    break;;
  }
  return;
 error: PrintError( line );
  return;
}

//-----------------------------------------------------------------------------
void TA2TOF::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup.
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  //                            name        pointer          type-spec
 //
  TA2Apparatus::LoadVariable();
  TA2DataManager::LoadVariable("Ncharged",  &fNcharged,      EISingleX);
  TA2DataManager::LoadVariable("CoincVetoE",  fCoincVetoE,      EDMultiX);
  TA2DataManager::LoadVariable("CoincToFE",  fCoincToFE,      EDMultiX);
  TA2DataManager::LoadVariable("CoincToFT",  fCoincToFT,      EDMultiX);
  TA2DataManager::LoadVariable("ToFE",  fToFE,      EDMultiX);
  TA2DataManager::LoadVariable("ToFT",  fToFT,      EDMultiX);
  TA2DataManager::LoadVariable("VetoTheta",  fVetoTheta,      EDMultiX);
  TA2DataManager::LoadVariable("LayerTheta",  fLayerTheta,      EDMultiX);
}


//-----------------------------------------------------------------------------
void TA2TOF::PostInit( )
{
  // Initialise arrays used to correlate hits in BaF2 and Veto detectors.
  // Load 2D cuts file and get the contained cuts classes
  // Demand particle ID class...if not there self destruct
  // Does not come back if called with EErrFatal
   if( !fParticleID )
    PrintError("",
	       "<Configuration aborted: ParticleID class MUST be specified>",
	       EErrFatal);
  // CB-PID Phi diff array
  fMaxParticle = fToFWall->GetNBar();

  fCoincVetoE=new Double_t[fToFWall->GetNBar()];
  fCoincToFE=new Double_t[fToFWall->GetNBar()];
  fToFT=new Double_t[fToFWall->GetNBar()];
  fToFE=new Double_t[fToFWall->GetNBar()];
  fCoincToFT=new Double_t[fToFWall->GetNBar()];
  fVetoTheta=new Double_t[fToFWall->GetNBar()];
  fLayerTheta=new Double_t[fToFWall->GetNBar()*fToFWall->GetNBar()];
  for(Int_t i=0;i<fToFWall->GetNBar();i++){
    fCoincVetoE[i]=EBufferEnd;
    fCoincToFE[i]=EBufferEnd;
    fCoincToFT[i]=EBufferEnd;
    fToFE[i]=EBufferEnd;
    fToFT[i]=EBufferEnd;
    fVetoTheta[i]=EBufferEnd;
    for(Int_t j=0;j<fToFWall->GetNBar();j++)fLayerTheta[i*fToFWall->GetNBar()+j]=EBufferEnd;
  }

  TA2Apparatus::PostInit();
}

void TA2TOF::Reconstruct( ){

  Int_t *bhits=fToFWall->GetBarHits();
  Int_t *vhits=fToFVeto->GetBarHits();
  
  TVector3 pos;
  TVector3 pos2;
  TVector3 vpos;
  //Double_t ToFE;
  Int_t nlaycoinc=0;
  Int_t donebar[fToFWall->GetNBar()*fToFWall->GetNBar()];
  Int_t done=0;
  fNcharged=0;
  fNparticle=0;
  for(Int_t i=0;i<fToFWall->GetNBarHits();i++){
    fToFE[i]=0;
    fToFT[i]=0;
    for(Int_t db=0;db<nlaycoinc;db++) if(bhits[i]==donebar[db])done=1;
    if(done==1) continue;
    pos=fToFWall->GetBarPos(bhits[i]);
    fToFE[i]=fToFWall->GetMeanEnergy(bhits[i]);
    fToFT[i]=fToFWall->GetMeanTime(bhits[i]);
    Int_t in=((bhits[i])%(fToFWall->GetNBar()/fNLayers));
    //loop over other bars to see if these are position correlated
    for(Int_t i2=i+1;i2<fToFWall->GetNBarHits();i2++){
      pos2=fToFWall->GetBarPos(bhits[i2]);
      if(nlaycoinc<fToFWall->GetNBar()*fToFWall->GetNBar())fLayerTheta[nlaycoinc]=pos.Angle(pos2)*TMath::RadToDeg();
      if(fLayerTheta[nlaycoinc]<fLayerThetaMax){//prob from the same particle
	fToFE[i]+=fToFWall->GetMeanEnergy(bhits[i2]);
	donebar[nlaycoinc++]=bhits[i2]; //delete hit so don't double count
      }
    }
    //loop over veto to see if coincidence
    for(Int_t j=0;j<fToFVeto->GetNBarHits();j++){
      if(fNcharged>fToFVeto->GetNBarHits()) {fNcharged=0; break;}//too many hits!!
      vpos=fToFVeto->GetBarPos(vhits[j]);
      fVetoTheta[fNcharged]=pos.Angle(vpos)*TMath::RadToDeg();
      if(fVetoTheta[fNcharged]<fVetoThetaMax){//coinc with veto
	fCoincVetoE[fNcharged]=fToFVeto->GetMeanEnergy(vhits[j]);
	fCoincToFE[fNcharged]=fToFE[i];
	fCoincToFT[fNcharged]=fToFWall->GetMeanTime(bhits[i]);
	fNcharged++;
	fPDG_ID[i] = kRootino;//set to rootino for just now
      }

    }
    //Now try some PID
    //Cuts set in .dat file follows TA2CrystalBall algorithm
    if( fPCut ){
      Int_t m = fPCutStart[bhits[i]]; //Start cut index for bar bhits[i]
      for( Int_t n=0; n<fNSectorCut[bhits[i]]; n++ ){ // Loop over specified cuts
	if( fPCut[m]->Test() ){           // Condition met?
	  fPDG_ID[i] = GetCutPDGIndex(m); // OK get associated PDG index
	  break;                          // cut OK so exit loop
	}
	m++;                              // for next try
      }
    }
    // Set 4 vector on the basis of ID above
    //Calculate enrgy from ToF
    Double_t partT=0;
    if(fPDG_ID[i]==kProton||fPDG_ID[i]==kNeutron){
      Double_t beta=pos.Mag()/100./(fToFT[i]*1E-9)/2.9979E8;
      Double_t gamma=1/sqrt(1-beta*beta);
      Double_t partT=(gamma-1)*fParticleID->GetMassMeV(fPDG_ID[i]);
    }
    else partT= fToFE[i];
    fParticleID->SetP4( fP4+i,fPDG_ID[i],partT,&pos );
    fP4tot += fP4[i];                    // accumulate total 4 mom.
    fNparticle++;
  }
  fCoincVetoE[fNcharged]=EBufferEnd;
  fCoincToFE[fNcharged]=EBufferEnd;
  fCoincToFT[fNcharged]=EBufferEnd;
  fToFE[fToFWall->GetNBarHits()]=EBufferEnd;
  fToFT[fToFWall->GetNBarHits()]=EBufferEnd;
  fVetoTheta[fNcharged]=EBufferEnd;
  fLayerTheta[nlaycoinc]=EBufferEnd;

}
