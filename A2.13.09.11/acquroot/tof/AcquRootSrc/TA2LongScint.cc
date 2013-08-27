//--Author	JRM Annand   27th Apr 2003
//--Rev 	JRM Annand...30th Sep 2003  New Detector bases
//--Rev 	JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004  3v8 compatible
//--Rev 	JRM Annand...25th Jul 2005  SetConfig hierarchy
//--Update	JRM Annand...18th Jan 2007  Inherit from TA2Detector
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2LongScint
//
// Decoding and calibration methods for a long bar of plastic scintillator
// Assumes one PMT at each end of bar -> position sensitivity, mean time etc.
//

#include "TA2LongScint.h"

// constants for command-line maps below
enum { 
  ELongScintMakeBar = 200, ELongScintBarSize, ELongScintNVetos
};

// Command-line key words which determine what to read in
static const Map_t kLongScintKeys[] = {
  {"Bar:",                  ELongScintMakeBar},
  {"Bar-Size:",             ELongScintBarSize},
  {"NVetos:",             ELongScintNVetos},
  //  {"All-Neighbour:",        ELongScintAllNeighbour},
  {NULL,          -1}
};

ClassImp(TA2LongScint)

//---------------------------------------------------------------------------
TA2LongScint::TA2LongScint( const char* name, TA2System* apparatus )
  :TA2Detector(name, apparatus)
{
  // Call Cluster-detector constructer and
  // Initialise local variables to undefined state
  AddCmdList( kLongScintKeys );   // add long-scint specific config
  fBar = NULL;
  fMeanEnergy = NULL;
  fMeanEnergyOR = NULL;
  fMeanTime = NULL;
  fMeanTimeOR = NULL;
  fTimeDiff = NULL;
  fTimeDiffOR = NULL;
  fBarPos = NULL;
  fBarHits = NULL;
  fBarHitCount = NULL;
  fElement2Bar = NULL;
  fNBar = fNbar = fNBarHits = 0;
  fNVetos=0;
}

//---------------------------------------------------------------------------
TA2LongScint::~TA2LongScint()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  // Start with arrays of TA2Detector class
  DeleteArrays();

  if( fBar ){
    for(Int_t i=0; i<fNBar; i++){ delete fBar[i]; }
    delete fBar;
  }
  if( fMeanEnergy ) delete fMeanEnergy;
  if( fMeanEnergyOR ) delete fMeanEnergyOR;
  if( fMeanTime ) delete fMeanTime;
  if( fMeanTimeOR ) delete fMeanTimeOR;
  if( fTimeDiff ) delete fTimeDiff;
  if( fTimeDiffOR ) delete fTimeDiffOR;
  if( fBarPos ) delete fBarPos;
  if( fBarHits ) delete fBarHits;
  if( fBarHitCount ) delete fBarHitCount;
  if( fElement2Bar ) delete fElement2Bar;
}

//---------------------------------------------------------------------------
void TA2LongScint::SetConfig( Char_t* line, Int_t key )
{
  // Load Ladder parameters from file or command line
  // Keywords which specify a type of command can be found in
  // the kLaddKeys array at the top of the source .cc file
  // The following are setup...
  //	1. # X-triggers, #elements in ladder array, global TDC offset
  //	2. X-trigger parameters if any, TDC, TDC calib, TDC cut
  //	3. Ladder parameters, TDC, calib, cut window, Eg calib, Scaler
  //	4. Any post initialisation
  //	5. Histograms...should be done after post-initialisation

  Int_t n;
  // Cluster specific configuration
  switch( key ){
  case ELongScintBarSize:
    // Setup scintillator bar arrays
    // Mean time/energy, time diff, hits, position etc.
    if( !(n = sscanf( line, "%d", &fNBar )) ) fNBar = fNelement/2;
    fBar = new LongBar_t*[fNBar];
    fMeanEnergy = new Double_t[fNBar];
    fMeanEnergyOR = new Double_t[fNBar+1];
    fMeanTime = new Double_t[fNBar];
    fMeanTimeOR = new Double_t[fNBar+1];
    fTimeDiff = new Double_t[fNBar];
    fTimeDiffOR = new Double_t[fNBar+1];
    fBarPos = new TVector3[fNBar];       // no (0,0,0) gcc 3.4
    fBarHits = new Int_t[fNBar];
    fBarHitCount = new Int_t[fNBar];
    fElement2Bar = new Int_t[fNelement];
    break;
  case ELongScintMakeBar:
    // Construct scintillator bar from 2 elements
    if( fNbar >= fNBar ){
      PrintError(line,"<Too many scintillators bars>");
      break;
    }
    fBar[fNbar] = new LongBar_t( line, fNbar, this );
    fNbar++;
    break;
  case ELongScintNVetos:
    if( !(n = sscanf( line, "%d", &fNVetos )) ) PrintError(line,"<How many TOF vetos?>");
    break;
  default:
    // Try standard cluster commands
    TA2Detector::SetConfig( line, key );
    break;;
  }
  return;
}

//-----------------------------------------------------------------------------
void TA2LongScint::ParseDisplay( char* line )
{
  // Input private histogram spec to general purpose parse
  // and setup routine

  //  const Name2Variable_t hist[] = {
  // Name          ->variable  single/mult    Fill-condition
  //    {"Nphoton_Minv", fM_Nphoton, EHistSingleX},
  //    {NULL,          0,         0}
  //  };
  // Do not remove the final NULL line

  TA2Detector::ParseDisplay(line);
  return;
}

//---------------------------------------------------------------------------
void TA2LongScint::PostInit()
{
  // Some further initialisation after all setup parameters read in
  // Start with alignment offsets
  // Create space for various output arrays
  TA2Detector::PostInit();

  // Construct element-to-bar cross reference
  for(Int_t i=0; i<(Int_t)fNelement; i++){
    for(Int_t j=0; j<fNBar; j++){
      if( (fBar[j]->GetEnd1() == i) || (fBar[j]->GetEnd2() == i) ){
	fElement2Bar[i] = j; 
	break; 
      }
    }
  }
}

//-----------------------------------------------------------------------------
void TA2LongScint::LoadVariable( )
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

  TA2DataManager::LoadVariable("MeanEnergy",   fMeanEnergy,      EDSingleX);
  TA2DataManager::LoadVariable("MeanEnergyOR", fMeanEnergyOR,    EDMultiX);
  TA2DataManager::LoadVariable("MeanTime",     fMeanTime,        EDSingleX);
  TA2DataManager::LoadVariable("MeanTimeOR",   fMeanTimeOR,      EDMultiX);
  TA2DataManager::LoadVariable("TimeDiff",     fTimeDiff,        EDSingleX);
  TA2DataManager::LoadVariable("TimeDiffOR",   fTimeDiffOR,      EDMultiX);
  TA2DataManager::LoadVariable("BarHits",      fBarHits,         EIMultiX);
  TA2DataManager::LoadVariable("NBarHits",     &fNBarHits,       EISingleX);
  TA2Detector::LoadVariable();
}
//---------------------------------------------------------------------------
void TA2LongScint::SaveDecoded( )
{
  // Save decoded info to Root Tree file
}

