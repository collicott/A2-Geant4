//--Author	JRM Annand   30th Aug 2003
//--Rev 	A.Starostin..24th Jan 2004  add theta,phi
//--Rev 	JRM Annand   30th Jan 2004  mod neighbour determination
//--Rev 	JRM Annand   21st Oct 2004  cluster time
//--Rev 	JRM Annand    9th Mar 2005  protected instead of private vars
//--Rev 	JRM Annand   14th Apr 2005  energy fraction central crystal
//--Rev 	JRM Annand    6th Jun 2005  up to 48 nearest neighbours
//--Rev 	JRM Annand   13th Jul 2005  add Merge function
//--Update	JRM Annand   19th Oct 2005  up to 64 nearest neighbours
//--Update	D Glazier    24th Aug 2007  Add IsTime check
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// HitCluster_t
// Specimen for hit cluster determination in a segmented calorimeter
// Mean cluster position obtained from the sqrt(E) weighted sum
//
//---------------------------------------------------------------------------

#include "HitCluster_t.h"
#include "TA2ClusterDetector.h"

//---------------------------------------------------------------------------
HitCluster_t::HitCluster_t( Char_t* line, UInt_t index, Int_t sizefactor )
{
  // store input parameters
  // # inner nearest neighbours (outer calculated from total read)
  // coordinates of center of front face of central element
  // List of nearest neighbours inner & outer
  UInt_t hit[64];
  fIndex = index;
  UInt_t n = 
    sscanf( line, 
	    "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
	    &fNNearNeighbour,
	    hit,   hit+1, hit+2, hit+3, hit+4, hit+5, hit+6, hit+7,
	    hit+8, hit+9, hit+10,hit+11,hit+12,hit+13,hit+14,hit+15,
	    hit+16,hit+17,hit+18,hit+19,hit+20,hit+21,hit+22,hit+23,
 	    hit+24,hit+25,hit+26,hit+27,hit+28,hit+29,hit+30,hit+31,
	    hit+32,hit+33,hit+34,hit+35,hit+36,hit+37,hit+38,hit+39,
	    hit+40,hit+41,hit+42,hit+43,hit+44,hit+45,hit+46,hit+47, 
	    hit+48,hit+49,hit+50,hit+51,hit+52,hit+53,hit+54,hit+55,
	    hit+56,hit+57,hit+58,hit+59,hit+60,hit+61,hit+62,hit+63 );

  // Consistency check...1st hit must be the index
  if( (n < (fNNearNeighbour + 1)) || (index != *hit) ){
    printf(" Error in nearest neighbour input at line:\n %s\n", line );
    return;
  }
  n -= 2;                         // # neighbours around central element
  fNNeighbour = n;
  fNeighbour = new UInt_t[n];
  fMaxHits = n * sizefactor;
  fHits = new UInt_t[ fMaxHits ];
  fHits[0] = ENullHit;
  fNhits = 0;
  fEnergy = (Double_t)ENullHit;
  for( UInt_t i=0; i<n; i++ ) fNeighbour[i] = hit[i+1];
  fMeanPosition = new TVector3(0.0, 0.0, 0.0);
}

//---------------------------------------------------------------------------
HitCluster_t::~HitCluster_t( )
{
  // delete arrays
  if( fHits ) delete fHits;
  if( fNeighbour ) delete fNeighbour;
  if( fMeanPosition ) delete fMeanPosition;
}

//---------------------------------------------------------------------------
void HitCluster_t::ClusterDetermine( TA2ClusterDetector* cl )
{
  // Determine which hits form the cluster.
  // Basic algorithm: supply list of nearest neighbours for
  // each element in calorimeter. Any (not previously claimed) hit in one
  // of these neighbours is added to the cluster. Clusters are initially
  // centered on elements where the energy signal is a local maximum.
  // The fIsIterate = kTRUE flag signals that a further search for cluster-
  // connected hits which are NOT in the nearest neighbours list
  // should be made.
  // The position of the cluster is taken as the
  // sqrt(energy)-weighted mean of the individual positions of the
  // hit elements.
 
  UInt_t i,j,k;

  Double_t* energy = cl->GetEnergy();
  Double_t* time = cl->GetTime();
  UInt_t* hits = cl->GetTempHits();
  TVector3** pos = cl->GetPosition();
  UInt_t nhits = cl->GetNhits();

  *fMeanPosition = *(pos[fIndex]);            // position = "centre"
  fEnergy = energy[fIndex];                   // energy in "central" element
  if(cl->IsTime())fTime = time[fIndex];                       // time in central element
  Double_t sqrtE = sqrt(fEnergy);
  fSqrtEtot = sqrtE;
  *fMeanPosition = *(pos[fIndex]) * sqrtE;   // position = "centre"
  fHits[0] = fIndex;
  
  // Accumulate weighted mean position
  for( i=0,k=1; i<nhits; i++ ){
    if( (j = hits[i]) == ENullHit ) continue; // was previously counted
    if( j == fIndex ){
      hits[i] = ENullHit;
      continue;
    }                                         // already got center
    if( IsNeighbour(j) ){                     // a neighbour of the center?
      hits[i] = ENullHit;                     // so its not double counted
      fHits[k] = j;                           // add to cluster hits collection
      sqrtE = sqrt(energy[j]);
      fEnergy += energy[j];
      fSqrtEtot += sqrtE;
      *fMeanPosition += ( *(pos[j]) * sqrtE );// root energy weighted pos
      k++;
    }
  }
  fNhits = k;
  // If the iterate flag is TRUE check if any unused hits (ie NOT in the near
  // neighbour array) are connected to existing hits
  if( cl->IsIterate() ) MoreNeighbours( cl );

  fHits[fNhits] = EBufferEnd;                  // mark no more hits
  // Normalise weighted mean, get fraction total energy in central crystal,
  // calc circular polar coordinates of cluster center
  *fMeanPosition = (*fMeanPosition) * (1./fSqrtEtot);// normalise weighted mean
  fCentralFrac = energy[fIndex]/fEnergy;      // fraction Etot in central elem
  fTheta = TMath::RadToDeg() * fMeanPosition->Theta();
  fPhi   = TMath::RadToDeg() * fMeanPosition->Phi();
}
    
//-----------------------------------------------------------------------------
void HitCluster_t::MoreNeighbours( TA2ClusterDetector* cl ){
  // Assume 1st go at cluster member search complete.
  // Now scan for any other close-proximity detector hits for potential
  // additional members of the cluster
  // Adapted from HitClusterTAPS_t (F.Zehr, Basle 2005)

  UInt_t i,j,k,n;
  Double_t distApart;                          // distance between elements
  TVector3** pos = cl->GetPosition();          // element position array
  UInt_t* hits = cl->GetTempHits();            // hits array
  UInt_t nhits = cl->GetNhits();               // # hits
  Double_t* energy = cl->GetEnergy();          // array of energies
  UInt_t* nextIter;                            // -> relevent part fHits array
  UInt_t* hitsEnd = fHits + fMaxHits - 1;      // end-stop of fHits array
  Double_t sqrtE;                              // square root energy weighting
  UInt_t pNhits = 0;                           // # previously processed hits

  // Iterate round while new cluster members still found
  do{
    nextIter = fHits + fNhits;                 // -> start new cluster members
    for( i=pNhits,n=0; i<fNhits; i++ ) {       // loop existing clust hits
      for ( j=0; j<nhits; j++ ) {              // loop total detector hits
	if ( (k = hits[j]) == ENullHit ) continue; // hit already spoken for
	// distance between cluster hit and potential new hit
	// if its within limits add the new hit to the cluster
	distApart = (*(pos[fHits[i]]) - *(pos[k])).Mag();
	if ( (distApart > cl->GetMinPosDiff()) &&
	     (distApart < cl->GetMaxPosDiff()) ) {
	  sqrtE = sqrt( energy[k] );
	  fEnergy += energy[ k];
	  fSqrtEtot += sqrtE;
	  *fMeanPosition += (*(pos[k])*sqrtE);// pos root energy weighted
	  hits[j] = ENullHit;                 // mark hit as spoken for
	  *nextIter++ = k;                    // add index to cluster hits list
	  n++;                                // update # new clust members
	  // Check if space for futher cluster members, finish if not
	  if( nextIter >= hitsEnd ){
	    fNhits += n;
	    return;
	  }
	} 
      }
    }
    pNhits = fNhits;         // update previously processed hits
    fNhits += n;             // update total processed hits
  }while( n );               // iterate while new cluster members found
  return;
}
