
#ifndef A2PolarizedTarget_h
#define A2PolarizedTarget_h 1

#include "A2Target.hh"
#include "A2MagneticField.hh"

class A2PolarizedTarget: public  A2Target
{

public:
  A2PolarizedTarget();
  ~A2PolarizedTarget();

  virtual G4VPhysicalVolume* Construct(G4LogicalVolume *MotherLogic);
  
  // Set magnetic field according to the field map
  virtual void SetMagneticField(G4String&);

private:
  A2MagneticField* fMagneticField;

};
#endif
