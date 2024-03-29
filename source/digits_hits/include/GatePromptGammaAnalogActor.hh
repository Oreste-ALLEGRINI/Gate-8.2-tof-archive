/*----------------------
  GATE version name: gate_v7
  Copyright (C): OpenGATE Collaboration
  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

#ifndef GATEPROMPTGAMMAPRODUCTIONANALOGACTOR_HH
#define GATEPROMPTGAMMAPRODUCTIONANALOGACTOR_HH

#include "GateConfiguration.h"
#include "GateVImageActor.hh"
#include "GateActorMessenger.hh"
#include "GatePromptGammaAnalogActorMessenger.hh"
#include "GateImageOfHistograms.hh"
#include "GatePromptGammaData.hh"


//-----------------------------------------------------------------------------
class GatePromptGammaAnalogActor: public GateVImageActor
{
public:
  virtual ~GatePromptGammaAnalogActor();

  FCT_FOR_AUTO_CREATOR_ACTOR(GatePromptGammaAnalogActor)

  virtual void Construct();
  virtual void UserPreTrackActionInVoxel(const int index, const G4Track* t);
  virtual void UserPostTrackActionInVoxel(const int index, const G4Track* t);
  virtual void UserSteppingActionInVoxel(const int index, const G4Step* step);
  virtual void BeginOfEventAction(const G4Event * e);/** Modif Oreste **/

  void SetInputDataFilename(std::string filename);
  virtual void SaveData();
  virtual void ResetData();

  void SetOutputCount(bool b) { mSetOutputCount = b; }  //output counts instead of yield

protected:
  GatePromptGammaAnalogActor(G4String name, G4int depth=0);
  GatePromptGammaAnalogActorMessenger * pMessenger;

  //we'll not use it, but extract the bins and binsizes for the PG output.
  std::string mInputDataFilename;
  GatePromptGammaData data;

  bool mSetOutputCount;
  bool alreadyHere;

  GateImageOfHistograms * mImageGamma;  //main output (yield)
  GateImageOfHistograms * mImagetof; /** Modif Oreste **/
  double startEvtTime;                  /** Modif Oreste **/
  TH1D * pTime = new TH1D("","",1000,0,5); /** Modif Oreste **/ //the source can be placed up to around 2.5 m upstream the target => Can be increased or decreased by modifying the TH size

};
//-----------------------------------------------------------------------------

MAKE_AUTO_CREATOR_ACTOR(PromptGammaAnalogActor,GatePromptGammaAnalogActor)

#endif // end GATEPROMPTGAMMAPRODUCTIONANALOGACTOR
