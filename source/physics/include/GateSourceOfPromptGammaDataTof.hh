/*----------------------
  GATE version name: gate_v6

  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

/*!
  \class  GateSourceOfPromptGammaDataTof

  Manage a 3D distribution of prompt gamma, with 1 time spectrum at
  each voxel, stored as TH1D. For voxels with yield==0, the TH1D is not
  instantiated/allocated.

*/

#ifndef GATEPROMPTGAMMASPATIALEMISSIONDISTRIBUTIONTOF_HH
#define GATEPROMPTGAMMASPATIALEMISSIONDISTRIBUTIONTOF_HH

#include "G4UnitsTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "GateConfiguration.h"
#include "GateImageOfHistograms.hh"
#include <random>
#include <iostream>
#include <fstream>

//------------------------------------------------------------------------
class GateSourceOfPromptGammaDataTof
{
public:
  GateSourceOfPromptGammaDataTof();
  ~GateSourceOfPromptGammaDataTof();

  void SampleRandomPositionToF(G4ThreeVector & position);
  void SampleRandomTime(double & time, int i, int j, int k);

  void LoadDataToF(std::string mFilename);
  void InitializeToF();
  double computesumtof;
  double ComputeSumToF() { return computesumtof; }

protected:
  // The 3D prompt gamma distribution
  GateImageOfHistograms * mImageTof; /** Modif Oreste **/
  std::vector<float> mDataCounts;

  //double * mean = new double[104125];
  //double * sigma = new double[104125];
  //std::vector<float> mean ;

  //std::vector<float> sigma;


  // Current pixel index for position in 3D space
  //int mCurrentIndex_i;
  //int mCurrentIndex_j;
  //int mCurrentIndex_k;

  // Physical coordinates from index coordinate
  std::vector<double> mIndexCoordX;
  std::vector<double> mIndexCoordY;
  std::vector<double> mIndexCoordZ;
  std::vector<TH1D*> mTimeGen; /** Modif Oreste **/

  G4SPSRandomGenerator mPositionXGenToF;
  std::vector<G4SPSRandomGenerator*> mPositionYGenToF;
  std::vector<std::vector<G4SPSRandomGenerator*> > mPositionZGenToF;

}; // end class
//------------------------------------------------------------------------

#endif /* end #define GATEPROMPTGAMMASPATIALEMISSIONDISTRIBUTION */
