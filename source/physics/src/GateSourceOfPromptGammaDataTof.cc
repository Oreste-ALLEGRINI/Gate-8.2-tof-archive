/*----------------------
  GATE version name: gate_v6

  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

#include "GateSourceOfPromptGammaDataTof.hh"
#include "GateMessageManager.hh"
#include "Randomize.hh" // needed for G4UniformRand
#include "G4Gamma.hh"
#include "GateRandomEngine.hh"
#include <random>
#include <iostream>
#include <fstream>

//------------------------------------------------------------------------
GateSourceOfPromptGammaDataTof::GateSourceOfPromptGammaDataTof()
{
    computesumtof = 0;
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
GateSourceOfPromptGammaDataTof::~GateSourceOfPromptGammaDataTof()
{
  for(unsigned int i=0; i< mPositionYGenToF.size(); i++)
    delete mPositionYGenToF[i];
  for(unsigned int i=0; i< mPositionZGenToF.size(); i++)
    for(unsigned int j=0; i< mPositionZGenToF[i].size(); j++)
      delete mPositionZGenToF[i][j];
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
void GateSourceOfPromptGammaDataTof::LoadDataToF(std::string mFilename)
{
  mImageTof = new GateImageOfHistograms("float");
  mImageTof->Read(G4String(removeExtension(mFilename))+"-tof."+G4String(getExtension(mFilename)));
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
void GateSourceOfPromptGammaDataTof::InitializeToF()
{
  unsigned int sizeX = mImageTof->GetResolution().x();
  unsigned int sizeY = mImageTof->GetResolution().y();
  unsigned int sizeZ = mImageTof->GetResolution().z();
  unsigned int nbOfBinsTof = mImageTof->GetNbOfBins();
  unsigned long nbOfValues = sizeX*sizeY*sizeZ;

  // Random generators for position: set sizes
  mPositionYGenToF.resize(sizeX);
  mPositionZGenToF.resize(sizeX);
  for(unsigned int i=0; i<sizeX; i++) {
    mPositionZGenToF[i].resize(sizeY);
  }

  // Build the scalar image with total number of counts at each pixel
  mImageTof->ComputeTotalOfCountsImageDataFloat(mDataCounts);
  computesumtof = mImageTof->ComputeSum();

  // Initialize random generator for position. Loop over the total
  // count scalar image (mDataCounts)
  mPositionXGenToF.SetXBias(G4ThreeVector(0., 0., 0.)); // important
  for(unsigned int i=0; i<sizeX; i++) {
    double sumYZ = 0.0;
    mPositionYGenToF[i] = new G4SPSRandomGenerator;
    mPositionYGenToF[i]->SetYBias(G4ThreeVector(0., 0., 0.)); // important
    for(unsigned int j=0; j<sizeY; j++) {
      double sumZ = 0.0;
      mPositionZGenToF[i][j] = new G4SPSRandomGenerator;
      mPositionZGenToF[i][j]->SetZBias(G4ThreeVector(0., 0., 0.)); // important
      for(unsigned int k=0; k<sizeZ; k++) {
        double val = mDataCounts[mImageTof->GetIndexFromPixelIndex(i, j, k)];
        sumZ += val;
        // Bias the Z component according to the voxel value
        mPositionZGenToF[i][j]->SetZBias(G4ThreeVector(k+1 ,val,0.));
      }
      sumYZ += sumZ;
      // Bias the Y component according to integration over Z
      mPositionYGenToF[i]->SetYBias(G4ThreeVector(j+1, sumZ,0.));
    }
    // Bias the X component according to integration over YZ plane
    mPositionXGenToF.SetXBias(G4ThreeVector(i+1, sumYZ,0.));
  }

  // Initialize time.
  mTimeGen.resize(nbOfValues);
  double timeStep  = (mImageTof->GetMaxValue()-mImageTof->GetMinValue())/nbOfBinsTof;
  double time = 0.0;
  long index_image_time = 0;
  long index_data_time = 0;
  float * data_time = mImageTof->GetDataFloatPointer();
  long nbNonZero_time = 0;

  // We only create TH1D for non zero pixel.
  for(unsigned int k=0; k<sizeZ; k++) {
    for(unsigned int j=0; j<sizeY; j++) {
      for(unsigned int i=0; i<sizeX; i++) {
        if (mDataCounts[index_image_time] == 0) { // FIXME
          index_data_time+=nbOfBinsTof;
        }
        else {
          time = mImageTof->GetMinValue();
          mTimeGen[index_image_time] = new TH1D;
          // This is much much faster to use the constructor without
          // param, the SetBins than using the following line with constructor :
          // new TH1D("", "", nbOfBins, mImage->GetMinValue(), mImage->GetMaxValue());
          TH1D * htime = mTimeGen[index_image_time];
          htime->SetBins(nbOfBinsTof, mImageTof->GetMinValue(), mImageTof->GetMaxValue());
          for(unsigned int l=0; l<nbOfBinsTof; l++) {
            htime->Fill(time, data_time[index_data_time]);
            index_data_time++;
            time += timeStep;
          }
          nbNonZero_time++;
        }
        index_image_time++;
      }
    }
  }

  // ATTENTION: THIS DELETES THE ON DISK DATA FROM MEMORY. ACCESSING THE DATA IN mImage
  // WILL SEGFAULT. METADATA IS KEPT.
  mImageTof->Deallocate();

}
//------------------------------------------------------------------------
void GateSourceOfPromptGammaDataTof::SampleRandomPositionToF(G4ThreeVector & position)
{
  // Random 3D position (in pixel). If size == 1, then bug in GenRand
  // (infinite loop), so we test and set random value [0:1]
  double x;
  if (mImageTof->GetResolution().x() == 1) x=G4UniformRand();
  else x = mPositionXGenToF.GenRandX();
  int i =  mCurrentIndex_i = floor(x);
  double y;
  if (mImageTof->GetResolution().y() == 1) y=G4UniformRand();
  else y = mPositionYGenToF[i]->GenRandY();
  int j = mCurrentIndex_j = floor(y);
  double z;
  if (mImageTof->GetResolution().z() == 1) z=G4UniformRand();
  else  z = mPositionZGenToF[i][j]->GenRandZ();
  mCurrentIndex_k = floor(z);

  // Offset according to image origin (and half voxel position)
  x = mImageTof->GetOrigin().x() + x*mImageTof->GetVoxelSize().x();
  y = mImageTof->GetOrigin().y() + y*mImageTof->GetVoxelSize().y();
  z = mImageTof->GetOrigin().z() + z*mImageTof->GetVoxelSize().z();

  // Vector
  position.setX(x);
  position.setY(y);
  position.setZ(z);
}

//------------------------------------------------------------------------

double GateSourceOfPromptGammaDataTof::SampleRandomTime(double & time)
{
    // Get time spectrum in the current pixel
    long index = mImageTof->GetIndexFromPixelIndex(mCurrentIndex_i, mCurrentIndex_j, mCurrentIndex_k);
    if (mDataCounts[index] != 0) {
      time = mTimeGen[index]->GetRandom();
      std::cout<<"Time = "<<time<<std::endl;
    }
    else time = 0.0;
}
