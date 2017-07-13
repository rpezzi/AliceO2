// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Support.cxx
/// \brief Class building the MFT PCB Supports
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, Rafael Pezzi <rafael.pezzi@cern.ch>

#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoBoolNode.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"

#include "FairLogger.h"

#include "MFTBase/Constants.h"
#include "MFTBase/Support.h"
#include "MFTBase/Geometry.h"

using namespace o2::MFT;

/// \cond CLASSIMP
ClassImp(o2::MFT::Support)
/// \endcond

//_____________________________________________________________________________
Support::Support() :
TNamed(),
mHalfDisk(nullptr),
mDiskGap(1.4),
mSupRad{17.5,17.5,17.5,25.5,25.5},
mSupThickness(.8),
mPhi0(0.),
mPhi1(180.),
mT_delta(0.00001),
mRaisedBoxHeight(0.305)
{

  initParameters();

}

//_____________________________________________________________________________
TGeoVolumeAssembly* Support::create(Int_t half, Int_t disk)
{

  Info("Create",Form("Creating Support_H%d_D%d", half,disk),0,0);

  mHalfDisk = new TGeoVolumeAssembly(Form("Support_H%d_D%d", half,disk));

  TGeoBBox *someBox;
  TGeoSubtraction *baseCut;
  TGeoUnion *baseUnion;
  TGeoCompositeShape *base_CS;
  TGeoTranslation *tCut;

  // ======= Support cuts (boxes) =========

  Int_t nBoxCuts[5]={7,7,7,7,7};
  Double_t BCuts[5][nBoxCuts[disk]][4]={   // {Width, Height, x_center, y_center}
                            {{mSupRad[0]+mT_delta, mDiskGap, 0, 0}, // <-Disk0
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[0],2)-15.5*15.5), (mSupRad[0]-15.5)/2., 0., (mSupRad[0]+15.5)/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[1]+mT_delta, mDiskGap, 0, 0},  // <-Disk1
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[1],2)-15.5*15.5), (mSupRad[1]-15.5)/2., 0., (mSupRad[1]+15.5)/2.},
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[2]+mT_delta, mDiskGap, 0 ,0},  // <-Disk2
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[2],2)-16.9*16.9), (mSupRad[2]-16.9)/2., 0., (mSupRad[2]+16.9)/2.},
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[3]+mT_delta, mDiskGap, 0, 0},  // <-Disk3
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[3],2)-22.9*22.9), (mSupRad[3]-22.9)/2., 0., (mSupRad[3]+22.9)/2.}, // Test values
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[4]+mT_delta, mDiskGap, 0, 0},  // <-Disk4
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[4],2)-22.9*22.9), (mSupRad[4]-22.9)/2., 0., (mSupRad[4]+22.9)/2.}, // Test values
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                          };

  auto *support = new TGeoVolumeAssembly(Form("Support_VA_H%d_D%d",half,disk));
  auto *base = new TGeoTubeSeg(Form("Base_H%d_D%d",half,disk), 0, mSupRad[disk], mSupThickness/2., mPhi0, mPhi1);

  for(Int_t cut=0 ; cut<nBoxCuts[disk]; cut++){
    auto *boxName =  Form("BoxCut_%d_H%d_D%d",cut, half, disk);
    auto *boxCSName = Form("BoxCS_%d_H%d_D%d",cut, half, disk);
    someBox = new TGeoBBox(boxName,BCuts[disk][cut][0],BCuts[disk][cut][1],  mSupThickness/2.+mT_delta);
    tCut = new TGeoTranslation(BCuts[disk][cut][2],BCuts[disk][cut][3], 0.);
    //The first subtraction needs a shape, the base tube
    if (cut ==0)  baseCut = new TGeoSubtraction(base, someBox, NULL,tCut);
      else    baseCut = new TGeoSubtraction(base_CS, someBox, NULL,tCut);
    base_CS = new TGeoCompositeShape(boxCSName, baseCut);
  }

  // ======= Support raised boxes =========
  // TODO: implement for each disk

  Int_t nRaisedBox=5;
  Double_t BRaised[nRaisedBox][4]={  // {Width, Height, x_center, y_center}
                            {(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.},
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}
                          };

  for(Int_t box=0 ; box<nRaisedBox; box++){
    someBox = new TGeoBBox(NULL,BRaised[box][0],BRaised[box][1],  mRaisedBoxHeight/2.);
    tCut = new TGeoTranslation(BRaised[box][2],BRaised[box][3], mRaisedBoxHeight/2.+mSupThickness/2.);
    baseUnion = new TGeoUnion(base_CS, someBox, NULL,tCut);
    base_CS = new TGeoCompositeShape(NULL, baseUnion);

    //For the backside
    tCut = new TGeoTranslation(-BRaised[box][2],BRaised[box][3], -(mRaisedBoxHeight/2.+mSupThickness/2.));
    baseUnion = new TGeoUnion(base_CS, someBox, 0,tCut);
    base_CS = new TGeoCompositeShape(NULL, baseUnion);
  }

  // =================  Holes ==================
  // TODO: Implement holes and opennings for each support


  // ======= Prepare support volume and add to HalfDisk =========

  auto *support_vol = new TGeoVolume(Form("Support_H%d_D%d",half,disk), base_CS, mSupportMedium);

  mHalfDisk->AddNode(support_vol, 0);

  return mHalfDisk;

}

//_____________________________________________________________________________
void Support::initParameters()
{
  mSupportMedium  = gGeoManager->GetMedium("MFT_PEEK$");


}
