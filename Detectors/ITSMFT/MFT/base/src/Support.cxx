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
mOuterCut{15.5,15.5,15.5,22.9,22.9},
mNumberOfBoxCuts{7,7,7,7,7},
mNumberOfRaixedBoxes{5,5,5,5,5},
mRaisedBoxHeight(0.305)
{

  initParameters();

}

//_____________________________________________________________________________
TGeoVolumeAssembly* Support::create(Int_t half, Int_t disk)
{

  Info("Create",Form("Creating Support_H%d_D%d", half,disk),0,0);

  mHalfDisk = new TGeoVolumeAssembly(Form("Support_H%d_D%d", half,disk));



  // ======= Support cuts (boxes) =========
  // TODO: Add real values for halfDisks 02 to 04
  Double_t BCuts[5][mNumberOfBoxCuts[disk]][4]={   // {Width, Height, x_center, y_center}
                            {{mSupRad[0]+mT_delta, mDiskGap, 0, 0}, // <-Disk0
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[0],2.)-pow(mOuterCut[0],2.)), (mSupRad[0]-mOuterCut[0])/2., 0., (mSupRad[0]+mOuterCut[0])/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[1]+mT_delta, mDiskGap, 0, 0},  // <-Disk1
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[1],2.)-pow(mOuterCut[1],2.)), (mSupRad[1]-mOuterCut[1])/2., 0., (mSupRad[1]+mOuterCut[1])/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[2]+mT_delta, mDiskGap, 0 ,0},  // <-Disk2
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[2],2.)-pow(mOuterCut[2],2.)), (mSupRad[2]-mOuterCut[2])/2., 0., (mSupRad[2]+mOuterCut[2])/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[3]+mT_delta, mDiskGap, 0, 0},  // <-Disk3
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[3],2.)-pow(mOuterCut[3],2.)), (mSupRad[1]-mOuterCut[3])/2., 0., (mSupRad[3]+mOuterCut[3])/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                            {{mSupRad[4]+mT_delta, mDiskGap, 0, 0},  // <-Disk4
                            {12.4, 6.91, 0, 0},
                            {7.95, 9.4, 0, 0},
                            {2.9, 11.885, 0, 0},
                            {sqrt(pow(mSupRad[4],2.)-pow(mOuterCut[4],2.)), (mSupRad[4]-mOuterCut[4])/2., 0., (mSupRad[4]+mOuterCut[4])/2.}, //External cut width: 2*sqrt(R²-x²)
                            {1.3875, 1.45, 16.1875, 7.9},
                            {1.3875, 1.45, -16.1875, 7.9}},
                          };

  auto *support = new TGeoVolumeAssembly(Form("Support_VA_H%d_D%d",half,disk));
  auto *base = new TGeoTubeSeg(Form("Base_H%d_D%d",half,disk), 0, mSupRad[disk], mSupThickness/2., mPhi0, mPhi1);

  for(Int_t cut=0 ; cut<mNumberOfBoxCuts[disk]; cut++){
    auto *boxName =  Form("BoxCut_%d_H%d_D%d",cut, half, disk);
    auto *boxCSName = Form("BoxCS_%d_H%d_D%d",cut, half, disk);
    mSomeBox = new TGeoBBox(boxName,BCuts[disk][cut][0],BCuts[disk][cut][1],  mSupThickness/2.+mT_delta);
    mSomeTranslation = new TGeoTranslation(BCuts[disk][cut][2],BCuts[disk][cut][3], 0.);
    //The first subtraction needs a shape, the base tube
    if (cut ==0)  mSomeSubtraction = new TGeoSubtraction(base, mSomeBox, NULL,mSomeTranslation);
      else    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeBox, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(boxCSName, mSomeSubtraction);
  }

  // ======= Support raised boxes =========
 // TODO: Add real values for halfDisks 02 to 04
  

  Double_t BRaised[5][mNumberOfRaixedBoxes[disk]][4]={  // {Width, Height, x_center, y_center}
                            {{(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 0
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}},
                            {{(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 1
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}},
                            {{(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 2
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}},
                            {{(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 3
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}},
                            {{(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 4
                            {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
                            {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
                            {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
                            {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}}
                          };

  for(Int_t box=0 ; box<mNumberOfRaixedBoxes[disk]; box++){
    mSomeBox = new TGeoBBox(NULL,BRaised[disk][box][0],BRaised[disk][box][1],  mRaisedBoxHeight/2.);
    mSomeTranslation = new TGeoTranslation(BRaised[disk][box][2],BRaised[disk][box][3], mRaisedBoxHeight/2.+mSupThickness/2.);
    mSomeUnion = new TGeoUnion(mSomeCS, mSomeBox, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeUnion);

    //For the backside
    mSomeTranslation = new TGeoTranslation(-BRaised[disk][box][2],BRaised[disk][box][3], -(mRaisedBoxHeight/2.+mSupThickness/2.));
    mSomeUnion = new TGeoUnion(mSomeCS, mSomeBox, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeUnion);
  }

  // =================  Holes ==================
// TODO: Add real values for halfDisks 02 to 04

  // ==== D6.5 H7 (6.5 mm diameter holes)
  Double_t D65 = .65/2.; //Radius
  auto *TubeD65 = new TGeoTube(Form("D65tube_H%d_D%d", half, disk),0, D65, mSupThickness/2.+2*mT_delta);
  Int_t nD65Holes[5]={2,2,2,2,2}; // Number of D6.5 mm Holes in each halfDisk support
  Double_t D65Holes[5][nD65Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-16.6, 13.5}, // #1     <- halfDisk 00
                            {16.6, 13.5}}, // #2
                            {{-16.6, 13.5}, // #1     <- halfDisk 01
                            {16.6, 13.5}}, // #2
                            {{-16.6, 14.9}, // #1    <- halfDisk 02
                            {16.6, 14.9}}, // #2
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };

  for(Int_t hole=0 ; hole<nD65Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D65Holes[disk][hole][0],
                               mOuterCut[disk]-D65Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD65, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D6 H7 (6 mm diameter holes)
  Double_t D6 = .6/2.; // Radius
  auto *TubeD6 = new TGeoTube(Form("D6tube_H%d_D%d", half, disk),0, D6, mSupThickness/2.+2*mT_delta);
  Int_t nD6Holes[5]={2,2,2,2,2}; // Number of D6 mm Holes in each halfDisk support
  Double_t D6Holes[5][nD6Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-16.6, 12.5}, // #3     <- halfDisk 00
                            {16.6, 12.5}}, // #4
                            {{-16.6, 12.5}, // #3     <- halfDisk 01
                            {16.6, 12.5}}, // #4
                            {{-16.6, 13.9}, // #3    <- halfDisk 02
                            {16.6, 13.9}}, // #4
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nD6Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D6Holes[disk][hole][0],
                               mOuterCut[disk]-D6Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD6, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }


  // ==== D8 H7 (8 mm diameter holes)
  Double_t D8 = .8/2.; // Radius
  auto *TubeD8 = new TGeoTube(Form("D8tube_H%d_D%d", half, disk),0, D8, mSupThickness/2.+2*mT_delta);
  Int_t nD8Holes[5]={2,2,2,2,2}; // Number of D8 mm Holes in each halfDisk support
  Double_t D8Holes[5][nD8Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-16.1, 10.3}, // #5     <- halfDisk 00
                            {16.1, 10.3}}, // #6
                            {{-16.1, 10.3}, // #5     <- halfDisk 01
                            {16.1, 10.3}}, // #6
                            {{-16.1, 11.7}, // #5    <- halfDisk 02
                            {16.1, 11.7}}, // #6
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nD8Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D8Holes[disk][hole][0],
                               mOuterCut[disk]-D8Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD8, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D3 H7 (3 mm diameter holes)
  Double_t D3 = .3/2.; // Radius
  auto *TubeD3 = new TGeoTube(Form("D3tube_H%d_D%d", half, disk),0, D3, mSupThickness/2.+2*mT_delta);
  Int_t nD3Holes[5]={2,2,2,2,2}; // Number of D3 mm Holes in each halfDisk support
  Double_t D3Holes[5][nD3Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-14.0, 6.0}, // #7     <- halfDisk 00
                            {14.0, 6.0}}, // #8
                            {{-14.0, 6.0}, // #7     <- halfDisk 01
                            {14.0, 6.0}}, // #8
                            {{-14.0, 7.4}, // #7    <- halfDisk 02
                            {14.0, 7.4}}, // #8
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nD3Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D3Holes[disk][hole][0],
                               mOuterCut[disk]-D3Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD3, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== M3 H7 (?? mm diameter holes)
  Double_t M3 = .3/2.; // Radius  TODO: Verify this!
  auto *TubeM3 = new TGeoTube(Form("M3tube_H%d_D%d", half, disk),0, M3, mSupThickness/2.+2*mT_delta);
  Int_t nM3Holes[5]={2,2,2,2,2}; // Number of M3 mm Holes in each halfDisk support
  Double_t M3Holes[5][nM3Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-11.2, 6.0}, // #11     <- halfDisk 00
                            {11.2, 6.0}}, // #12
                            {{-11.2, 6.0}, // #11     <- halfDisk 01
                            {11.2, 6.0}}, // #12
                            {{-12.0, 7.4}, // #7    <- halfDisk 02
                            {12.0, 7.4}}, // #8
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nM3Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-M3Holes[disk][hole][0],
                               mOuterCut[disk]-M3Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeM3, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // D4.5 H9
  Double_t D45 = .45/2.; // Radius
  auto *TubeD45 = new TGeoTube(Form("D45tube_H%d_D%d", half, disk),0, D45, mSupThickness/2.+2*mT_delta);
  Int_t nD45Holes[5]={2,2,2,2,2}; // Number of D45 mm Holes in each halfDisk support
  Double_t D45Holes[5][nD45Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-11.0, 4.0}, // #13     <- halfDisk 00
                            {11.0, 4.0}}, // #14
                            {{-11.0, 4.0}, // #13     <- halfDisk 01
                            {11.0, 4.0}}, // #14
                            {{-11.0, 5.4}, // #13    <- halfDisk 02
                            {11.0, 5.4}}, // #14
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nD45Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D45Holes[disk][hole][0],
                               mOuterCut[disk]-D45Holes[disk][hole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD45, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D2 H7 - 4 mm deep (on lower surface)
  Double_t D2 = .2/2., // Radius
           height_D2 = .4/2.; //
  auto *TubeD2 = new TGeoTube(Form("D2tube_H%d_D%d", half, disk),0, D2, mSupThickness/2.+2*mT_delta);
  Int_t nD2Holes[5]={2,2,2,2,2}; // Number of D2 mm Holes in each halfDisk support
  Double_t D2Holes[5][nD2Holes[disk]][2]={  // {x_center, y_center} // Blueprint REF.
                            {{-12.2, 8.295}, // #9     <- halfDisk 00
                            {12.2, 8.295}}, // #10
                            {{-12.2, 8.295}, // #9     <- halfDisk 01
                            {12.2, 8.295}}, // #10
                            {{-12.6, 9.695}, // #9    <- halfDisk 02
                            {12.6, 9.695}}, // #10
                            {{-16.6, 22.9}, // #     <- halfDisk 03
                            {16.6, 22.9}}, // #
                            {{-16.6, 22.9}, // #     <- halfDisk 04
                            {16.6, 22.9}} // #
                          };
  for(Int_t hole=0 ; hole<nD2Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D2Holes[disk][hole][0],
                               mOuterCut[disk]-D2Holes[disk][hole][1],
                               mSupThickness/2.-height_D2);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD2, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== M2 6mm deep (?? mm diameter holes)

  Double_t rad_M2 = .156657/2., // TODO: Check Radius of M2 holes
           height_M2 = .6/2;

  auto *TubeM2 = new TGeoTube(Form("sc_tube1_a_H%d_D%d", half, disk),0, rad_M2, height_M2+2.*mT_delta);

  //Int_t mNumberOfBoxCuts[5]={7,7,7,7,7};
  //Double_t BCuts[5][mNumberOfBoxCuts[disk]][4]

  Int_t nM2Holes[5]={12,12,12,12,12}; // Number of M2 Holes in each halfDisk support
  Double_t M2Holes[5][nM2Holes[disk]][2]={  // {x_center, y_center}
                            {{-8.75, 7.29}, // #16     <- halfDisk 00
                            {-7.05, 4.275}, // #18
                            {-5.35, 4.275}, // #20
                            {-3.65, 4.275}, // #22
                            {-1.95, 2.429}, // #24
                            {-0.25, 2.175}, // #26
                            {1.45, 2.269}, // #28
                            {3.15, 4.275}, // #30
                            {4.85, 4.275}, // #32
                            {6.55, 4.275}, // #34
                            {8.25, 7.29}, // #36
                            {9.95, 7.29}}, // #38
                            {{-8.75, 7.29}, // #16     <- halfDisk 01
                            {-7.05, 4.275}, // #18
                            {-5.53, 4.275}, // #20
                            {-3.65, 4.275}, // #22
                            {-1.95, 2.429}, // #24
                            {-0.25, 2.175}, // #26
                            {1.45, 2.269}, // #28
                            {3.15, 4.275}, // #30
                            {4.85, 4.275}, // #32
                            {6.55, 4.275}, // #34
                            {8.25, 7.29}, // #36
                            {9.95, 7.29}}, // #38
                            {{-8.75, 7.29}, // #     <- halfDisk 02
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}}, // #
                            {{-8.75, 7.29}, // #     <- halfDisk 03
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}}, // #
                            {{-8.75, 7.29}, // #     <- halfDisk 04
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}} // #
                          };

  for(Int_t hole=0 ; hole<nM2Holes[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-M2Holes[disk][hole][0],
                                mOuterCut[disk]-M2Holes[disk][hole][1],
                                mRaisedBoxHeight+mSupThickness/2.-height_M2);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeM2, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);

    //For the backside
    mSomeTranslation = new TGeoTranslation(M2Holes[disk][hole][0],mOuterCut[disk]-M2Holes[disk][hole][1], -(mRaisedBoxHeight+mSupThickness/2.-height_M2));
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeM2, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }


  // ==== D2 H7 - 4 mm deep (on higher surface)

  Double_t rad_D2_h = .2/2.,
           height_D2_h = .4/2;

  auto *TubeD2_h = new TGeoTube(Form("sc_tube1_a_H%d_D%d", half, disk),0, rad_D2_h, height_D2_h+2.*mT_delta);

  //Int_t mNumberOfBoxCuts[5]={7,7,7,7,7};
  //Double_t BCuts[5][mNumberOfBoxCuts[disk]][4]

  Int_t nD2_hHoles[5]={12,12,12,12,12}; // Number of D2_h Holes in each halfDisk support
  Double_t D2_hHoles[5][nD2_hHoles[disk]][2]={  // {x_center, y_center}
                            {{-8.75, 8.09}, // #15     <- halfDisk 00
                            {-7.05, 5.075}, // #17
                            {-5.35, 5.075}, // #19
                            {-3.65, 5.075}, // #21
                            {-1.95, 3.229}, // #23
                            {-0.25, 2.975}, // #25
                            {1.45, 3.069}, // #27
                            {3.15, 5.075}, // #29
                            {4.85, 5.075}, // #31
                            {6.55, 5.075}, // #33
                            {8.25, 8.09}, // #35
                            {9.95, 8.09}}, // #37
                            {{-8.75, 8.09}, // #15     <- halfDisk 01
                            {-7.05, 5.075}, // #17
                            {-5.35, 5.075}, // #19
                            {-3.65, 5.075}, // #21
                            {-1.95, 3.229}, // #23
                            {-0.25, 2.975}, // #25
                            {1.45, 3.069}, // #27
                            {3.15, 5.075}, // #29
                            {4.85, 5.075}, // #31
                            {6.55, 5.075}, // #33
                            {8.25, 8.09}, // #35
                            {9.95, 8.09}}, // #37
                            {{-8.75, 7.29}, // #     <- halfDisk 02
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}}, // #
                            {{-8.75, 7.29}, // #     <- halfDisk 03
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}}, // #
                            {{-8.75, 7.29}, // #     <- halfDisk 04
                            {-7.05, 4.275}, // #
                            {-5.53, 4.275}, // #
                            {-3.65, 4.275}, // #
                            {-1.95, 2.429}, // #
                            {-0.25, 2.175}, // #
                            {1.45, 2.269}, // #
                            {3.15, 4.275}, // #
                            {4.85, 4.275}, // #
                            {6.55, 4.275}, // #
                            {8.25, 7.29}, // #
                            {9.95, 7.29}} // #
                          };

  for(Int_t hole=0 ; hole<nD2_hHoles[disk]; hole++){
    mSomeTranslation = new TGeoTranslation(-D2_hHoles[disk][hole][0],
                                mOuterCut[disk]-D2_hHoles[disk][hole][1],
                                mRaisedBoxHeight+mSupThickness/2.-height_D2_h);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD2_h, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);

    //For the backside
    mSomeTranslation = new TGeoTranslation(D2_hHoles[disk][hole][0],mOuterCut[disk]-D2_hHoles[disk][hole][1], -(mRaisedBoxHeight+mSupThickness/2.-height_D2_h));
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, TubeD2_h, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }




  // ======= Prepare support volume and add to HalfDisk =========

  auto *support_vol = new TGeoVolume(Form("Support_H%d_D%d",half,disk), mSomeCS, mSupportMedium);

  mHalfDisk->AddNode(support_vol, 0);

  return mHalfDisk;

}

//_____________________________________________________________________________
void Support::initParameters()
{
 mSupportMedium  = gGeoManager->GetMedium("MFT_PEEK$");

 mSupportMedium  = gGeoManager->GetMedium("MFT_PEEK$");


}
