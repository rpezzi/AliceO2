//Info// Copyright CERN and copyright holders of ALICE O2. This software is
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
mRaisedBoxHeight(0.305),
mRad_M2(.156657/2.), // TODO: Check Radius of M2 holes
mHeight_M2(.6/2), // Height of M2 holes on raised boxes
mRad_D2_h(.2/2.),
mHeight_D2_h(.4/2),
mTwoHoles(2), // Number of D6.5 mm Holes in each halfDisk support
mD65(.65/2.), //Radius
mD6(.6/2.), // Radius
mD8(.8/2.), // Radius
mD3(.3/2.), // Radius
mM3(.3/2.), // Radius   TODO: Verify this!
mD45(.45/2.), // Radius
mD2(.2/2.), // Radius
mHeight_D2(.4/2.)
{

  initParameters();

}

//_____________________________________________________________________________
TGeoVolumeAssembly* Support::create(Int_t half, Int_t disk)
{

  Info("Create",Form("Creating Support_H%d_D%d", half,disk),0,0);
  mHalfDisk = new TGeoVolumeAssembly(Form("Support_H%d_D%d", half,disk));
  auto *support = new TGeoVolumeAssembly(Form("Support_VA_H%d_D%d",half,disk));
  auto *base = new TGeoTubeSeg(Form("Base_H%d_D%d",half,disk), 0, mSupRad[disk], mSupThickness/2., mPhi0, mPhi1);

  // Cutting boxes
  //Info("Create",Form("Cutting Boxes Support_H%d_D%d", half,disk),0,0);
  for(Int_t cut = 0 ; cut<mNumberOfBoxCuts[disk]; cut++){
    auto *boxName =  Form("BoxCut_%d_H%d_D%d",cut, half, disk);
    auto *boxCSName = Form("BoxCS_%d_H%d_D%d",cut, half, disk);
    mSomeBox = new TGeoBBox(boxName,mBoxCuts[disk][cut][0],mBoxCuts[disk][cut][1],  mSupThickness/2.+mT_delta);
    mSomeTranslation = new TGeoTranslation(mBoxCuts[disk][cut][2],mBoxCuts[disk][cut][3], 0.);
    //The first subtraction needs a shape, the base tube
    if (cut ==0)  mSomeSubtraction = new TGeoSubtraction(base, mSomeBox, NULL,mSomeTranslation);
      else    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeBox, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(boxCSName, mSomeSubtraction);
  }

  // Adding raisedBoxes
  //Info("Create",Form("Adding raised boxes Support_H%d_D%d", half,disk),0,0);
  for(Int_t box=0 ; box<mNumberOfRaixedBoxes[disk]; box++){
    mSomeBox = new TGeoBBox(NULL,mBRaised[disk][box][0], mBRaised[disk][box][1], mRaisedBoxHeight/2.);
    mSomeTranslation = new TGeoTranslation(mBRaised[disk][box][2],
      mBRaised[disk][box][3],
      mRaisedBoxHeight/2.+mSupThickness/2.);
    mSomeUnion = new TGeoUnion(mSomeCS, mSomeBox, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeUnion);
      //For the backside
    mSomeTranslation = new TGeoTranslation(-mBRaised[disk][box][2],
      mBRaised[disk][box][3],
      -(mRaisedBoxHeight/2.+mSupThickness/2.));
    mSomeUnion = new TGeoUnion(mSomeCS, mSomeBox, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeUnion);
  }

  // =================  Holes ==================

  // ==== M2 6mm deep holes)
  //Info("Create",Form("Cutting M2 6 mm deep holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("sc_tube1_a_H%d_D%d", half, disk),0, mRad_M2, mHeight_M2+2.*mT_delta);
  for(Int_t iHole=0 ; iHole<mNumberOfM2Holes[disk]; iHole++){
    mSomeTranslation = new TGeoTranslation(-mM2Holes[disk][iHole][0],
                                mOuterCut[disk]-mM2Holes[disk][iHole][1],
                                mRaisedBoxHeight+mSupThickness/2.-mHeight_M2);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);

    //For the backside
    mSomeTranslation = new TGeoTranslation(mM2Holes[disk][iHole][0],
                                mOuterCut[disk]-mM2Holes[disk][iHole][1],
                                -(mRaisedBoxHeight+mSupThickness/2.-mHeight_M2));
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D2 H7 - 4 mm deep (on raisedBoxes)
  //Info("Create",Form("Cutting D2 mm holes on raisedboxes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("sc_tube1_a_H%d_D%d", half, disk),0, mRad_D2_h, mHeight_D2_h+2.*mT_delta);

  for(Int_t iHole=0 ; iHole<mNumberOfD2_hHoles[disk]; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD2_hHoles[disk][iHole][0],
                                mOuterCut[disk]-mD2_hHoles[disk][iHole][1],
                                mRaisedBoxHeight+mSupThickness/2.-mHeight_D2_h);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);

    //For the backside
    mSomeTranslation = new TGeoTranslation(mD2_hHoles[disk][iHole][0],
                                mOuterCut[disk]-mD2_hHoles[disk][iHole][1],
                                -(mRaisedBoxHeight+mSupThickness/2.-mHeight_D2_h));
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, 0,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D6.5 H7 (6.5 mm diameter holes)
  //Info("Create",Form("Cutting 6.5 holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D65tube_H%d_D%d", half, disk),0, mD65, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole= 0 ; iHole < mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD65Holes[disk][iHole][0],
                               mOuterCut[disk]-mD65Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL, mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D6 H7 (6 mm diameter holes)
  //Info("Create",Form("Cutting 6 mm holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D6tube_H%d_D%d", half, disk),0, mD6, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD6Holes[disk][iHole][0],
                               mOuterCut[disk]-mD6Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D8 H7 (8 mm diameter holes)
  //Info("Create",Form("Cutting 8 mm holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D8tube_H%d_D%d", half, disk),0, mD8, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD8Holes[disk][iHole][0],
                               mOuterCut[disk]-mD8Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D3 H7 (3 mm diameter holes)
  //Info("Create",Form("Cutting 3 mm holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D3tube_H%d_D%d", half, disk),0, mD3, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD3Holes[disk][iHole][0],
                               mOuterCut[disk]-mD3Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== M3 H7 (?? mm diameter holes)
  //Info("Create",Form("Cutting M3 H7 holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("M3tube_H%d_D%d", half, disk),0, mM3, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mM3Holes[disk][iHole][0],
                               mOuterCut[disk]-mM3Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D4.5 H9
  //Info("Create",Form("Cutting 4.5 mm holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D45tube_H%d_D%d", half, disk),0, mD45, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD45Holes[disk][iHole][0],
                               mOuterCut[disk]-mD45Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ==== D2 H7 - 4 mm deep (on lower surface)
  //Info("Create",Form("Cutting D2 holes Support_H%d_D%d", half,disk),0,0);
  mSomeTube = new TGeoTube(Form("D2tube_H%d_D%d", half, disk),0, mD2, mSupThickness/2.+2*mT_delta);
  for(Int_t iHole=0 ; iHole<mTwoHoles; iHole++){
    mSomeTranslation = new TGeoTranslation(-mD2Holes[disk][iHole][0],
                               mOuterCut[disk]-mD2Holes[disk][iHole][1],
                               0.);
    mSomeSubtraction = new TGeoSubtraction(mSomeCS, mSomeTube, NULL,mSomeTranslation);
    mSomeCS = new TGeoCompositeShape(NULL, mSomeSubtraction);
  }

  // ======= Prepare support volume and add to HalfDisk =========

  auto *support_vol = new TGeoVolume(Form("Support_H%d_D%d",half,disk), mSomeCS, mSupportMedium);
  //mHalfDisk->AddNode(support_vol, 0);

  if (disk==0 && half==0) { // Adding test.gdml in place of support 0, half 0
    TGDMLParse parser;
    TGeoVolume* gdmlInput;
    gdmlInput = parser.GDMLReadFile("test.gdml"); // test.gdml at current/working runtime directory
    gdmlInput->SetMedium(mSupportMedium);
    mHalfDisk->AddNode(gdmlInput,1);
  }
    else mHalfDisk->AddNode(support_vol, 0);

  return mHalfDisk;

}

//_____________________________________________________________________________
void Support::initParameters()
{

  mSupportMedium  = gGeoManager->GetMedium("MFT_PEEK$");

  // # Support parametrization =====
  // TODO: Add real values for halfDisks 02 to 04

  // ================================================
  // ## Cut boxes (squares)
  // ### halfDisks 00
  mNumberOfBoxCuts[0]=7;
  // Cut boxes {Width, Height, x_center, y_center}
  mBoxCuts[00] = new Double_t[mNumberOfBoxCuts[0]][4]{
    {mSupRad[0]+mT_delta, mDiskGap, 0, 0},
    {sqrt(pow(mSupRad[0],2.)-pow(mOuterCut[0],2.)),
      (mSupRad[0]-mOuterCut[0])/2.,
      0.,
      (mSupRad[0]+mOuterCut[0])/2.}, //External cut width: 2*sqrt(R²-x²)
    {12.4, 6.91, 0, 0},
    {7.95, 9.4, 0, 0},
    {2.9, 11.885, 0, 0},
    {1.3875, 1.45, 16.1875, 7.9},
    {1.3875, 1.45, -16.1875, 7.9}
  };

  // ### halfDisks 01
  mNumberOfBoxCuts[1]=mNumberOfBoxCuts[0];
  mBoxCuts[01]=mBoxCuts[00];

  // ### halfDisk 02
  mNumberOfBoxCuts[2]=7;
  mBoxCuts[02] = new Double_t[mNumberOfBoxCuts[2]][4]{
    {mSupRad[2]+mT_delta, mDiskGap, 0, 0},
    {sqrt(pow(mSupRad[2],2.)-pow(mOuterCut[2],2.)),
      (mSupRad[2]-mOuterCut[2])/2.,
      0.,
      (mSupRad[2]+mOuterCut[2])/2.}, //External cut width: 2*sqrt(R²-x²)
    {12.4, 6.91, 0, 0},
    {7.95, 9.4, 0, 0},
    {2.9, 11.885, 0, 0},
    {1.3875, 1.45, 16.1875, 7.9},
    {1.3875, 1.45, -16.1875, 7.9}
  };

  // ### halfDisk 03
  mNumberOfBoxCuts[3]=7;
  mBoxCuts[03] = new Double_t[mNumberOfBoxCuts[3]][4]{
    {mSupRad[3]+mT_delta, mDiskGap, 0, 0},
    {sqrt(pow(mSupRad[3],2.)-pow(mOuterCut[3],2.)),
      (mSupRad[3]-mOuterCut[3])/2.,
      0.,
      (mSupRad[3]+mOuterCut[3])/2.}, //External cut width: 2*sqrt(R²-x²)
    {12.4, 6.91, 0, 0},
    {7.95, 9.4, 0, 0},
    {2.9, 11.885, 0, 0},
    {1.3875, 1.45, 16.1875, 7.9},
    {1.3875, 1.45, -16.1875, 7.9}
  };

  // ### halfDisk 04
  mNumberOfBoxCuts[4]=7;
  mBoxCuts[04] = new Double_t[mNumberOfBoxCuts[4]][4]{
    {mSupRad[4]+mT_delta, mDiskGap, 0, 0},
    {sqrt(pow(mSupRad[4],2.)-pow(mOuterCut[4],2.)),
      (mSupRad[4]-mOuterCut[4])/2.,
      0.,
      (mSupRad[4]+mOuterCut[4])/2.}, //External cut width: 2*sqrt(R²-x²)
    {12.4, 6.91, 0, 0},
    {7.95, 9.4, 0, 0},
    {2.9, 11.885, 0, 0},
    {1.3875, 1.45, 16.1875, 7.9},
    {1.3875, 1.45, -16.1875, 7.9}
  };

  // ================================================
  // ## Raised boxes
  // ### halfDisks 00
  mNumberOfRaixedBoxes[0]=5;
  // Raised Boxes {Width, Height, x_center, y_center}
  mBRaised[0] = new Double_t[mNumberOfRaixedBoxes[0]][4]{
    {(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.},
    {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
    {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
    {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
    {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}
  };

  // ### halfDisks 01
  mNumberOfRaixedBoxes[1]=mNumberOfRaixedBoxes[0];
  mBRaised[1]=mBRaised[0];

  // ### halfDisk 02
  mNumberOfRaixedBoxes[2]=5;
  mBRaised[02] = new Double_t[mNumberOfRaixedBoxes[2]][4]{
    {(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.},
    {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
    {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
    {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
    {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}
  };

  // ### halfDisk 03
  mNumberOfRaixedBoxes[3]=5;
  mBRaised[03] = new Double_t[mNumberOfRaixedBoxes[3]][4]{
    {(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.},
    {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
    {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
    {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
    {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}
  };

  // ### halfDisk 04
  mNumberOfRaixedBoxes[4]=5;
  mBRaised[04] = new Double_t[mNumberOfRaixedBoxes[4]][4] {
    {(9.35-7.95)/2.,(8.81-6.91)/2.,(9.35+7.95)/2.,(8.81+6.91)/2.}, // <- halfDisk 4
    {(7.65-2.9)/2.,(11.82-9.4)/2.,(7.65+2.9)/2.,(11.82+9.4)/2.},
    {(2.55+2.05)/2.,(13.92-11.885)/2.,(2.55-2.05)/2.,(13.92+11.885)/2.},
    {(7.15-2.9)/2.,(11.82-9.4)/2.,(-7.152-2.92)/2.,(11.82+9.4)/2},
    {(10.55-7.95)/2.,(8.81-6.91)/2.,(-10.55-7.95)/2.,(8.81+6.91)/2.}
  };

  // ================================================
  // ## M2 6mm deep
  // ### halfDisk00
  mNumberOfM2Holes[0] = 12;
  mM2Holes[0] = new Double_t[mNumberOfM2Holes[0]][2]{
    {-8.75, 7.29}, // #16
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
    {9.95, 7.29} // #38
  };

  // ### halfDisk01
  mNumberOfM2Holes[1] = mNumberOfM2Holes[0];
  mM2Holes[1]=mM2Holes[0];

  // ### halfDisk02
  mNumberOfM2Holes[2]=12;
  mM2Holes[2] = new Double_t[mNumberOfM2Holes[2]][2]{
    {-8.75, 7.29},
    {-7.05, 4.275},
    {-5.35, 4.275},
    {-3.65, 4.275},
    {-1.95, 2.429},
    {-0.25, 2.175},
    {1.45, 2.269},
    {3.15, 4.275},
    {4.85, 4.275},
    {6.55, 4.275},
    {8.25, 7.29},
    {9.95, 7.29}
  };

  // ### halfDisk03
  mNumberOfM2Holes[3]=12;
  mM2Holes[3] = new Double_t[mNumberOfM2Holes[3]][2] {
    {-8.75, 7.29}, //
    {-7.05, 4.275}, //
    {-5.35, 4.275}, //
    {-3.65, 4.275}, //
    {-1.95, 2.429}, //
    {-0.25, 2.175}, //
    {1.45, 2.269}, //
    {3.15, 4.275}, //
    {4.85, 4.275}, //
    {6.55, 4.275}, //
    {8.25, 7.29}, //
    {9.95, 7.29} //
  };

  // ### halfDisk04
  mNumberOfM2Holes[4]=12;
  mM2Holes[4] = new Double_t[mNumberOfM2Holes[4]][2] {
    {-8.75, 7.29}, //
    {-7.05, 4.275}, //
    {-5.35, 4.275}, //
    {-3.65, 4.275}, //
    {-1.95, 2.429}, //
    {-0.25, 2.175}, //
    {1.45, 2.269}, //
    {3.15, 4.275}, //
    {4.85, 4.275}, //
    {6.55, 4.275}, //
    {8.25, 7.29}, //
    {9.95, 7.29} //
  };

  // ================================================
  // ## D2 H7 - 4 mm deep (on raisedBoxes)
  // ### halfDisk00
  mNumberOfD2_hHoles[0] = 12;
  mD2_hHoles[0] = new Double_t[mNumberOfD2_hHoles[0]][2]
    {{-8.75, 8.09}, // #15
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
    {9.95, 8.09}
  };

  // ### halfDisk01
  mNumberOfD2_hHoles[1] = mNumberOfD2_hHoles[0];
  mD2_hHoles[1]=mD2_hHoles[0];

  // ### halfDisk02
  mNumberOfD2_hHoles[2] = 12;
  mD2_hHoles[2] = new Double_t[mNumberOfD2_hHoles[2]][2]{
    {-8.75, 7.29}, // #
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
    {9.95, 7.29}
  };

  // ### halfDisk03
  mNumberOfD2_hHoles[3] = 12;
  mD2_hHoles[3] = new Double_t[mNumberOfD2_hHoles[3]][2]{
    {-8.75, 7.29},  // #
    {-7.05, 4.275}, // #
    {-5.53, 4.275}, // #
    {-3.65, 4.275}, // #
    {-1.95, 2.429}, // #
    {-0.25, 2.175}, // #
    {1.45, 2.269},  // #
    {3.15, 4.275},  // #
    {4.85, 4.275},  // #
    {6.55, 4.275},  // #
    {8.25, 7.29},   // #
    {9.95, 7.29}    // #
  };

  // ### halfDisk04
  mNumberOfD2_hHoles[4] = 12;
  mD2_hHoles[4] = new Double_t[mNumberOfD2_hHoles[4]][2]{
    {-8.75, 7.29}, // #
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
    {9.95, 7.29} // #
  };

  // ================================================
  // ## D6.5 mm holes
  // ### halfDisk00
  mD65Holes[0] = new Double_t[mTwoHoles][2]{
    {-16.6, 13.5},
    {16.6, 13.5}
  };

  // ### halfDisk01
  mD65Holes[1] = mD65Holes[0];

  // ### halfDisk02
  mD65Holes[2] = new Double_t[mTwoHoles][2]{
    {-16.6, 14.9}, // #1
    {16.6, 14.9}   // #2
  };

  // ### halfDisk03
  mD65Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9}   // #
  };

  // ### halfDisk02
  mD65Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #     <- halfDisk 04
    {16.6, 22.9}    // #2
  };

  // ================================================
  // ## D6 mm holes
  // ### halfDisk00
  mD6Holes[0] = new Double_t[mTwoHoles][2]{
    {-16.6, 12.5}, // #3
    {16.6, 12.5} // #4
  };

  // ### halfDisk01
  mD6Holes[1] = mD6Holes[0];


  // ### halfDisk02
  mD6Holes[2] = new Double_t[mTwoHoles][2]{
    {-16.6, 13.9}, // #3
    {16.6, 13.9} // #4
  };

  // ### halfDisk03
  mD6Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mD6Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ================================================
  // ## D8 mm holes
  // ### halfDisk00
  mD8Holes[0] = new Double_t[mTwoHoles][2]{
    {-16.1, 10.3}, // #5
    {16.1, 10.3} // #6
  };

  // ### halfDisk01
  mD8Holes[1] = mD8Holes[0];

  // ### halfDisk02
  mD8Holes[2] = new Double_t[mTwoHoles][2]{
    {-16.1, 11.7}, // #5
    {16.1, 11.7} // #6
  };

  // ### halfDisk03
  mD8Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mD8Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ================================================
  // ## D3 mm holes
  // ### halfDisk00
  mD3Holes[0] = new Double_t[mTwoHoles][2]{
    {-14.0, 6.0}, // #7
    {14.0, 6.0}   // #8
  };

  // ### halfDisk01
  mD3Holes[1] = mD3Holes[0];

  // ### halfDisk02
  mD3Holes[2] = new Double_t[mTwoHoles][2]{
    {-14.0, 7.4}, // #7
    {14.0, 7.4} // #8
  };

  // ### halfDisk03
  mD3Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mD3Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ================================================
  // ## M3 H7 mm holes
  // ### halfDisk00
  mM3Holes[0] = new Double_t[mTwoHoles][2]{
    {-11.2, 6.0}, // #11
    {11.2, 6.0} // #12
  };

  // ### halfDisk01
  mM3Holes[1] = mM3Holes[0];

  // ### halfDisk02
  mM3Holes[2] = new Double_t[mTwoHoles][2]{
    {-12.0, 7.4}, // #7
    {12.0, 7.4} // #8
  };

  // ### halfDisk03
  mM3Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mM3Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ================================================
  // ## D45 H7 mm holes
  // ### halfDisk00
  mD45Holes[0] = new Double_t[mTwoHoles][2]{
    {-11.0, 4.0}, // #13
    {11.0, 4.0} // #14
  };

  // ### halfDisk01
  mD45Holes[1] = mD45Holes[0];

  // ### halfDisk02
  mD45Holes[2] = new Double_t[mTwoHoles][2]{
    {-11.0, 5.4}, // #13
    {11.0, 5.4} // #14
  };

  // ### halfDisk03
  mD45Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mD45Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ================================================
  // ## D2 H7 mm holes - 4 mm deep (on lower surface)
  // ### halfDisk00
  mD2Holes[0] = new Double_t[mTwoHoles][2]{
    {-12.2, 8.295}, // #9
    {12.2, 8.295} // #10
  };

  // ### halfDisk01
  mD2Holes[1] = mD2Holes[0];

  // ### halfDisk02
  mD2Holes[2] = new Double_t[mTwoHoles][2]{
    {-12.6, 9.695}, // #9
    {12.6, 9.695} // #10
  };

  // ### halfDisk03
  mD2Holes[3] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

  // ### halfDisk04
  mD2Holes[4] = new Double_t[mTwoHoles][2]{
    {-16.6, 22.9}, // #
    {16.6, 22.9} // #
  };

}
