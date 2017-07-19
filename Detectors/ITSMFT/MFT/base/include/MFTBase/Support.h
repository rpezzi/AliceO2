// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Support.h
/// \brief MFT PCB support builder
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, Rafael Pezzi <rafael.pezzi@cern.ch>

#ifndef ALICEO2_MFT_SUPPORT_H_
#define ALICEO2_MFT_SUPPORT_H_

#include "TNamed.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoBoolNode.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "FairLogger.h"

namespace o2 {
namespace MFT {

class Support : public TNamed {

 public:

  Support();
  ~Support() override = default;
  TGeoVolumeAssembly* create(Int_t kHalf, Int_t disk);

 private:

  void initParameters();
  TGeoVolumeAssembly *mHalfDisk;
  TGeoMedium *mSupportMedium;  TGeoBBox *mSomeBox;
  TGeoSubtraction *mSomeSubtraction;
  TGeoUnion *mSomeUnion;
  TGeoTranslation *mSomeTranslation;
  TGeoCompositeShape *mSomeCS;

  Double_t mSupThickness; //Support Thickness
  Double_t mSupRad[5]; // Radius of each support disk
  Double_t mDiskGap; //gap between half disks
  Double_t mPhi0;
  Double_t mPhi1;
  Double_t mT_delta; //Excess to remove to avoid coplanar surfaces that causes visualization glitches
  Double_t mRaisedBoxHeight;
  Double_t mOuterCut[5]; //Distance of external disk cuts (oposite to beam pipe)
  Int_t mNumberOfBoxCuts[5]; // Number of box cuts in each half disk support
  Int_t mNumberOfRaixedBoxes[5]; //Number of Raised boxes in each halfDisk support



  ClassDefOverride(Support, 2);


};

}
}

#endif
