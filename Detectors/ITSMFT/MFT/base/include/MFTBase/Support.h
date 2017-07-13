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
  TGeoMedium *mSupportMedium;

  Double_t mSupThickness; //Support Thickness
  Double_t mSupRad[5]; // Radius of each support disk
  Double_t mDiskGap; //gap between half disks
  Double_t mPhi0;
  Double_t mPhi1;
  Double_t mT_delta; //Excess to remove to avoid coplanar surfaces that causes visualization glitches
  Double_t mRaisedBoxHeight;

  ClassDefOverride(Support, 2);


};

}
}

#endif
