// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PCBSupport.h
/// \brief Class describing geometry of one MFT half-disk PCBsupport
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>

#ifndef ALICEO2_MFT_PCBSUPPORT_H_
#define ALICEO2_MFT_PCBSUPPORT_H_

#include "TNamed.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoArb8.h"
#include "TGeoBoolNode.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "FairLogger.h"

class TGeoVolume;
class TGeoCompositeShape;

namespace o2
{
namespace MFT
{

class PCBSupport : public TNamed
{

 public:
  PCBSupport();
  ~PCBSupport() override = default;
  TGeoVolumeAssembly* create(Int_t kHalf, Int_t disk);

 private:

  void initParameters();
  TGeoVolumeAssembly *mHalfDisk;
  TGeoMedium *mSupportMedium;
  TGeoBBox *mSomeBox;
  TGeoTube *mSomeTube;
  TGeoArb8 *mSomeArb;

  TGeoSubtraction *mSomeSubtraction;
  TGeoUnion *mSomeUnion;
  TGeoTranslation *mSomeTranslation;
  TGeoCompositeShape *mSomeCS;

  Double_t mPCBThickness; //Support Thickness
  Double_t mPCBRad[5]; // Radius of each PCB disk
  Double_t mDiskGap; //gap between half disks
  Double_t mPhi0;
  Double_t mPhi1;
  Double_t mT_delta; //Excess to remove to avoid coplanar surfaces that causes visualization glitches
  Double_t mOuterCut[5]; //Distance of external disk cuts (oposite to beam pipe)
                         // this is the y origin on Guillamet's PDF blueprints

  Int_t mNumberOfBoxCuts[5]; // Number of box cuts in each half disk support
  Double_t (*mBoxCuts[5])[4];// Box cuts on each disk

  ClassDefOverride(PCBSupport, 2);

};
}
}

#endif
