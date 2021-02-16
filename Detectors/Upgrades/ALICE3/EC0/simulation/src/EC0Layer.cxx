// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file EC0Layer.cxx
/// \brief Implementation of the EC0Layer class
/// \author Mario Sitta <sitta@to.infn.it>
/// \author Chinorat Kobdaj (kobdaj@g.sut.ac.th)

#include "EC0Simulation/EC0Layer.h"
#include "EC0Base/GeometryTGeo.h"
#include "EC0Simulation/Detector.h"

#include "FairLogger.h" // for LOG

#include <TGeoManager.h>        // for TGeoManager, gGeoManager
#include <TGeoMatrix.h>         // for TGeoCombiTrans, TGeoRotation, etc
#include <TGeoTube.h>           // for TGeoTube, TGeoTubeSeg
#include <TGeoVolume.h>         // for TGeoVolume, TGeoVolumeAssembly
#include <TGeoCompositeShape.h> // for TGeoCompositeShape
#include "TMathBase.h"          // for Abs
#include <TMath.h>              // for Sin, RadToDeg, DegToRad, Cos, Tan, etc

#include <cstdio> // for snprintf

class TGeoMedium;

using namespace TMath;
using namespace o2::ec0;
using namespace o2::itsmft;

ClassImp(EC0Layer);

EC0Layer::EC0Layer()
{
}

EC0Layer::~EC0Layer() = default;

EC0Layer::EC0Layer(Int_t layerNumber, std::string layerName, Float_t eta_in, Float_t eta_out, Float_t z, Float_t passive_x2X0)
{
  // Creates a simple parametrized EndCap layer covering the given
  // pseudorapidity range at the z layer position
  mLayerNumber = layerNumber;
  mLayerName = layerName;
  mZ = z;
  mPassive_x2X0 = passive_x2X0;
  mSensorThickness = 50 * 1.0E-4; // Micron
  mInnerRadius = std::abs(mZ * std::tan(2.f * std::atan(std::exp(-eta_in))));
  mOuterRadius = std::abs(mZ * std::tan(2.f * std::atan(std::exp(-eta_out))));
  mChipThickness = passive_x2X0; // TODO: calculate ChipThickness to match x/X0
}

void EC0Layer::createLayer(TGeoVolume* motherVolume)
{
  if (mLayerNumber >= 0) {
    // Create tube, set sensitive volume, add to mother volume

    std::string chipName = "EC0Chip_" + std::to_string(mLayerNumber),
                sensName = o2::ec0::GeometryTGeo::getEC0SensorPattern() + std::to_string(mLayerNumber);

    TGeoTube* sensor = new TGeoTube(mInnerRadius, mOuterRadius, mSensorThickness / 2);
    TGeoTube* chip = new TGeoTube(mInnerRadius, mOuterRadius, mChipThickness / 2);
    TGeoTube* layer = new TGeoTube(mInnerRadius, mOuterRadius, mChipThickness / 2);

    TGeoMedium* medSi = gGeoManager->GetMedium("EC0_SI$");
    TGeoMedium* medAir = gGeoManager->GetMedium("EC0_AIR$");

    TGeoVolume* sensVol = new TGeoVolume(sensName.c_str(), sensor, medSi);
    TGeoVolume* chipVol = new TGeoVolume(chipName.c_str(), chip, medAir);
    TGeoVolume* layerVol = new TGeoVolume(mLayerName.c_str(), layer, medAir);

    LOG(INFO) << "Inserting " << sensVol->GetName() << " inside " << chipVol->GetName();
    chipVol->AddNode(sensVol, 1, nullptr);

    LOG(INFO) << "Inserting " << chipVol->GetName() << " inside " << layerVol->GetName();
    layerVol->AddNode(chipVol, 0, nullptr);

    // Finally put everything in the mother volume
    auto* FwdDiskRotation = new TGeoRotation("FwdDiskRotation", 0, 0, 180);
    auto* FwdDiskCombiTrans = new TGeoCombiTrans(0, 0, mZ, FwdDiskRotation);

    LOG(INFO) << "Inserting " << layerVol->GetName() << " inside " << motherVolume->GetName();
    motherVolume->AddNode(layerVol, 1, FwdDiskCombiTrans);

    return;
  }
}
