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

void EC0Layer::createLayer(TGeoVolume* motherVolume)
{
}
