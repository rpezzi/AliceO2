// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file V3Layer.cxx
/// \brief Implementation of the V3Layer class
/// \author Mario Sitta <sitta@to.infn.it>
/// \author Chinorat Kobdaj (kobdaj@g.sut.ac.th)

#include "EC0Simulation/V3Layer.h"
#include "ECLayersBase/GeometryTGeo.h"
#include "EC0Simulation/Detector.h"

#include "FairLogger.h" // for LOG

#include <TGeoArb8.h>           // for TGeoArb8
#include <TGeoBBox.h>           // for TGeoBBox
#include <TGeoCone.h>           // for TGeoConeSeg, TGeoCone
#include <TGeoPcon.h>           // for TGeoPcon
#include <TGeoManager.h>        // for TGeoManager, gGeoManager
#include <TGeoMatrix.h>         // for TGeoCombiTrans, TGeoRotation, etc
#include <TGeoTrd1.h>           // for TGeoTrd1
#include <TGeoTube.h>           // for TGeoTube, TGeoTubeSeg
#include <TGeoVolume.h>         // for TGeoVolume, TGeoVolumeAssembly
#include <TGeoXtru.h>           // for TGeoXtru
#include <TGeoCompositeShape.h> // for TGeoCompositeShape
#include "TMathBase.h"          // for Abs
#include <TMath.h>              // for Sin, RadToDeg, DegToRad, Cos, Tan, etc

#include <cstdio> // for snprintf

class TGeoMedium;

using namespace TMath;
using namespace o2::ecl;
using namespace o2::endcaps;

ClassImp(V3Layer);

#define SQ(A) (A) * (A)

V3Layer::V3Layer()
  : V11Geometry()
{
}

V3Layer::~V3Layer() = default;

void V3Layer::createLayer(TGeoVolume* motherVolume)
{
}

void V3Layer::createLayerTurbo(TGeoVolume* motherVolume)
{
}

TGeoCombiTrans* V3Layer::createCombiTrans(const char* name, Double_t dy, Double_t dz, Double_t dphi, Bool_t planeSym)
{
  TGeoTranslation t1(dy * cosD(90. + dphi), dy * sinD(90. + dphi), dz);
  TGeoRotation r1("", 0., 0., dphi);
  TGeoRotation r2("", 90, 180, -90 - dphi);

  TGeoCombiTrans* combiTrans1 = new TGeoCombiTrans(name);
  combiTrans1->SetTranslation(t1);
  if (planeSym) {
    combiTrans1->SetRotation(r1);
  } else {
    combiTrans1->SetRotation(r2);
  }
  return combiTrans1;
}

void V3Layer::addTranslationToCombiTrans(TGeoCombiTrans* ct, Double_t dx, Double_t dy, Double_t dz) const
{
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  const Double_t* vect = ct->GetTranslation();
  Double_t newVect[3] = {vect[0] + dx, vect[1] + dy, vect[2] + dz};
  ct->SetTranslation(newVect);
}
