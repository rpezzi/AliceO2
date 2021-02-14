// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GeometryTGeo.h
/// \brief Definition of the GeometryTGeo class
/// \author cvetan.cheshkov@cern.ch - 15/02/2007
/// \author ruben.shahoyan@cern.ch - adapted to ITSupg 18/07/2012
/// \author rafael.pezzi@cern.ch - adapted to PostLS4EndCaps 25/06/2020

#ifndef ALICEO2_ENDCAPSLAYERS_GEOMETRYTGEO_H_
#define ALICEO2_ENDCAPSLAYERS_GEOMETRYTGEO_H_

#include <TGeoMatrix.h> // for TGeoHMatrix
#include <TObject.h>    // for TObject
#include <array>
#include <string>
#include <vector>
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "EndCapsBase/GeometryTGeo.h"
#include "MathUtils/Utils.h"
#include "Rtypes.h" // for Int_t, Double_t, Bool_t, UInt_t, etc

class TGeoPNEntry;

namespace o2
{
namespace ecl
{
/// GeometryTGeo is a simple interface class to TGeoManager. It is used in the simulation
/// in order to query the TGeo EC0 geometry.
/// RS: In order to preserve the static character of the class but make it dynamically access
/// geometry, we need to check in every method if the structures are initialized. To be converted
/// to singleton at later stage.

class GeometryTGeo : public o2::endcaps::GeometryTGeo
{
 public:
  typedef o2::math_utils::Transform3D Mat3D;
  using DetMatrixCache::getMatrixL2G;
  using DetMatrixCache::getMatrixT2GRot;
  using DetMatrixCache::getMatrixT2L;
  // this method is not advised for ITS: for barrel detectors whose tracking frame is just a rotation
  // it is cheaper to use T2GRot
  using DetMatrixCache::getMatrixT2G;

  static GeometryTGeo* Instance()
  {
    // get (create if needed) a unique instance of the object
    if (!sInstance)
      sInstance = std::unique_ptr<GeometryTGeo>(new GeometryTGeo(true, 0));
    return sInstance.get();
  }

  // adopt the unique instance from external raw pointer (to be used only to read saved instance from file)
  static void adopt(GeometryTGeo* raw);

  // constructor
  // ATTENTION: this class is supposed to behave as a singleton, but to make it root-persistent
  // we must define public default constructor.
  // NEVER use it, it will throw exception if the class instance was already created
  // Use GeometryTGeo::Instance() instead
  GeometryTGeo(bool build = kFALSE, int loadTrans = 0
               /*o2::base::utils::bit2Mask(o2::TransformType::T2L, // default transformations to load
           o2::TransformType::T2G,
           o2::TransformType::L2G)*/
  );

  /// Default destructor
  ~GeometryTGeo() override = default;

  GeometryTGeo(const GeometryTGeo& src) = delete;
  GeometryTGeo& operator=(const GeometryTGeo& geom) = delete;

  // implement filling of the matrix cache
  using o2::endcaps::GeometryTGeo::fillMatrixCache;
  void fillMatrixCache(int mask) override;

  /// Exract EC0 parameters from TGeo
  void Build(int loadTrans = 0) override;

  void Print(Option_t* opt = "") const;

 protected:
  static constexpr int MAXLAYERS = 15; ///< max number of active layers

  Int_t mNumberOfLayers;                 ///< number of layers
  static std::string sVolumeName;        ///< Mother volume name
  static std::string sLayerName;         ///< Layer name
  static std::string sWrapperVolumeName; ///< Wrapper volume name

 private:
  static std::unique_ptr<o2::ecl::GeometryTGeo> sInstance; ///< singletone instance

  ClassDefOverride(GeometryTGeo, 1); // EC0 geometry based on TGeo
};
} // namespace ecl
} // namespace o2

#endif
