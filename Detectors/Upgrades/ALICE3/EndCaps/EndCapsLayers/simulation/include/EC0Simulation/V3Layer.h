// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file V3Layer.h
/// \brief Definition of the V3Layer class
/// \author Mario Sitta <sitta@to.infn.it>
/// \author Chinorat Kobdaj (kobdaj@g.sut.ac.th)

#ifndef ALICEO2_ENDCAPSLAYERS_UPGRADEV3LAYER_H_
#define ALICEO2_ENDCAPSLAYERS_UPGRADEV3LAYER_H_

#include <TGeoManager.h>               // for gGeoManager
#include "Rtypes.h"                    // for Double_t, Int_t, Bool_t, etc
#include "EC0Simulation/V11Geometry.h" // for V11Geometry
#include "EC0Simulation/Detector.h"    // for Detector, Detector::Model

class TGeoXtru;

class TGeoCombiTrans;

class TGeoVolume; // lines 15-15

namespace o2
{
namespace ecl
{

/// This class defines the Geometry for the ITS  using TGeo. This is a work class used
/// to study different configurations during the development of the new ITS structure
class V3Layer : public V11Geometry
{

  // Default constructor
  V3Layer();

  /// Constructor setting layer number and debugging level
  /// for a "turbo" layer (i.e. where staves overlap in phi)
  V3Layer(Int_t lay, Bool_t turbo = kFALSE, Int_t debug = 0);

  /// Copy constructor
  V3Layer(const V3Layer&) = default;

  /// Assignment operator
  V3Layer& operator=(const V3Layer&) = default;

  /// Default destructor
  ~V3Layer() override;

  Bool_t isTurbo() const { return mIsTurbo; };

  /// Creates the actual Layer and places inside its mother volume
  /// \param motherVolume the TGeoVolume owing the volume structure
  virtual void createLayer(TGeoVolume* motherVolume);

 private:
  /// Creates the actual Layer and places inside its mother volume
  /// A so-called "turbo" layer is a layer where staves overlap in phi
  /// User can set width and tilt angle, no check is performed here
  /// to avoid volume overlaps
  /// \param motherVolume The TGeoVolume owing the volume structure
  void createLayerTurbo(TGeoVolume* motherVolume);

  /// Computes the inner radius of the air container for the Turbo configuration
  /// as the radius of either the circle tangent to the stave or the circle
  /// passing for the stave's lower vertex. Returns the radius of the container
  /// if >0, else flag to use the lower vertex

  /// Help method to create a TGeoCombiTrans matrix from a similar method with same name and
  /// function in V11GeometrySDD class by L.Gaudichet)
  /// Returns the TGeoCombiTrans which make a translation in y and z and a rotation in phi
  /// in the global coord system. If planeSym = true, the rotation places the object
  /// symetrically (with respect to the transverse plane) to its position in the
  /// case planeSym = false
  TGeoCombiTrans* createCombiTrans(const char* name, Double_t dy, Double_t dz, Double_t dphi, Bool_t planeSym = kFALSE);

  /// Help method to add a translation to a TGeoCombiTrans matrix (from a similar method
  /// with same name and function in V11GeometrySDD class by L.Gaudichet)
  void addTranslationToCombiTrans(TGeoCombiTrans* ct, Double_t dx = 0, Double_t dy = 0, Double_t dz = 0) const;

  Int_t mLayerNumber;        ///< Current layer number
  Double_t mPhi0;            ///< lab phi of 1st stave, in degrees!!!
  Double_t mLayerRadius;     ///< Inner radius of this layer
  Double_t mSensorThickness; ///< Sensor thickness
  Double_t mChipThickness;   ///< Chip thickness
  ///< container)
  Int_t mNumberOfChips; ///< Number chips per container (module, HalfStave, Stave, whatever is
  /// container)

  UInt_t mChipTypeID; ///< detector type id
  Bool_t mIsTurbo;    ///< True if this layer is a "turbo" layer
  Int_t mBuildLevel;  ///< Used for material studies

  // Parameters for the  geometry

  ClassDefOverride(V3Layer, 0); // ITS v3 geometry
};
} // namespace ecl
} // namespace o2

#endif
