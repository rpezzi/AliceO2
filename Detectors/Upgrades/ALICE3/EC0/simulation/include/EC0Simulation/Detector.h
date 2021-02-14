// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Detector.h
/// \brief Definition of the Detector class

#ifndef ALICEO2_EC0_DETECTOR_H_
#define ALICEO2_EC0_DETECTOR_H_

#include <vector>                             // for vector
#include "DetectorsBase/GeometryManager.h"    // for getSensID
#include "DetectorsBase/Detector.h"           // for Detector
#include "DetectorsCommonDataFormats/DetID.h" // for Detector
#include "ITSMFTSimulation/Hit.h"             // for Hit
#include "Rtypes.h"                           // for Int_t, Double_t, Float_t, Bool_t, etc
#include "TArrayD.h"                          // for TArrayD
#include "TGeoManager.h"                      // for gGeoManager, TGeoManager (ptr only)
#include "TLorentzVector.h"                   // for TLorentzVector
#include "TVector3.h"                         // for TVector3

class FairVolume;
class TGeoVolume;

class TParticle;

class TString;

namespace o2
{
namespace ec0
{
class GeometryTGeo;
}
} // namespace o2
namespace o2
{
namespace ec0
{
class V3Layer;
}
} // namespace o2

namespace o2
{
namespace ec0
{
class V3Layer;

class Detector : public o2::base::DetImpl<Detector>
{
 public:
  enum Model {
    kIBModelDummy = 0,
    kIBModel0 = 1,
    kIBModel1 = 2,
    kIBModel21 = 3,
    kIBModel22 = 4,
    kIBModel3 = 5,
    kIBModel4 = 10,
    kOBModelDummy = 6,
    kOBModel0 = 7,
    kOBModel1 = 8,
    kOBModel2 = 9
  };

  static constexpr Int_t sNumberLayers = 7;           ///< Number of layers in ITSU
  static constexpr Int_t sNumberInnerLayers = 3;      ///< Number of inner layers in ITSU
  static constexpr Int_t sNumberOfWrapperVolumes = 3; ///< Number of wrapper volumes

  /// Name : Detector Name
  /// Active: kTRUE for active detectors (ProcessHits() will be called)
  ///         kFALSE for inactive detectors
  Detector(Bool_t active);

  /// Default constructor
  Detector();

  /// Default destructor
  ~Detector() override;

  /// Initialization of the detector is done here
  void InitializeO2Detector() override;

  /// This method is called for each step during simulation (see FairMCApplication::Stepping())
  Bool_t ProcessHits(FairVolume* v = nullptr) override;

  /// Registers the produced collections in FAIRRootManager
  void Register() override;

  /// Gets the produced collections
  std::vector<o2::itsmft::Hit>* getHits(Int_t iColl) const
  {
    if (iColl == 0) {
      return mHits;
    }
    return nullptr;
  }

 public:
  /// Has to be called after each event to reset the containers
  void Reset() override;

  /// Base class to create the detector geometry
  void ConstructGeometry() override;

  /// This method is an example of how to add your own point of type Hit to the clones array
  o2::itsmft::Hit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                          const TVector3& startMom, double startE, double endTime, double eLoss,
                          unsigned char startStatus, unsigned char endStatus);

  /// Set per wrapper volume parameters
  void defineWrapperVolume(Int_t id, Double_t rmin, Double_t rmax, Double_t zspan) override;

  /// Add alignable top volumes
  void addAlignableVolumes() const override;

  void EndOfEvent() override;

  void FinishPrimary() override { ; }
  virtual void finishRun() { ; }
  void BeginPrimary() override { ; }
  void PostTrack() override { ; }
  void PreTrack() override { ; }
  /// Prints out the content of this class in ASCII format
  /// \param ostream *os The output stream
  void Print(std::ostream* os) const;

  /// Reads in the content of this class in the format of Print
  /// \param istream *is The input stream
  void Read(std::istream* is);

  /// Returns the number of layers
  Int_t getNumberOfLayers() const { return sNumberLayers; }

  GeometryTGeo* mGeometryTGeo; //! access to geometry details

 protected:
  Int_t mLayerID[sNumberLayers];     //! [sNumberLayers] layer identifier
  TString mLayerName[sNumberLayers]; //! [sNumberLayers] layer identifier

 private:
  /// this is transient data about track passing the sensor
  struct TrackData {               // this is transient
    bool mHitStarted;              //! hit creation started
    unsigned char mTrkStatusStart; //! track status flag
    TLorentzVector mPositionStart; //! position at entrance
    TLorentzVector mMomentumStart; //! momentum
    double mEnergyLoss;            //! energy loss
  } mTrackData;                    //!

  Int_t mNumberOfDetectors;

  Bool_t mModifyGeometry;

  Double_t mWrapperMinRadius[sNumberOfWrapperVolumes]; //! Min radius of wrapper volume
  Double_t mWrapperMaxRadius[sNumberOfWrapperVolumes]; //! Max radius of wrapper volume
  Double_t mWrapperZSpan[sNumberOfWrapperVolumes];     //! Z span of wrapper volume
  Int_t mWrapperLayerId[sNumberLayers];                //! Id of wrapper layer to which layer belongs (-1 if not wrapped)

  Bool_t mTurboLayer[sNumberLayers];          //! True for "turbo" layers
  Double_t mLayerPhi0[sNumberLayers];         //! Vector of layer's 1st stave phi in lab
  Double_t mLayerRadii[sNumberLayers];        //! Vector of layer radii
  Double_t mChipThickness[sNumberLayers];     //! Vector of chip thicknesses
  Double_t mDetectorThickness[sNumberLayers]; //! Vector of detector thicknesses
  UInt_t mChipTypeID[sNumberLayers];          //! Vector of detector type id
  Int_t mBuildLevel[sNumberLayers];           //! Vector of Material Budget Studies

  /// Container for hit data
  std::vector<o2::itsmft::Hit>* mHits;

  /// Create the detector materials
  virtual void createMaterials();

  /// Construct the detector geometry
  void constructDetectorGeometry();

  /// Define the sensitive volumes of the geometry
  void defineSensitiveVolumes();

  Detector(const Detector&);

  Detector& operator=(const Detector&);

  V3Layer* mGeometry[sNumberLayers]; //! Geometry

  template <typename Det>
  friend class o2::base::DetImpl;
  ClassDefOverride(Detector, 1);
};

// Input and output function for standard C++ input/output.
std::ostream& operator<<(std::ostream& os, Detector& source);

std::istream& operator>>(std::istream& os, Detector& source);
} // namespace ec0
} // namespace o2

#ifdef USESHM
namespace o2
{
namespace base
{
template <>
struct UseShm<o2::ec0::Detector> {
  static constexpr bool value = true;
};
} // namespace base
} // namespace o2
#endif

#endif
