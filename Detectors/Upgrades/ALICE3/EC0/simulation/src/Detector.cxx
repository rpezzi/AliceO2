// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Detector.cxx
/// \brief Implementation of the Detector class

#include "ITSMFTSimulation/Hit.h"
#include "EC0Base/GeometryTGeo.h"
#include "EC0Simulation/Detector.h"
#include "EC0Simulation/V3Layer.h"

#include "SimulationDataFormat/Stack.h"
#include "SimulationDataFormat/TrackReference.h"

// FairRoot includes
#include "FairDetector.h"    // for FairDetector
#include "FairLogger.h"      // for LOG, LOG_IF
#include "FairRootManager.h" // for FairRootManager
#include "FairRun.h"         // for FairRun
#include "FairRuntimeDb.h"   // for FairRuntimeDb
#include "FairVolume.h"      // for FairVolume
#include "FairRootManager.h"

#include "TGeoManager.h"     // for TGeoManager, gGeoManager
#include "TGeoTube.h"        // for TGeoTube
#include "TGeoPcon.h"        // for TGeoPcon
#include "TGeoVolume.h"      // for TGeoVolume, TGeoVolumeAssembly
#include "TString.h"         // for TString, operator+
#include "TVirtualMC.h"      // for gMC, TVirtualMC
#include "TVirtualMCStack.h" // for TVirtualMCStack

#include <cstdio> // for NULL, snprintf

class FairModule;

class TGeoMedium;

class TParticle;

using std::cout;
using std::endl;

using namespace o2::ec0;
using o2::itsmft::Hit;

Detector::Detector()
  : o2::base::DetImpl<Detector>("EC0", kTRUE),
    mTrackData(),
    /*
    mHitStarted(false),
    mTrkStatusStart(),
    mPositionStart(),
    mMomentumStart(),
    mEnergyLoss(),
    */
    mNumberOfDetectors(-1),
    mModifyGeometry(kFALSE),
    mHits(o2::utils::createSimVector<Hit>())
{
}

static double radii2Turbo(double rMin, double rMid, double rMax, double sensW)
{
  // compute turbo angle from radii and sensor width
  return TMath::ASin((rMax * rMax - rMin * rMin) / (2 * rMid * sensW)) * TMath::RadToDeg();
}

static void configEC0(Detector* ec0)
{
  // build EC0 upgrade detector
}

Detector::Detector(Bool_t active)
  : o2::base::DetImpl<Detector>("EC0", active),
    mTrackData(),
    mNumberOfDetectors(-1),
    mModifyGeometry(kFALSE),
    mHits(o2::utils::createSimVector<Hit>())
{

  configEC0(this);
}

Detector::Detector(const Detector& rhs)
  : o2::base::DetImpl<Detector>(rhs),
    mTrackData(),
    /*
    mHitStarted(false),
    mTrkStatusStart(),
    mPositionStart(),
    mMomentumStart(),
    mEnergyLoss(),
    */
    mNumberOfDetectors(rhs.mNumberOfDetectors),
    mModifyGeometry(rhs.mModifyGeometry),

    /// Container for data points
    mHits(o2::utils::createSimVector<Hit>())
{
}

Detector::~Detector()
{

  if (mHits) {
    // delete mHits;
    o2::utils::freeSimVector(mHits);
  }
}

Detector& Detector::operator=(const Detector& rhs)
{
  // The standard = operator
  // Inputs:
  //   Detector   &h the sourse of this copy
  // Outputs:
  //   none.
  // Return:
  //  A copy of the sourse hit h

  if (this == &rhs) {
    return *this;
  }

  // base class assignment
  base::Detector::operator=(rhs);

  mNumberOfDetectors = rhs.mNumberOfDetectors;

  mModifyGeometry = rhs.mModifyGeometry;

  /// Container for data points
  mHits = nullptr;

  return *this;
}

void Detector::InitializeO2Detector()
{
  // Define the list of sensitive volumes
  defineSensitiveVolumes();

  for (int i = 0; i < sNumberLayers; i++) {
    mLayerID[i] = gMC ? TVirtualMC::GetMC()->VolId(mLayerName[i]) : 0;
  }

  mGeometryTGeo = GeometryTGeo::Instance();
  //  FairRuntimeDb* rtdb= FairRun::Instance()->GetRuntimeDb();
  //  O2itsGeoPar* par=(O2itsGeoPar*)(rtdb->getContainer("O2itsGeoPar"));
}

Bool_t Detector::ProcessHits(FairVolume* vol)
{
  // This method is called from the MC stepping
  if (!(fMC->TrackCharge())) {
    return kFALSE;
  }

  Int_t lay = 0, volID = vol->getMCid();

  // FIXME: Determine the layer number. Is this information available directly from the FairVolume?
  bool notSens = false;
  while ((lay < sNumberLayers) && (notSens = (volID != mLayerID[lay]))) {
    ++lay;
  }
  if (notSens)
    return kFALSE; // RS: can this happen? This method must be called for sensors only?

  // Is it needed to keep a track reference when the outer EC0 volume is encountered?
  auto stack = (o2::data::Stack*)fMC->GetStack();
  if (fMC->IsTrackExiting() && (lay == 0 || lay == 6)) {
    // Keep the track refs for the innermost and outermost layers only
    o2::TrackReference tr(*fMC, GetDetId());
    tr.setTrackID(stack->GetCurrentTrackNumber());
    tr.setUserId(lay);
    stack->addTrackReference(tr);
  }
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (fMC->IsTrackEntering()) {
    status |= Hit::kTrackEntering;
  }
  if (fMC->IsTrackInside()) {
    status |= Hit::kTrackInside;
  }
  if (fMC->IsTrackExiting()) {
    status |= Hit::kTrackExiting;
  }
  if (fMC->IsTrackOut()) {
    status |= Hit::kTrackOut;
  }
  if (fMC->IsTrackStop()) {
    status |= Hit::kTrackStopped;
  }
  if (fMC->IsTrackAlive()) {
    status |= Hit::kTrackAlive;
  }

  // track is entering or created in the volume
  if ((status & Hit::kTrackEntering) || (status & Hit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (Hit::kTrackExiting | Hit::kTrackOut | Hit::kTrackStopped))) {
    stopHit = true;
  }

  // increment energy loss at all steps except entrance
  if (!startHit)
    mTrackData.mEnergyLoss += fMC->Edep();
  if (!(startHit | stopHit))
    return kFALSE; // do noting

  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    fMC->TrackMomentum(mTrackData.mMomentumStart);
    fMC->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop;
    fMC->TrackPosition(positionStop);
    // Retrieve the indices with the volume path
    int stave(0), halfstave(0), chipinmodule(0), module(0);
    int chipindex = 0;

    Hit* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                    mTrackData.mMomentumStart.Vect(), mTrackData.mMomentumStart.E(), positionStop.T(),
                    mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    // p->SetTotalEnergy(vmc->Etot());

    // RS: not sure this is needed
    // Increment number of Detector det points in TParticle
    stack->addHit(GetDetId());
  }

  return kTRUE;
}

void Detector::createMaterials()
{
  Int_t ifield = 2;
  Float_t fieldm = 10.0;
  o2::base::Detector::initFieldTrackingParams(ifield, fieldm);
  ////////////

  Float_t tmaxfd = 0.1;   // 1.0; // Degree
  Float_t stemax = 1.0;   // cm
  Float_t deemax = 0.1;   // 30.0; // Fraction of particle's energy 0<deemax<=1
  Float_t epsil = 1.0E-4; // 1.0; // cm
  Float_t stmin = 0.0;    // cm "Default value used"

  Float_t tmaxfdSi = 0.1;    // .10000E+01; // Degree
  Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
  Float_t deemaxSi = 0.1;    // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
  Float_t epsilSi = 1.0E-4;  // .10000E+01;
  Float_t stminSi = 0.0;     // cm "Default value used"

  Float_t tmaxfdAir = 0.1;        // .10000E+01; // Degree
  Float_t stemaxAir = .10000E+01; // cm
  Float_t deemaxAir = 0.1;        // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
  Float_t epsilAir = 1.0E-4;      // .10000E+01;
  Float_t stminAir = 0.0;         // cm "Default value used"

  // AIR
  Float_t aAir[4] = {12.0107, 14.0067, 15.9994, 39.948};
  Float_t zAir[4] = {6., 7., 8., 18.};
  Float_t wAir[4] = {0.000124, 0.755267, 0.231781, 0.012827};
  Float_t dAir = 1.20479E-3;

  // Water
  Float_t aWater[2] = {1.00794, 15.9994};
  Float_t zWater[2] = {1., 8.};
  Float_t wWater[2] = {0.111894, 0.888106};
  Float_t dWater = 1.0;

  // PEEK CF30
  Float_t aPEEK[3] = {12.0107, 1.00794, 15.9994};
  Float_t zPEEK[3] = {6., 1., 8.};
  Float_t wPEEK[3] = {19., 12., 3};
  Float_t dPEEK = 1.32;

  // Kapton
  Float_t aKapton[4] = {1.00794, 12.0107, 14.010, 15.9994};
  Float_t zKapton[4] = {1., 6., 7., 8.};
  Float_t wKapton[4] = {0.026362, 0.69113, 0.07327, 0.209235};
  Float_t dKapton = 1.42;

  // Tungsten Carbide
  Float_t aWC[2] = {183.84, 12.0107};
  Float_t zWC[2] = {74, 6};
  Float_t wWC[2] = {0.5, 0.5};
  Float_t dWC = 15.63;

  // BEOL (Metal interconnection stack in Si sensors)
  Float_t aBEOL[3] = {26.982, 28.086, 15.999};
  Float_t zBEOL[3] = {13, 14, 8}; // Al, Si, O
  Float_t wBEOL[3] = {0.170, 0.388, 0.442};
  Float_t dBEOL = 2.28;

  // Inox 304
  Float_t aInox304[4] = {12.0107, 51.9961, 58.6928, 55.845};
  Float_t zInox304[4] = {6., 24., 28, 26};       // C, Cr, Ni, Fe
  Float_t wInox304[4] = {0.0003, 0.18, 0.10, 0}; // [3] will be computed
  Float_t dInox304 = 7.85;

  // Ceramic (for IB capacitors) (BaTiO3)
  Float_t aCeramic[3] = {137.327, 47.867, 15.999};
  Float_t zCeramic[3] = {56, 22, 8}; // Ba, Ti, O
  Float_t wCeramic[3] = {1, 1, 3};   // Molecular composition
  Float_t dCeramic = 6.02;

  // Rohacell (C9 H13 N1 O2)
  Float_t aRohac[4] = {12.01, 1.01, 14.010, 16.};
  Float_t zRohac[4] = {6., 1., 7., 8.};
  Float_t wRohac[4] = {9., 13., 1., 2.};
  Float_t dRohac = 0.05;

  o2::base::Detector::Mixture(1, "AIR$", aAir, zAir, dAir, 4, wAir);
  o2::base::Detector::Medium(1, "AIR$", 1, 0, ifield, fieldm, tmaxfdAir, stemaxAir, deemaxAir, epsilAir, stminAir);

  o2::base::Detector::Mixture(2, "WATER$", aWater, zWater, dWater, 2, wWater);
  o2::base::Detector::Medium(2, "WATER$", 2, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  o2::base::Detector::Material(3, "SI$", 0.28086E+02, 0.14000E+02, 0.23300E+01, 0.93600E+01, 0.99900E+03);
  o2::base::Detector::Medium(3, "SI$", 3, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  o2::base::Detector::Material(4, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7); // From AliPIPEv3
  o2::base::Detector::Medium(4, "BERILLIUM$", 4, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  o2::base::Detector::Material(5, "COPPER$", 0.63546E+02, 0.29000E+02, 0.89600E+01, 0.14300E+01, 0.99900E+03);
  o2::base::Detector::Medium(5, "COPPER$", 5, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // needed for STAVE , Carbon, kapton, Epoxy, flexcable

  // AliceO2::Base::Detector::Material(6,"CARBON$",12.0107,6,2.210,999,999);
  o2::base::Detector::Material(6, "CARBON$", 12.0107, 6, 2.210 / 1.3, 999, 999);
  o2::base::Detector::Medium(6, "CARBON$", 6, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  o2::base::Detector::Mixture(7, "KAPTON(POLYCH2)$", aKapton, zKapton, dKapton, 4, wKapton);
  o2::base::Detector::Medium(7, "KAPTON(POLYCH2)$", 7, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // values below modified as compared to source AliEC0v11 !

  // BEOL (Metal interconnection stack in Si sensors)
  o2::base::Detector::Mixture(29, "METALSTACK$", aBEOL, zBEOL, dBEOL, 3, wBEOL);
  o2::base::Detector::Medium(29, "METALSTACK$", 29, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // Glue between IB chip and FPC: density reduced to take into account
  // empty spaces (160 glue spots/chip , diam. 1 spot = 1 mm)
  o2::base::Detector::Material(30, "GLUE_IBFPC$", 12.011, 6, 1.05 * 0.3, 999, 999);
  o2::base::Detector::Medium(30, "GLUE_IBFPC$", 30, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // Ceramic for IB capacitors (nmat < 0 => wmat contains number of atoms)
  o2::base::Detector::Mixture(31, "CERAMIC$", aCeramic, zCeramic, dCeramic, -3, wCeramic);
  o2::base::Detector::Medium(31, "CERAMIC$", 31, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // All types of carbon
  // Unidirectional prepreg
  o2::base::Detector::Material(8, "K13D2U2k$", 12.0107, 6, 1.643, 999, 999);
  o2::base::Detector::Medium(8, "K13D2U2k$", 8, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  o2::base::Detector::Material(17, "K13D2U120$", 12.0107, 6, 1.583, 999, 999);
  o2::base::Detector::Medium(17, "K13D2U120$", 17, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // Carbon prepreg woven
  o2::base::Detector::Material(18, "F6151B05M$", 12.0107, 6, 2.133, 999, 999);
  o2::base::Detector::Medium(18, "F6151B05M$", 18, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // Impregnated thread
  o2::base::Detector::Material(9, "M60J3K$", 12.0107, 6, 2.21, 999, 999);
  o2::base::Detector::Medium(9, "M60J3K$", 9, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // Impregnated thread
  o2::base::Detector::Material(10, "M55J6K$", 12.0107, 6, 1.63, 999, 999);
  o2::base::Detector::Medium(10, "M55J6K$", 10, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // Fabric(0/90)
  o2::base::Detector::Material(11, "T300$", 12.0107, 6, 1.725, 999, 999);
  o2::base::Detector::Medium(11, "T300$", 11, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // AMEC Thermasol
  o2::base::Detector::Material(12, "FGS003$", 12.0107, 6, 1.6, 999, 999);
  o2::base::Detector::Medium(12, "FGS003$", 12, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
  // Carbon fleece
  o2::base::Detector::Material(13, "CarbonFleece$", 12.0107, 6, 0.4, 999, 999);
  o2::base::Detector::Medium(13, "CarbonFleece$", 13, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi,
                             stminSi);
  // Rohacell
  o2::base::Detector::Mixture(32, "ROHACELL$", aRohac, zRohac, dRohac, -4, wRohac);
  o2::base::Detector::Medium(32, "ROHACELL$", 32, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  // PEEK CF30
  o2::base::Detector::Mixture(19, "PEEKCF30$", aPEEK, zPEEK, dPEEK, -3, wPEEK);
  o2::base::Detector::Medium(19, "PEEKCF30$", 19, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  // Flex cable
  Float_t aFCm[5] = {12.0107, 1.00794, 14.0067, 15.9994, 26.981538};
  Float_t zFCm[5] = {6., 1., 7., 8., 13.};
  Float_t wFCm[5] = {0.520088819984, 0.01983871336, 0.0551367996, 0.157399667056, 0.247536};
  // Float_t dFCm = 1.6087;  // original
  // Float_t dFCm = 2.55;   // conform with STAR
  Float_t dFCm = 2.595; // conform with Corrado

  o2::base::Detector::Mixture(14, "FLEXCABLE$", aFCm, zFCm, dFCm, 5, wFCm);
  o2::base::Detector::Medium(14, "FLEXCABLE$", 14, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  // AliceO2::Base::Detector::Material(7,"GLUE$",0.12011E+02,0.60000E+01,0.1930E+01/2.015,999,999);
  // // original
  o2::base::Detector::Material(15, "GLUE$", 12.011, 6, 1.93 / 2.015, 999, 999); // conform with ATLAS, Corrado, Stefan
  o2::base::Detector::Medium(15, "GLUE$", 15, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  o2::base::Detector::Material(16, "ALUMINUM$", 0.26982E+02, 0.13000E+02, 0.26989E+01, 0.89000E+01, 0.99900E+03);
  o2::base::Detector::Medium(16, "ALUMINUM$", 16, 0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);

  o2::base::Detector::Mixture(20, "TUNGCARB$", aWC, zWC, dWC, 2, wWC);
  o2::base::Detector::Medium(20, "TUNGCARB$", 20, 0, ifield, fieldm, tmaxfd, stemax, deemaxSi, epsilSi, stminSi);

  wInox304[3] = 1. - wInox304[0] - wInox304[1] - wInox304[2];
  o2::base::Detector::Mixture(21, "INOX304$", aInox304, zInox304, dInox304, 4, wInox304);
  o2::base::Detector::Medium(21, "INOX304$", 21, 0, ifield, fieldm, tmaxfd, stemax, deemaxSi, epsilSi, stminSi);

  // Tungsten (for gamma converter rods)
  o2::base::Detector::Material(28, "TUNGSTEN$", 183.84, 74, 19.25, 999, 999);
  o2::base::Detector::Medium(28, "TUNGSTEN$", 28, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
}

void Detector::EndOfEvent() { Reset(); }

void Detector::Register()
{
  // This will create a branch in the output tree called Hit, setting the last
  // parameter to kFALSE means that this collection will not be written to the file,
  // it will exist only during the simulation

  if (FairRootManager::Instance()) {
    FairRootManager::Instance()->RegisterAny(addNameTo("Hit").data(), mHits, kTRUE);
  }
}

void Detector::Reset()
{
  if (!o2::utils::ShmManager::Instance().isOperational()) {
    mHits->clear();
  }
}

void Detector::defineWrapperVolume(Int_t id, Double_t rmin, Double_t rmax, Double_t zspan)
{
  // set parameters of id-th wrapper volume
  if (id >= sNumberOfWrapperVolumes || id < 0) {
    LOG(FATAL) << "id " << id << " of wrapper volume is not in 0-" << sNumberOfWrapperVolumes - 1 << " range";
  }

  mWrapperMinRadius[id] = rmin;
  mWrapperMaxRadius[id] = rmax;
  mWrapperZSpan[id] = zspan;
}

void Detector::ConstructGeometry()
{
  // Create the detector materials
  createMaterials();

  // Construct the detector geometry
  constructDetectorGeometry();
}

void Detector::constructDetectorGeometry()
{
  // Create the geometry and insert it in the mother volume EC0V
  TGeoManager* geoManager = gGeoManager;

  TGeoVolume* vALIC = geoManager->GetVolume("barrel");

  if (!vALIC) {
    LOG(FATAL) << "Could not find the top volume";
  }
}

void Detector::addAlignableVolumes() const
{

  return;
}

void Detector::defineSensitiveVolumes()
{
  TGeoManager* geoManager = gGeoManager;
  TGeoVolume* v;

  TString volumeName;
}

Hit* Detector::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                      const TVector3& startMom, double startE, double endTime, double eLoss, unsigned char startStatus,
                      unsigned char endStatus)
{
  mHits->emplace_back(trackID, detID, startPos, endPos, startMom, startE, endTime, eLoss, startStatus, endStatus);
  return &(mHits->back());
}

void Detector::Print(std::ostream* os) const
{
  // Standard output format for this class.
  // Inputs:
  //   ostream *os   The output stream
  // Outputs:
  //   none.
  // Return:
  //   none.

#if defined __GNUC__
#if __GNUC__ > 2
  std::ios::fmtflags fmt;
#else
  Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
  ios::fmtflags fmt;
#else
  Int_t fmt;
#endif
#endif
  // RS: why do we need to pring this garbage?

  // fmt = os->setf(std::ios::scientific); // set scientific floating point output
  // fmt = os->setf(std::ios::hex); // set hex for mStatus only.
  // fmt = os->setf(std::ios::dec); // every thing else decimel.
  //  *os << mModule << " ";
  //  *os << mEnergyDepositionStep << " " << mTof;
  //  *os << " " << mStartingStepX << " " << mStartingStepY << " " << mStartingStepZ;
  //    *os << " " << endl;
  // os->flags(fmt); // reset back to old formating.
  return;
}

void Detector::Read(std::istream* is)
{
  // Standard input format for this class.
  // Inputs:
  //   istream *is  the input stream
  // Outputs:
  //   none.
  // Return:
  //   none.
  // RS no need to read garbage
  return;
}

std::ostream& operator<<(std::ostream& os, Detector& p)
{
  // Standard output streaming function.
  // Inputs:
  //   ostream os  The output stream
  //   Detector p The his to be printed out
  // Outputs:
  //   none.
  // Return:
  //   The input stream

  p.Print(&os);
  return os;
}

std::istream& operator>>(std::istream& is, Detector& r)
{
  // Standard input streaming function.
  // Inputs:
  //   istream is  The input stream
  //   Detector p The Detector class to be filled from this input stream
  // Outputs:
  //   none.
  // Return:
  //   The input stream

  r.Read(&is);
  return is;
}

ClassImp(o2::ec0::Detector);
