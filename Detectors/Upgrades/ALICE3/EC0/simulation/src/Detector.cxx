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
#include "EC0Simulation/EC0Layer.h"

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
    mHits(o2::utils::createSimVector<Hit>())
{
}

static void configEC0(Detector* ec0)
{
  // build EC0 upgrade detector
}

Detector::Detector(Bool_t active)
  : o2::base::DetImpl<Detector>("EC0", active),
    mTrackData(),
    mHits(o2::utils::createSimVector<Hit>())
{

  configEC0(this);
}

Detector::Detector(const Detector& rhs)
  : o2::base::DetImpl<Detector>(rhs),
    mTrackData(),

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
