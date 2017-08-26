// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Digitizer.cxx
/// \brief Implementation of the ITS digitizer

#include "ITSMFTBase/Digit.h"
#include "ITSMFTBase/SegmentationPixel.h"
#include "ITSMFTSimulation/Point.h"
#include "ITSSimulation/Digitizer.h"

#include "FairLogger.h"   // for LOG
#include "TClonesArray.h" // for TClonesArray
#include <TRandom.h>
#include <climits>

ClassImp(o2::ITS::Digitizer)

using o2::ITSMFT::Point;
using o2::ITSMFT::Chip;
using o2::ITSMFT::SimulationAlpide;
using o2::ITSMFT::SegmentationPixel;
using o2::ITSMFT::Digit;

using namespace o2::ITS;

//_______________________________________________________________________
Digitizer::Digitizer() : mGeometry(), mSimulations() {}


//_______________________________________________________________________
Digitizer::~Digitizer() = default;

//_______________________________________________________________________
void Digitizer::init(Bool_t build)
{
  mGeometry.Build(build);
  const Int_t numOfChips = mGeometry.getNumberOfChips();
  
  for (Int_t i = 0; i < numOfChips; i++) {
    mSimulations.emplace_back(&mParams, i, mGeometry.getMatrixSensor(i));
  }
}

//_______________________________________________________________________
void Digitizer::process(TClonesArray* points, TClonesArray* digits)
{
  // digitize single event

  
  const Int_t numOfChips = mGeometry.getNumberOfChips();  
  const SegmentationPixel* seg = (SegmentationPixel*)mGeometry.getSegmentationById(0);
  
  // estimate the smalles RO Frame this event may have
  double hTime0 = mEventTime - mParams.getTimeOffset();
  if (hTime0 > UINT_MAX) {
    LOG(WARNING) << "min Hit RO Frame undefined: time: " << hTime0 << " is in far future: "
		 << " EventTime: " << mEventTime << " TimeOffset: "
		 << mParams.getTimeOffset() << FairLogger::endl;
    return;
  }
  
  UInt_t minNewROFrame = static_cast<UInt_t>(hTime0)/mParams.getROFrameLenght();

  LOG(INFO) << "Digitizing ITS event at time " << mEventTime
	    << " (TOffset= " << mParams.getTimeOffset() << " ROFrame= " << minNewROFrame << ")"
	    << " cont.mode: " << isContinuous() << " current Min/Max RO Frames "
	    << mROFrameMin << "/" << mROFrameMax << FairLogger::endl ;
  
  if (mParams.isContinuous() && minNewROFrame>mROFrameMin) {
    // if there are already digits cached for previous RO Frames AND the new event
    // cannot contribute to these digits, move them to the output container
    if (mROFrameMax<minNewROFrame) mROFrameMax = minNewROFrame-1;
    for (auto rof=mROFrameMin; rof<minNewROFrame; rof++) {
      fillOutputContainer(digits, rof);
    }
    //    fillOutputContainer(digits, minNewROFrame-1);
  }
  
  // accumulate points for every chip
  TIter nextPoint(points);
  Point* point = nullptr;
  while ( (point = (Point*)nextPoint()) ) {
    mSimulations[point->GetDetectorID()].InsertPoint(point);
  }
    
  // Convert points to digits  
  for (auto &simulation : mSimulations) {
    simulation.Points2Digits(seg, mEventTime, mROFrameMin, mROFrameMax);
    simulation.ClearPoints();
  }

  // in the triggered mode store digits after every MC event
  if (!mParams.isContinuous()) {
    fillOutputContainer(digits, mROFrameMax);
  }
}

//_______________________________________________________________________
void Digitizer::setEventTime(double t)
{
  // assign event time, it should be in a strictly increasing order
  // convert to ns
  t *= mCoeffToNanoSecond;
  
  if (t<mEventTime && mParams.isContinuous()) {
    LOG(FATAL) << "New event time (" << t << ") is < previous event time (" << mEventTime << ")" << FairLogger::endl;
  }
  mEventTime = t;
  // to limit the range of RO Frame IDs we subtract the meaningless offset
  if (mParams.isContinuous()) { // in continuous mode we set the offset only in the very beginning
    if (!mParams.isTimeOffsetSet()) { // offset is initially at -inf
      mParams.setTimeOffset(mEventTime + mParams.getROFrameLenght()*(gRandom->Rndm()-0.5));
    }
  }
  else { // in the triggered mode we start from 0 ROFrame in every event
    mParams.setTimeOffset(+0.);
    mROFrameMin = 0;  // so we reset the frame counters
    mROFrameMax = 0;
  }
}

//_______________________________________________________________________
void Digitizer::fillOutputContainer(TClonesArray* digits, UInt_t maxFrame)
{
  // fill output with digits ready to be stored, generating the noise beforehand
  if (maxFrame>mROFrameMax) maxFrame = mROFrameMax;
  const SegmentationPixel* seg = (SegmentationPixel*)mGeometry.getSegmentationById(0);

  LOG(INFO) << "Filling ITS digits output for RO frames " << mROFrameMin << ":" << maxFrame << FairLogger::endl ;

  for (auto &simulation : mSimulations) {
    // add the random noise to all ROFrame being stored
    simulation.addNoise(seg,mROFrameMin,maxFrame);
  }

  // we have to write chips in RO increasing order, therefore have to loop over the frames here
  for (auto rof=mROFrameMin;rof<=maxFrame;rof++) {
    for (auto &simulation : mSimulations) {
      simulation.fillOutputContainer(digits,rof);
    }
  }
  mROFrameMin = maxFrame+1;
}
