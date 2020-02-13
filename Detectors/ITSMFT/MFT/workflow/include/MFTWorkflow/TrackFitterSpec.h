// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackFitterSpec.h
/// \brief Definition of a data processor to read, refit and send tracks with attached clusters
///
/// \author Philippe Pillot, Subatech

#ifndef ALICEO2_MFT_TRACKFITTERSPEC_H_
#define ALICEO2_MFT_TRACKFITTERSPEC_H_

#include "MFTTracking/TrackFitter.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"


namespace o2
{
namespace mft
{

class TrackFitterTask : public o2::framework::Task
{
 public:
  TrackFitterTask() {}
  ~TrackFitterTask() override = default;
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;

 private:
  int mState = 0;
  std::unique_ptr<o2::mft::TrackFitter> mTrackFitter = nullptr;
};


o2::framework::DataProcessorSpec getTrackFitterSpec();

} // end namespace mft
} // end namespace o2

#endif // ALICEO2_MFT_TRACKFITTERSPEC_H_
