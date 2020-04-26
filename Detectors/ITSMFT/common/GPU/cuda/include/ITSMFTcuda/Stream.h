// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Stream.h
/// \brief
///

#ifndef ALICEO2_ITSMFT_GPU_STREAM_H_
#define ALICEO2_ITSMFT_GPU_STREAM_H_

#include "ITSMFTcuda/Definitions.h"

namespace o2
{
namespace itsmft
{
namespace GPU
{

class Stream final
{

 public:
  Stream();
  ~Stream();

  Stream(const Stream&) = delete;
  Stream& operator=(const Stream&) = delete;

  const GPUStream& get() const;

 private:
  GPUStream mStream;
};
} // namespace GPU
} // namespace itsmft
} // namespace o2

#endif /* ALICEO2_ITSMFT_GPU_STREAM_H_ */
