// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef FRAMEWORK_DATAPROCESSING_DEVICE_H
#define FRAMEWORK_DATAPROCESSING_DEVICE_H

#include "Framework/AlgorithmSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ContextRegistry.h"
#include "Framework/DataAllocator.h"
#include "Framework/DataRelayer.h"
#include "Framework/DeviceSpec.h"
#include "Framework/ExpirationHandler.h"
#include "Framework/MessageContext.h"
#include "Framework/RootObjectContext.h"
#include "Framework/ArrowContext.h"
#include "Framework/StringContext.h"
#include "Framework/RawBufferContext.h"
#include "Framework/ServiceRegistry.h"
#include "Framework/InputRoute.h"
#include "Framework/ForwardRoute.h"
#include "Framework/TimingInfo.h"

#include <fairmq/FairMQDevice.h>
#include <fairmq/FairMQParts.h>

#include <memory>

namespace o2
{
namespace framework
{

class DataProcessingDevice : public FairMQDevice
{
 public:
  DataProcessingDevice(const DeviceSpec& spec, ServiceRegistry&);
  void Init() final;
  void PreRun() final;
  void PostRun() final;
  void Reset() final;
  bool ConditionalRun() final;

 protected:
  bool handleData(FairMQParts&);
  bool tryDispatchComputation();
  void error(const char* msg);

 private:
  DeviceSpec const& mSpec;
  AlgorithmSpec::InitCallback mInit;
  AlgorithmSpec::ProcessCallback mStatefulProcess;
  AlgorithmSpec::ProcessCallback mStatelessProcess;
  AlgorithmSpec::ErrorCallback mError;
  std::unique_ptr<ConfigParamRegistry> mConfigRegistry;
  ServiceRegistry& mServiceRegistry;
  TimingInfo mTimingInfo;
  MessageContext mFairMQContext;
  RootObjectContext mRootContext;
  StringContext mStringContext;
  ArrowContext mDataFrameContext;
  RawBufferContext mRawBufferContext;
  ContextRegistry mContextRegistry;
  DataAllocator mAllocator;
  DataRelayer mRelayer;
  std::vector<ExpirationHandler> mExpirationHandlers;

  int mErrorCount;
  int mProcessingCount;
  uint64_t mLastMetricSentTimestamp = 0; /// The timestamp of the last time we sent non-realtime metrics
  uint64_t mBeginIterationTimestamp = 0; /// The timestamp of when the current ConditionalRun was started
};

} // namespace framework
} // namespace o2
#endif
