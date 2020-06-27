// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "EC0Workflow/RecoWorkflow.h"
#include "CommonUtils/ConfigurableParam.h"
#include "EC0tracking/TrackingConfigParam.h"
#include "EC0tracking/Configuration.h"

#include "GPUO2Interface.h"
#include "GPUReconstruction.h"
#include "GPUChainEC0.h"

using namespace o2::framework;

// ------------------------------------------------------------------

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  // option allowing to set parameters
  std::vector<o2::framework::ConfigParamSpec> options{
    {"digits-from-upstream", o2::framework::VariantType::Bool, false, {"digits will be provided from upstream, skip digits reader"}},
    {"clusters-from-upstream", o2::framework::VariantType::Bool, false, {"clusters will be provided from upstream, skip clusterizer"}},
    {"disable-root-output", o2::framework::VariantType::Bool, false, {"do not write output root files"}},
    {"disable-mc", o2::framework::VariantType::Bool, false, {"disable MC propagation even if available"}},
    {"trackerCA", o2::framework::VariantType::Bool, false, {"use trackerCA (default: trackerCM)"}},
    {"entropy-encoding", o2::framework::VariantType::Bool, false, {"produce entropy encoded data"}},
    {"gpuDevice", o2::framework::VariantType::Int, 1, {"use gpu device: CPU=1,CUDA=2,HIP=3 (default: CPU)"}}};

  std::swap(workflowOptions, options);

  std::string keyvaluehelp("Semicolon separated key=value strings (e.g.: 'EC0DigitizerParam.roFrameLength=6000.;...')");
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {keyvaluehelp}});
}

// ------------------------------------------------------------------

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& configcontext)
{
  // Update the (declared) parameters if changed from the command line
  o2::conf::ConfigurableParam::updateFromString(configcontext.options().get<std::string>("configKeyValues"));
  // write the configuration used for the digitizer workflow
  o2::conf::ConfigurableParam::writeINI("o2eclrecoflow_configuration.ini");

  auto useMC = !configcontext.options().get<bool>("disable-mc");
  auto useCAtracker = configcontext.options().get<bool>("trackerCA");
  auto gpuDevice = static_cast<o2::gpu::GPUDataTypes::DeviceType>(configcontext.options().get<int>("gpuDevice"));
  auto extDigits = configcontext.options().get<bool>("digits-from-upstream");
  auto extClusters = configcontext.options().get<bool>("clusters-from-upstream");
  auto disableRootOutput = configcontext.options().get<bool>("disable-root-output");
  auto eencode = configcontext.options().get<bool>("entropy-encoding");
  return std::move(o2::ecl::reco_workflow::getWorkflow(useMC, useCAtracker, gpuDevice, extDigits, extClusters, disableRootOutput, eencode));
}
