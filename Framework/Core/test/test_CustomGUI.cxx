// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include <iostream>
#include <boost/algorithm/string.hpp>

#include "Framework/InputSpec.h"
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/ParallelContext.h"
#include "Framework/runDataProcessing.h"
#include "Framework/DebugGUI.h"
#include "DebugGUI/imgui.h"

using namespace o2::framework;
using DataHeader = o2::header::DataHeader;

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    { "source",
      Inputs{},
      {
        OutputSpec{ { "test" }, "TST", "A" },
      },
      AlgorithmSpec{
        [](ProcessingContext& ctx) {
          sleep(1);
          auto out = ctx.outputs().make<int>(OutputRef{ "test", 0 });
        } } },
    { "dest",
      Inputs{
        { "test", "TST", "A" } },
      Outputs{},
      AlgorithmSpec{
        [](InitContext& ic) {
          void* window = initGUI("A test window");
          auto count = std::make_shared<int>(0);
          auto guiCallback = [count]() {
            ImGui::Begin("Some sub-window");
            ImGui::Text("Counter value: %i", *count);
            ImGui::End();
          };
          ic.services().get<CallbackService>().set(CallbackService::Id::ClockTick,
                                                   [count, window, guiCallback]() {
                    (*count)++; window ? pollGUI(window, guiCallback) : false; });
          return [count](ProcessingContext& ctx) {
            if (*count > 1000) {
              ctx.services().get<ControlService>().readyToQuit(true);
            }
          };
        } } }
  };
}
