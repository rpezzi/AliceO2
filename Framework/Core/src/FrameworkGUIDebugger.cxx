// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/FrameworkGUIDebugger.h"
#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include "Framework/ConfigContext.h"
#include "Framework/ConfigParamRegistry.h"
#include "DebugGUI/imgui.h"
#include "DriverControl.cxx"
#include "DriverInfo.cxx"
#include "FrameworkGUIDeviceInspector.h"
#include "Framework/FrameworkGUIDevicesGraph.h"
#include "Framework/FrameworkGUIDataRelayerUsage.h"
#include "Framework/PaletteHelpers.h"
#include "Framework/FrameworkGUIState.h"

static inline ImVec2 operator+(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x + rhs.x, lhs.y + rhs.y); }
static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }

namespace o2
{
namespace framework
{

static gui::WorkspaceGUIState gState;

ImVec4 colorForLogLevel(LogParsingHelpers::LogLevel logLevel)
{
  switch (logLevel) {
    case LogParsingHelpers::LogLevel::Info:
      return PaletteHelpers::GREEN;
    case LogParsingHelpers::LogLevel::Debug:
      return ImVec4(153. / 255, 61. / 255, 61. / 255, 255. / 255);
    case LogParsingHelpers::LogLevel::Warning:
      return PaletteHelpers::DARK_YELLOW;
    case LogParsingHelpers::LogLevel::Error:
      return PaletteHelpers::RED;
    case LogParsingHelpers::LogLevel::Unknown:
      return ImVec4(153. / 255, 61. / 255, 61. / 255, 255. / 255);
    default:
      return PaletteHelpers::WHITE;
  };
}

bool startsWith(std::string mainStr, std::string toMatch)
{
  if (mainStr.find(toMatch) == 0) {
    return true;
  } else {
    return false;
  }
}

void displayHistory(const DeviceInfo& info, DeviceControl& control)
{
  if (info.history.empty()) {
    return;
  }
  int startPos = info.historyPos;
  const int historySize = info.history.size();

  int triggerStartPos = startPos + 1 % historySize;
  int triggerStopPos = startPos % historySize;

  int j = startPos;
  // We look for a stop trigger, so that we know where to stop the search for
  // out start search. If no stop trigger is found, we search until the end
  if (control.logStopTrigger[0]) {
    while ((j % historySize) != ((startPos + 1) % historySize)) {
      assert(j >= 0);
      assert(j < historySize);
      auto& line = info.history[j];
      if (strstr(line.c_str(), control.logStopTrigger)) {
        triggerStopPos = (j + 1) % historySize;
        break;
      }
      // Wrap in case we end up below 0
      j = (j == 0) ? historySize - 1 : j - 1;
    }
  }

  // Look for the last instance of the start trigger before the
  // last stop trigger.
  j = startPos + 1;
  if (control.logStartTrigger[0]) {
    while ((j % historySize) != triggerStopPos) {
      assert(historySize > j);
      assert(historySize == 1000);
      auto& line = info.history[j];
      if (strstr(line.c_str(), control.logStartTrigger)) {
        triggerStartPos = j;
      }
      j = (j + 1) % historySize;
    }
  }

  // We start from the last trigger found. Eventually this is the first
  // line in the ring buffer, if no trigger is specified.
  size_t ji = triggerStartPos % historySize;
  size_t je = triggerStopPos % historySize;
  size_t iterations = 0;
  while (historySize && ((ji % historySize) != je)) {
    assert(iterations < historySize);
    iterations++;
    assert(historySize == 1000);
    assert(ji < historySize);
    auto& line = info.history[ji];
    auto logLevel = info.historyLevel[ji];

    // Skip empty lines
    if (line.empty()) {
      ji = (ji + 1) % historySize;
      continue;
    }
    // Print matching lines
    if (strstr(line.c_str(), control.logFilter) != 0) {
      auto color = colorForLogLevel(logLevel);
      // We filter twice, once on input, to reduce the
      // stream, a second time at display time, to avoid
      // showing unrelevant messages from past.
      if (logLevel >= control.logLevel) {
        ImGui::TextColored(color, line.c_str(), line.c_str() + line.size());
      }
    }
    ji = (ji + 1) % historySize;
  }
}

template <typename T>
struct HistoData {
  int mod;
  size_t first;
  size_t size;
  const T* points;
};

void historyBar(gui::DeviceGUIState& state, const DeviceSpec& spec, const DeviceMetricsInfo& metricsInfo)
{
  bool open = ImGui::TreeNode(state.label.c_str());
  if (open) {
    ImGui::Text("# channels: %lu", spec.outputChannels.size() + spec.inputChannels.size());
    ImGui::TreePop();
  }
  ImGui::NextColumn();

  if (gState.selectedMetric == -1) {
    ImGui::NextColumn();
    return;
  }

  auto currentMetricName = gState.availableMetrics[gState.selectedMetric];

  size_t i = DeviceMetricsHelper::metricIdxByName(currentMetricName, metricsInfo);
  // We did not find any plot, skipping this.
  if (i == metricsInfo.metricLabelsIdx.size()) {
    ImGui::NextColumn();
    return;
  }
  auto& metric = metricsInfo.metrics[i];

  switch (metric.type) {
    case MetricType::Int: {
      HistoData<int> data;
      data.mod = metricsInfo.timestamps[i].size();
      data.first = metric.pos - data.mod;
      data.size = metricsInfo.intMetrics[metric.storeIdx].size();
      data.points = metricsInfo.intMetrics[metric.storeIdx].data();

      auto getter = [](void* hData, int idx) -> float {
        auto histoData = reinterpret_cast<HistoData<int>*>(hData);
        size_t pos = (histoData->first + static_cast<size_t>(idx)) % histoData->mod;
        assert(pos >= 0 && pos < 1024);
        return histoData->points[pos];
      };
      ImGui::PlotLines(currentMetricName.c_str(), getter, &data, data.size);
      ImGui::NextColumn();
    } break;
    case MetricType::Float: {
      HistoData<float> data;
      data.mod = metricsInfo.timestamps[i].size();
      data.first = metric.pos - data.mod;
      data.size = metricsInfo.floatMetrics[metric.storeIdx].size();
      data.points = metricsInfo.floatMetrics[metric.storeIdx].data();

      auto getter = [](void* hData, int idx) -> float {
        auto histoData = reinterpret_cast<HistoData<float>*>(hData);
        size_t pos = (histoData->first + static_cast<size_t>(idx)) % histoData->mod;
        assert(pos >= 0 && pos < 1024);
        return histoData->points[pos];
      };
      ImGui::PlotLines(currentMetricName.c_str(), getter, &data, data.size);
      ImGui::NextColumn();
    } break;
    default:
      ImGui::NextColumn();
      return;
      break;
  }
}

void displayDeviceHistograms(gui::WorkspaceGUIState &state,
                             const std::vector<DeviceInfo>& infos, const std::vector<DeviceSpec>& devices,
                             std::vector<DeviceControl>& controls, const std::vector<DeviceMetricsInfo>& metricsInfos)
{
  showTopologyNodeGraph(state, infos, devices, controls, metricsInfos);
  if (state.bottomPaneVisible == true) {
    ImGui::SetNextWindowPos(ImVec2(0, ImGui::GetIO().DisplaySize.y - state.bottomPaneSize), 0);
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, state.bottomPaneSize), 0);

    // Calculate the unique set of metrics, as available in the metrics service
    std::set<std::string> allMetricsNames;
    for (const auto& metricsInfo : metricsInfos) {
      for (const auto& labelsPairs : metricsInfo.metricLabelsIdx) {
        allMetricsNames.insert(labelsPairs.first);
      }
    }
    using NamesIndex = std::vector<std::string>;
    gState.availableMetrics.clear();
    std::copy(allMetricsNames.begin(), allMetricsNames.end(), std::back_inserter(gState.availableMetrics));

    ImGui::Begin("Devices", nullptr, ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
    ImGui::Combo("Select metric", &gState.selectedMetric,
                 [](void* data, int idx, const char** outText) -> bool {
                   NamesIndex* v = reinterpret_cast<NamesIndex*>(data);
                   if (idx >= v->size()) {
                     return false;
                   }
                   *outText = v->at(idx).c_str();
                   return true;
                 },
                 &gState.availableMetrics, gState.availableMetrics.size());

    // Calculate the full timestamp range for the selected metric
    if (gState.selectedMetric >= 0) {
      auto currentMetricName = gState.availableMetrics[gState.selectedMetric];
      size_t minTime = -1;
      size_t maxTime = 0;
      for (auto& metricInfo : metricsInfos) {
        size_t mi = DeviceMetricsHelper::metricIdxByName(currentMetricName, metricInfo);
        if (mi == metricInfo.metricLabelsIdx.size()) {
          continue;
        }
        auto& metric = metricInfo.metrics[mi];
        size_t minRangePos = metric.pos % metricInfo.timestamps.size();
        size_t maxRangePos = (size_t)(metric.pos) - 1 % metricInfo.timestamps.size();
        auto& timestamps = metricInfo.timestamps[mi];
        size_t curMinTime = timestamps[minRangePos];
        size_t curMaxTime = timestamps[maxRangePos];
        minTime = minTime < curMinTime ? minTime : curMinTime;
       maxTime = maxTime > curMaxTime ? maxTime : curMaxTime;
      }
      ImGui::Text("min timestamp: %zu, max timestamp: %zu", minTime, maxTime);
    }
    ImGui::BeginChild("ScrollingRegion", ImVec2(0, -ImGui::GetItemsLineHeightWithSpacing()), false,
                      ImGuiWindowFlags_HorizontalScrollbar);
    ImGui::Columns(2);
    ImGui::SetColumnOffset(1, 300);
    for (size_t i = 0; i < gState.devices.size(); ++i) {
      gui::DeviceGUIState& guiState = gState.devices[i];
      const DeviceSpec& spec = devices[i];
      const DeviceMetricsInfo& metricsInfo = metricsInfos[i];

      historyBar(guiState, spec, metricsInfo);
    }
    ImGui::Columns(1);
    ImGui::EndChild();
    ImGui::End();
  }
}

void pushWindowColorDueToStatus(const DeviceInfo& info)
{
  using LogLevel = LogParsingHelpers::LogLevel;
  ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 4.);
  if (info.active == false) {
     ImGui::PushStyleColor(ImGuiCol_TitleBg, PaletteHelpers::DARK_RED);
     ImGui::PushStyleColor(ImGuiCol_TitleBgActive, PaletteHelpers::RED);
     ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, PaletteHelpers::RED);
     return;
  }
  switch (info.maxLogLevel) {
    case LogLevel::Error:
      ImGui::PushStyleColor(ImGuiCol_TitleBg, PaletteHelpers::SHADED_RED);
      ImGui::PushStyleColor(ImGuiCol_TitleBgActive, PaletteHelpers::RED);
      ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, PaletteHelpers::SHADED_RED);
      break;
    case LogLevel::Warning:
      ImGui::PushStyleColor(ImGuiCol_TitleBg, PaletteHelpers::SHADED_YELLOW);
      ImGui::PushStyleColor(ImGuiCol_TitleBgActive, PaletteHelpers::YELLOW);
      ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, PaletteHelpers::SHADED_YELLOW);
      break;
    case LogLevel::Info:
      ImGui::PushStyleColor(ImGuiCol_TitleBg, PaletteHelpers::SHADED_GREEN);
      ImGui::PushStyleColor(ImGuiCol_TitleBgActive, PaletteHelpers::GREEN);
      ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, PaletteHelpers::SHADED_GREEN);
      break;
    default:
      ImGui::PushStyleColor(ImGuiCol_TitleBg, PaletteHelpers::SHADED_BLUE);
      ImGui::PushStyleColor(ImGuiCol_TitleBgActive, PaletteHelpers::BLUE);
      ImGui::PushStyleColor(ImGuiCol_TitleBgCollapsed, PaletteHelpers::SHADED_BLUE);
      break;
  }
}

void popWindowColorDueToStatus()
{
  ImGui::PopStyleColor(3);
  ImGui::PopStyleVar(1);
}

struct DriverHelper {
  static char const* stateToString(enum DriverState state)
  {
    static const char* names[static_cast<int>(DriverState::LAST)] = {
      "INIT",             //
      "SCHEDULE",         //
      "RUNNING",          //
      "GUI",              //
      "REDEPLOY_GUI",     //
      "QUIT_REQUESTED",   //
      "HANDLE_CHILDREN",  //
      "EXIT",             //
      "UNKNOWN"           //
      "PERFORM_CALLBACKS" //
    };
    return names[static_cast<int>(state)];
  }
};

/// Display information window about the driver
/// and its state.
void displayDriverInfo(DriverInfo const& driverInfo, DriverControl& driverControl)
{
  ImGui::Begin("Driver information");
  ImGui::Text("Numer of running devices: %lu", driverInfo.socket2DeviceInfo.size() / 2);

  if (driverControl.state == DriverControlState::STEP) {
    driverControl.state = DriverControlState::PAUSE;
  }
  auto state = reinterpret_cast<int*>(&driverControl.state);
  ImGui::RadioButton("Play", state, static_cast<int>(DriverControlState::PLAY));
  ImGui::SameLine();
  ImGui::RadioButton("Pause", state, static_cast<int>(DriverControlState::PAUSE));
  ImGui::SameLine();
  ImGui::RadioButton("Step", state, static_cast<int>(DriverControlState::STEP));

  if (driverControl.state == DriverControlState::PAUSE) {
    driverControl.forcedTransitions.push_back(DriverState::GUI);
  }

  auto &registry = driverInfo.configContext->options();
  ImGui::TextUnformatted("Workflow options:");
  ImGui::Columns(2);
  for (auto &option : driverInfo.workflowOptions) {
    ImGui::TextUnformatted(option.name.c_str());
    ImGui::NextColumn();
    switch (option.type) {
      case ConfigParamSpec::ParamType::Int64:
      case ConfigParamSpec::ParamType::Int:
        ImGui::Text("%d", registry.get<int>(option.name.c_str()));
        break;
      case ConfigParamSpec::ParamType::Float:
        ImGui::Text("%f", registry.get<float>(option.name.c_str()));
        break;
      case ConfigParamSpec::ParamType::Double:
        ImGui::Text("%f", registry.get<double>(option.name.c_str()));
        break;
      case ConfigParamSpec::ParamType::String:
        ImGui::Text("%s", registry.get<std::string>(option.name.c_str()).c_str());
        break;
      case ConfigParamSpec::ParamType::Bool:
        ImGui::TextUnformatted(registry.get<bool>(option.name.c_str()) ? "true" : "false");
        break;
      case ConfigParamSpec::ParamType::Empty:
      case ConfigParamSpec::ParamType::Unknown:
        break;
    }
    ImGui::NextColumn();
  }
  ImGui::Columns();

  ImGui::Text("State stack (depth %lu)", driverInfo.states.size());

  for (size_t i = 0; i < driverInfo.states.size(); ++i) {
    ImGui::Text("#%lu: %s", i, DriverHelper::stateToString(driverInfo.states[i]));
  }

  ImGui::End();
}

// FIXME: return empty function in case we were not built
// with GLFW support.
///
std::function<void(void)> getGUIDebugger(const std::vector<DeviceInfo>& infos, const std::vector<DeviceSpec>& devices,
                                         const std::vector<DeviceMetricsInfo>& metricsInfos,
                                         const DriverInfo& driverInfo, std::vector<DeviceControl>& controls,
                                         DriverControl& driverControl)
{
  gState.selectedMetric = -1;
  // FIXME: this should probaly have a better mapping between our window state and
  gState.devices.resize(infos.size());
  for (size_t i = 0; i < gState.devices.size(); ++i) {
    gui::DeviceGUIState& state = gState.devices[i];
    state.label = devices[i].id + "(" + std::to_string(infos[i].pid) + ")";
  }
  gState.bottomPaneSize = 300;
  gState.leftPaneSize = 100;
  gState.rightPaneSize = 300;

  // Show all the panes by default.
  gState.bottomPaneVisible = true;
  gState.leftPaneVisible = true;
  gState.rightPaneVisible = true;

  return [&infos, &devices, &controls, &metricsInfos, &driverInfo, &driverControl]() {
    ImGuiStyle& style = ImGui::GetStyle();
    style.FrameRounding = 0.;
    style.WindowRounding = 0.;
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0x1b / 255.f, 0x1b / 255.f, 0x1b / 255.f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0x1b / 255.f, 0x1b / 255.f, 0x1b / 255.f, 1.00f);

    displayDeviceHistograms(gState, infos, devices, controls, metricsInfos);
    displayDriverInfo(driverInfo, driverControl);

    int windowPosStepping = (ImGui::GetIO().DisplaySize.y - 500) / gState.devices.size();

    for (size_t i = 0; i < gState.devices.size(); ++i) {
      gui::DeviceGUIState& state = gState.devices[i];
      assert(i < infos.size());
      assert(i < devices.size());
      const DeviceInfo& info = infos[i];
      const DeviceSpec& spec = devices[i];
      const DeviceMetricsInfo& metrics = metricsInfos[i];

      assert(controls.size() == devices.size());
      DeviceControl& control = controls[i];

      pushWindowColorDueToStatus(info);
      if (control.logVisible) {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x / 3 * 2, i * windowPosStepping), ImGuiSetCond_Once);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x / 3, ImGui::GetIO().DisplaySize.y - 300),
                                 ImGuiSetCond_Once);
        ImGui::Begin(state.label.c_str(), &control.logVisible);

        ImGui::InputText("Log filter", control.logFilter, sizeof(control.logFilter));
        ImGui::InputText("Log start trigger", control.logStartTrigger, sizeof(control.logStartTrigger));
        ImGui::InputText("Log stop trigger", control.logStopTrigger, sizeof(control.logStopTrigger));
        ImGui::Checkbox("Stop logging", &control.quiet);
        ImGui::SameLine();
        ImGui::Combo("Log level", reinterpret_cast<int*>(&control.logLevel), LogParsingHelpers::LOG_LEVELS,
                     (int)LogParsingHelpers::LogLevel::Size, 5);

        ImGui::Separator();
        ImGui::BeginChild("ScrollingRegion", ImVec2(0, -ImGui::GetItemsLineHeightWithSpacing()), false,
                          ImGuiWindowFlags_HorizontalScrollbar|ImGuiWindowFlags_NoMove);
        displayHistory(info, control);
        ImGui::EndChild();
        ImGui::End();
      }
      popWindowColorDueToStatus();
    }
  };
}

} // namespace framework
} // namespace o2
