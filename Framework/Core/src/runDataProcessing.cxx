// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "FairMQDevice.h"
#include "Framework/ChannelConfigurationPolicy.h"
#include "Framework/ChannelMatching.h"
#include "Framework/ConfigParamsHelper.h"
#include "Framework/DataProcessingDevice.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSourceDevice.h"
#include "Framework/DebugGUI.h"
#include "Framework/DeviceControl.h"
#include "Framework/DeviceExecution.h"
#include "Framework/DeviceInfo.h"
#include "Framework/DeviceMetricsInfo.h"
#include "Framework/DeviceSpec.h"
#include "Framework/FrameworkGUIDebugger.h"
#include "Framework/LocalRootFileService.h"
#include "Framework/LogParsingHelpers.h"
#include "Framework/ParallelContext.h"
#include "Framework/RawDeviceService.h"
#include "Framework/SimpleMetricsService.h"
#include "Framework/SimpleRawDeviceService.h"
#include "Framework/TextControlService.h"
#include "Framework/WorkflowSpec.h"

#include "DDSConfigHelpers.h"
#include "DeviceSpecHelpers.h"
#include "DriverInfo.h"
#include "DriverControl.h"
#include "GraphvizHelpers.h"
#include "options/FairMQProgOptions.h"

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <type_traits>

#include <sys/resource.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <csignal>

#include <fairmq/DeviceRunner.h>
#include <fairmq/FairMQLogger.h>

using namespace o2::framework;

std::vector<DeviceInfo> gDeviceInfos;
std::vector<DeviceMetricsInfo> gDeviceMetricsInfos;
std::vector<DeviceControl> gDeviceControls;
std::vector<DeviceExecution> gDeviceExecutions;

namespace bpo = boost::program_options;

// FIXME: probably find a better place
// these are the device options added by the framework, but they can be
// overloaded in the config spec
bpo::options_description gHiddenDeviceOptions("Hidden child options");

// Read from a given fd and print it.
// return true if we can still read from it,
// return false if we need to close the input pipe.
//
// FIXME: We should really print full lines.
bool getChildData(int infd, DeviceInfo& outinfo)
{
  char buffer[1024];
  int bytes_read;
  // NOTE: do not quite understand read ends up blocking if I read more than
  //        once. Oh well... Good enough for now.
  // do {
  bytes_read = read(infd, buffer, 1024);
  if (bytes_read == 0) {
    return false;
  }
  if (bytes_read == -1) {
    switch (errno) {
      case EAGAIN:
        return true;
      default:
        return false;
    }
  }
  assert(bytes_read > 0);
  outinfo.unprinted += std::string(buffer, bytes_read);
  //  } while (bytes_read != 0);
  return true;
}

/// Return true if all the DeviceInfo in \a infos are
/// ready to quit. false otherwise.
/// FIXME: move to an helper class
bool checkIfCanExit(std::vector<DeviceInfo> const& infos) {
  if (infos.empty()) {
    return false;
  }
  for (auto& info : infos) {
    if (info.readyToQuit == false) {
       return false;
    }
  }
  return true;
}

// Kill all the active children. Exit code
// is != 0 if any of the children had an error.
int killChildren(std::vector<DeviceInfo>& infos)
{
  int exitCode = 0;
  for (auto& info : infos) {
    if (exitCode == 0 && info.maxLogLevel >= LogParsingHelpers::LogLevel::Error) {
      LOG(ERROR) << "Child " << info.pid << " had at least one "
                 << "message above severity ERROR: " << info.lastError;
      exitCode = 1;
    }
    if (!info.active) {
      continue;
    }
    kill(info.pid, SIGKILL);
    int status;
    waitpid(info.pid, &status, 0);
  }
  return exitCode;
}

int createPipes(int maxFd, int* pipes)
{
  auto p = pipe(pipes);
  maxFd = maxFd > pipes[0] ? maxFd : pipes[0];
  maxFd = maxFd > pipes[1] ? maxFd : pipes[1];

  if (p == -1) {
    std::cerr << "Unable to create PIPE: ";
    switch (errno) {
      case EFAULT:
        assert(false && "EFAULT while reading from pipe");
        break;
      case EMFILE:
        std::cerr << "Too many active descriptors";
        break;
      case ENFILE:
        std::cerr << "System file table is full";
        break;
      default:
        std::cerr << "Unknown PIPE" << std::endl;
    };
    // Kill immediately both the parent and all the children
    kill(-1 * getpid(), SIGKILL);
  }
  return maxFd;
}

// FIXME: I should really do this gracefully, by doing the following:
// - Kill all the children
// - Set a sig_atomic_t to say we did.
// - Wait for all the children to exit
// - Return gracefully.
static void handle_sigint(int signum)
{
  auto exitCode = killChildren(gDeviceInfos);
  // We kill ourself after having killed all our children (SPOOKY!)
  signal(SIGINT, SIG_DFL);
  kill(getpid(), SIGINT);
}

/// This will start a new device by forking and executing a
/// new child
void spawnDevice(DeviceSpec const& spec,
                 std::map<int, size_t> &socket2DeviceInfo,
                 DeviceControl &control,
                 DeviceExecution &execution,
                 std::vector<DeviceInfo> &deviceInfos,
                 int &maxFd, fd_set &childFdset) {
  int childstdout[2];
  int childstderr[2];

  maxFd = createPipes(maxFd, childstdout);
  maxFd = createPipes(maxFd, childstderr);

  // If we have a framework id, it means we have already been respawned
  // and that we are in a child. If not, we need to fork and re-exec, adding
  // the framework-id as one of the options.
  pid_t id = 0;
  id = fork();
  // We are the child: prepare options and reexec.
  if (id == 0) {
    // We allow being debugged and do not terminate on SIGTRAP
    signal(SIGTRAP, SIG_IGN);

    // We do not start the process if control.noStart is set.
    if (control.stopped) {
      kill(getpid(), SIGSTOP);
    }

    // This is the child. We close the read part of the pipe, stdout
    // and dup2 the write part of the pipe on it. Then we can restart.
    close(childstdout[0]);
    close(childstderr[0]);
    close(STDOUT_FILENO);
    close(STDERR_FILENO);
    dup2(childstdout[1], STDOUT_FILENO);
    dup2(childstderr[1], STDERR_FILENO);
    execvp(execution.args[0], execution.args.data());
  }

  // This is the parent. We close the write end of
  // the child pipe and and keep track of the fd so
  // that we can later select on it.
  struct sigaction sa_handle_int;
  sa_handle_int.sa_handler = handle_sigint;
  sigemptyset(&sa_handle_int.sa_mask);
  sa_handle_int.sa_flags = SA_RESTART;
  if (sigaction(SIGINT, &sa_handle_int, nullptr) == -1) {
    perror("Unable to install signal handler");
    exit(1);
  }

  std::cout << "Starting " << spec.id << " on pid " << id << "\n";
  DeviceInfo info;
  info.pid = id;
  info.active = true;
  info.readyToQuit = false;
  info.historySize = 1000;
  info.historyPos = 0;
  info.maxLogLevel = LogParsingHelpers::LogLevel::Debug;

  socket2DeviceInfo.insert(std::make_pair(childstdout[0], deviceInfos.size()));
  socket2DeviceInfo.insert(std::make_pair(childstderr[0], deviceInfos.size()));
  deviceInfos.emplace_back(info);
  // Let's add also metrics information for the given device
  gDeviceMetricsInfos.emplace_back(DeviceMetricsInfo{});

  close(childstdout[1]);
  close(childstderr[1]);
  FD_SET(childstdout[0], &childFdset);
  FD_SET(childstderr[0], &childFdset);
}

void doRunning(DriverInfo &driverInfo,
               std::vector<DeviceInfo> &infos,
               std::vector<DeviceSpec> const &specs,
               std::vector<DeviceControl> &controls,
               std::vector<DeviceMetricsInfo> &metricsInfos) {
  // Wait for children to say something. When they do
  // print it.
  fd_set fdset;
  timeval timeout;
  timeout.tv_sec = 0;
  timeout.tv_usec = 16666; // This should be enough to allow 60 HZ redrawing.
  memcpy(&fdset, &driverInfo.childFdset, sizeof(fd_set));
  int numFd = select(driverInfo.maxFd, &fdset, nullptr, nullptr, &timeout);
  if (numFd == 0) {
    return;
  }
  for (int si = 0; si < driverInfo.maxFd; ++si) {
    if (FD_ISSET(si, &fdset)) {
      assert(driverInfo.socket2DeviceInfo.find(si) != driverInfo.socket2DeviceInfo.end());
      auto& info = infos[driverInfo.socket2DeviceInfo[si]];

      bool fdActive = getChildData(si, info);
      // If the pipe was closed due to the process exiting, we
      // can avoid the select.
      if (!fdActive) {
        info.active = false;
        close(si);
        FD_CLR(si, &driverInfo.childFdset);
      }
      --numFd;
    }
    // FIXME: no need to check after numFd gets to 0.
  }
  // Display part. All you need to display should actually be in
  // `infos`.
  // TODO: split at \n
  // TODO: update this only once per 1/60 of a second or
  // things like this.
  // TODO: have multiple display modes
  // TODO: graphical view of the processing?
  assert(infos.size() == controls.size());
  std::smatch match;
  std::string token;
  const std::string delimiter("\n");
  for (size_t di = 0, de = infos.size(); di < de; ++di) {
    DeviceInfo& info = infos[di];
    DeviceControl& control = controls[di];
    DeviceMetricsInfo& metrics = metricsInfos[di];

    if (info.unprinted.empty()) {
      continue;
    }

    auto s = info.unprinted;
    size_t pos = 0;
    info.history.resize(info.historySize);
    info.historyLevel.resize(info.historySize);

    while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      auto logLevel = LogParsingHelpers::parseTokenLevel(token);

      // Check if the token is a metric from SimpleMetricsService
      // if yes, we do not print it out and simply store it to be displayed
      // in the GUI.
      // Then we check if it is part of our Poor man control system
      // if yes, we execute the associated command.
      if (parseMetric(token, match)) {
        LOG(DEBUG) << "Found metric with key " << match[2] << " and value " << match[4];
        processMetric(match, metrics);
      } else if (parseControl(token, match)) {
        auto command = match[1];
        auto validFor = match[2];
        LOG(DEBUG) << "Found control command " << command << " from pid " << info.pid << " valid for " << validFor;
        if (command == "QUIT") {
          if (validFor == "ALL") {
            for (auto& deviceInfo : infos) {
              deviceInfo.readyToQuit = true;
            }
          }
        }
      } else if (!control.quiet && (strstr(token.c_str(), control.logFilter) != nullptr) &&
                 logLevel >= control.logLevel) {
        assert(info.historyPos >= 0);
        assert(info.historyPos < info.history.size());
        info.history[info.historyPos] = token;
        info.historyLevel[info.historyPos] = logLevel;
        info.historyPos = (info.historyPos + 1) % info.history.size();
        std::cout << "[" << info.pid << "]: " << token << std::endl;
      }
      // We keep track of the maximum log error a
      // device has seen.
      if (logLevel > info.maxLogLevel && logLevel > LogParsingHelpers::LogLevel::Info
          && logLevel != LogParsingHelpers::LogLevel::Unknown) {
        info.maxLogLevel = logLevel;
      }
      if (logLevel == LogParsingHelpers::LogLevel::Error) {
        info.lastError = token;
      }
      s.erase(0, pos + delimiter.length());
    }
    info.unprinted = s;
  }
  // FIXME: for the gui to work correctly I would actually need to
  //        run the loop more often and update whenever enough time has
  //        passed.
}

void handle_sigchld(int sig)
{
  int saved_errno = errno;
  std::vector<pid_t> pids;
  while (true) {
    pid_t pid = waitpid((pid_t)(-1), nullptr, WNOHANG);
    if (pid > 0) {
      pids.push_back(pid);
      continue;
    } else {
      break;
    }
  }
  errno = saved_errno;
  for (auto& pid : pids) {
    printf("Child exited: %d\n", pid);
    for (auto& info : gDeviceInfos) {
      if (info.pid == pid) {
        info.active = false;
      }
    }
    fflush(stdout);
  }
}
// This is the handler for the parent inner loop.
// So far the only responsibility for it are:
//
// - Echo children output in a sensible manner
//
//
// - TODO: allow single child view?
// - TODO: allow last line per child mode?
// - TODO: allow last error per child mode?
void doParent(std::vector<DeviceInfo> &infos,
              std::vector<DeviceSpec> &deviceSpecs,
              std::vector<DeviceControl> &controls,
              std::vector<DeviceMetricsInfo> &metricsInfos,
              std::vector<DeviceExecution> &deviceExecutions,
              std::vector<DriverState> const &startProgram,
              bool batch)
{
  DriverInfo driverInfo;
  driverInfo.maxFd = 0;
  driverInfo.states.reserve(10);
  driverInfo.states = startProgram;
  DriverControl driverControl;

  void* window = nullptr;
  decltype(getGUIDebugger(infos, deviceSpecs, metricsInfos, driverInfo, controls, driverControl)) debugGUICallback;

  if (batch == false) {
    window = initGUI("O2 Framework debug GUI");
  }
  if (batch == false && window == nullptr) {
    LOG(WARN) << "Could not create GUI. Switching to batch mode. Do you have GLFW on your system?";
    batch = true;
  }
  bool guiQuitRequested = false;
  // FIXME: I should really have some way of exiting the
  // parent..
  DriverState current;
  while (true) {
    if (driverInfo.states.empty() == false) {
      current = driverInfo.states.back();
    } else {
      current = DriverState::UNKNOWN;
    }
    driverInfo.states.pop_back();
    switch(current) {
      case DriverState::INIT:
        LOG(INFO) << "Initialising O2 Data Processing Layer";

        // Install signal handler for quitting children.
        driverInfo.sa_handle_child.sa_handler = &handle_sigchld;
        sigemptyset(&driverInfo.sa_handle_child.sa_mask);
        driverInfo.sa_handle_child.sa_flags = SA_RESTART | SA_NOCLDSTOP;
        if (sigaction(SIGCHLD, &driverInfo.sa_handle_child, nullptr) == -1) {
          perror(nullptr);
          exit(1);
        }
        FD_ZERO(&(driverInfo.childFdset));

        /// After INIT we go into RUNNING and eventually to SCHEDULE from
        /// there and back into running. This is because the general case
        /// would be that we start an application and then we wait for
        /// resource offers from DDS or whatever resource manager we use.
        driverInfo.states.push_back(DriverState::RUNNING);
        driverInfo.states.push_back(DriverState::GUI);
//        driverInfo.states.push_back(DriverState::REDEPLOY_GUI);
        LOG(INFO) << "O2 Data Processing Layer initialised. We brake for nobody.";
        break;
      case DriverState::REDEPLOY_GUI:
        // The callback for the GUI needs to be recalculated every time
        // the deployed configuration changes, e.g. a new device 
        // has been added to the topology.
        // We need to recreate the GUI callback every time we reschedule
        // because getGUIDebugger actually recreates the GUI state.
        if (window) {
          debugGUICallback = getGUIDebugger(infos,
              deviceSpecs,
              metricsInfos,
              driverInfo,
              controls,
              driverControl);
        }
        break;
      case DriverState::SCHEDULE:
        LOG(INFO) << "Redeployment of configuration asked.";
        for (size_t di = 0; di < deviceSpecs.size(); ++di) {
          spawnDevice(deviceSpecs[di],
                      driverInfo.socket2DeviceInfo,
                      controls[di],
                      deviceExecutions[di],
                      infos,
                      driverInfo.maxFd, driverInfo.childFdset);
        }
        driverInfo.maxFd += 1;
        assert(infos.empty() == false);
        // Go to RUNNING next.
        driverInfo.states.push_back(DriverState::RUNNING);
        driverInfo.states.push_back(DriverState::GUI);
        driverInfo.states.push_back(DriverState::REDEPLOY_GUI);
        LOG(INFO) << "Redeployment of configuration done.";
        break;
      case DriverState::RUNNING:
        // Calculate what we should do next and eventually
        // show the GUI
        if (batch
            || guiQuitRequested
            || (checkIfCanExit(infos) == true)) {
          // Something requested to quit. Let's update the GUI
          // one more time and then EXIT.
          LOG(INFO) << "Quitting";
          driverInfo.states.push_back(DriverState::EXIT);
          driverInfo.states.push_back(DriverState::GUI);
        } else if (infos.size() != deviceSpecs.size()) {
          // If the number of deviceSpecs is different from
          // the DeviceInfos it means the speicification
          // does not match what is running, so we need to do
          // further scheduling.
          driverInfo.states.push_back(DriverState::SCHEDULE);
          driverInfo.states.push_back(DriverState::GUI);
        }
        else if (deviceSpecs.size() == 0) {
          LOG(INFO) << "No device resulting from the workflow. Quitting.";
          // If there are no deviceSpecs, we exit. 
          driverInfo.states.push_back(DriverState::EXIT);
        }
        else {
          driverInfo.states.push_back(DriverState::RUNNING);
          driverInfo.states.push_back(DriverState::GUI);
        }
        // Update information about devices
        doRunning(driverInfo, infos, deviceSpecs, controls, metricsInfos);
        break;
      case DriverState::GUI:
        guiQuitRequested = (pollGUI(window, debugGUICallback) == false);
        break;
      case DriverState::EXIT:
        killChildren(infos);
        return;
      default:
        LOG(ERROR) << "Driver transitioned in an unknown state. Shutting down.";
        killChildren(infos);
        return;
    }
  }
}


// Magic to support both old and new FairRoot logger behavior
// is_define<fair::Logger> only works if Logger is fully defined,
// not forward declared.
namespace fair
{
// This is required to be declared (but not defined) for the SFINAE below
// to work
class Logger;
namespace mq
{
namespace logger
{
// This we can leave it empty, as the sole purpose is to
// allow the using namespace below.
}
} // namespace mq
} // namespace fair

template <class, class, class = void>
struct is_defined : std::false_type {
};

template <class S, class T>
struct is_defined<S, T, std::enable_if_t<std::is_object<T>::value && !std::is_pointer<T>::value && (sizeof(T) > 0)>>
  : std::true_type {
  using logger = T;
};

using namespace fair::mq::logger;

struct TurnOffColors {

  template <typename T>
  static auto apply(T value, int) -> decltype(ReinitLogger(value), void())
  {
    ReinitLogger(value);
  }
  template <typename T>
  static auto apply(T value, long) -> typename std::enable_if<is_defined<T, fair::Logger>::value == true, void>::type
  {
    // By using is_defined<T, fair::Logger>::logger rather than fair::Logger
    // directly, we avoid the error about using an incomplete type in a nested
    // expression.
    is_defined<T, fair::Logger>::logger::SetConsoleColor(value);
  }
};

int doChild(int argc, char** argv, const o2::framework::DeviceSpec& spec)
{
  TurnOffColors::apply(false, 0);
  LOG(INFO) << "Spawing new device " << spec.id << " in process with pid " << getpid();

  try {
    fair::mq::DeviceRunner runner{ argc, argv };

    // Populate options from the command line. Notice that only the options
    // declared in the workflow definition are allowed.
    runner.AddHook<fair::mq::hooks::SetCustomCmdLineOptions>([&spec](fair::mq::DeviceRunner& r) {
      boost::program_options::options_description optsDesc;
      populateBoostProgramOptions(optsDesc, spec.options, gHiddenDeviceOptions);
      r.fConfig.AddToCmdLineOptions(optsDesc, true);
    });

    // We initialise this in the driver, because different drivers might have
    // different versions of the service
    ServiceRegistry serviceRegistry;
    serviceRegistry.registerService<MetricsService>(new SimpleMetricsService());
    serviceRegistry.registerService<RootFileService>(new LocalRootFileService());
    serviceRegistry.registerService<ControlService>(new TextControlService());
    serviceRegistry.registerService<ParallelContext>(new ParallelContext(spec.rank, spec.nSlots));

    std::unique_ptr<FairMQDevice> device;
    serviceRegistry.registerService<RawDeviceService>(new SimpleRawDeviceService(nullptr));

    if (spec.inputs.empty()) {
      LOG(DEBUG) << spec.id << " is a source\n";
      device.reset(new DataSourceDevice(spec, serviceRegistry));
    } else {
      LOG(DEBUG) << spec.id << " is a processor\n";
      device.reset(new DataProcessingDevice(spec, serviceRegistry));
    }

    serviceRegistry.get<RawDeviceService>().setDevice(device.get());

    runner.AddHook<fair::mq::hooks::InstantiateDevice>([&device](fair::mq::DeviceRunner& r) {
      r.fDevice = std::shared_ptr<FairMQDevice>{ std::move(device) };
      TurnOffColors::apply(false, 0);
    });

    return runner.Run();
  } catch (std::exception& e) {
    LOG(ERROR) << "Unhandled exception reached the top of main: " << e.what() << ", device shutting down.";
    return 1;
  } catch (...) {
    LOG(ERROR) << "Unknown exception reached the top of main.\n";
    return 1;
  }
  return 0;
}

// This is a toy executor for the workflow spec
// What it needs to do is:
//
// - Print the properties of each DataProcessorSpec
// - Fork one process per DataProcessorSpec
//   - Parent -> wait for all the children to complete (eventually
//     killing them all on ctrl-c).
//   - Child, pick the data-processor ID and start a O2DataProcessorDevice for
//     each DataProcessorSpec
int doMain(int argc, char** argv, const o2::framework::WorkflowSpec& specs,
           std::vector<ChannelConfigurationPolicy> const& channelPolicies)
{
  bpo::options_description executorOptions("Executor options");
  executorOptions.add_options()((std::string("help") + ",h").c_str(), "print this help")(
    (std::string("quiet") + ",q").c_str(), bpo::value<bool>()->zero_tokens()->default_value(false), "quiet operation")(
    (std::string("stop") + ",s").c_str(), bpo::value<bool>()->zero_tokens()->default_value(false),
    "stop before device start")((std::string("batch") + ",b").c_str(),
                                bpo::value<bool>()->zero_tokens()->default_value(false), "batch processing mode")(
    (std::string("graphviz") + ",g").c_str(), bpo::value<bool>()->zero_tokens()->default_value(false),
    "produce graph output")((std::string("dds") + ",D").c_str(),
                            bpo::value<bool>()->zero_tokens()->default_value(false), "create DDS configuration");

  // some of the options must be forwarded by default to the device
  executorOptions.add(DeviceSpecHelpers::getForwardedDeviceOptions());

  gHiddenDeviceOptions.add_options()((std::string("id") + ",i").c_str(), bpo::value<std::string>(),
                                     "device id for child spawning")(
    "channel-config", bpo::value<std::vector<std::string>>(), "channel configuration")("control", "control plugin")(
    "log-color", "logging color scheme")("color", "logging color scheme");

  bpo::options_description visibleOptions;
  visibleOptions.add(executorOptions);
  // Use the hidden options as veto, all config specs matching a definition
  // in the hidden options are skipped in order to avoid duplicate definitions
  // in the main parser. Note: all config specs are forwarded to devices
  visibleOptions.add(prepareOptionDescriptions(specs, gHiddenDeviceOptions));

  bpo::options_description od;
  od.add(visibleOptions);
  od.add(gHiddenDeviceOptions);

  // FIXME: decide about the policy for handling unrecognized arguments
  // command_line_parser with option allow_unregistered() can be used
  bpo::variables_map varmap;
  bpo::store(bpo::parse_command_line(argc, argv, od), varmap);

  bool defaultQuiet = varmap["quiet"].as<bool>();
  bool defaultStopped = varmap["stop"].as<bool>();
  bool batch = varmap["batch"].as<bool>();
  bool graphViz = varmap["graphviz"].as<bool>();
  bool generateDDS = varmap["dds"].as<bool>();
  std::string frameworkId;
  if (varmap.count("id"))
    frameworkId = varmap["id"].as<std::string>();
  if (varmap.count("help")) {
    bpo::options_description helpOptions;
    helpOptions.add(executorOptions);
    // this time no veto is applied, so all the options are added for printout
    helpOptions.add(prepareOptionDescriptions(specs));
    std::cout << helpOptions << std::endl;
    exit(0);
  }

  std::vector<DeviceSpec> deviceSpecs;

  try {
    DeviceSpecHelpers::dataProcessorSpecs2DeviceSpecs(specs, channelPolicies, deviceSpecs);
    // This should expand nodes so that we can build a consistent DAG.
  } catch (std::runtime_error& e) {
    std::cerr << "Invalid workflow: " << e.what() << std::endl;
    return 1;
  }

  // Up to here, parent and child need to do exactly the same thing. After, we
  // distinguish between something which has a framework id (the children) and something
  // which does not, the parent.
  if (frameworkId.empty() == false) {
    for (auto& spec : deviceSpecs) {
      if (spec.id == frameworkId) {
        return doChild(argc, argv, spec);
      }
    }
    LOG(ERROR) << "Unable to find component with id " << frameworkId;
  }

  assert(frameworkId.empty());

  gDeviceControls.resize(deviceSpecs.size());
  gDeviceExecutions.resize(deviceSpecs.size());

  DeviceSpecHelpers::prepareArguments(argc, argv, defaultQuiet, defaultStopped, deviceSpecs, gDeviceExecutions,
                                      gDeviceControls);

  std::vector<DriverState> startProgram;
  if (graphViz) {
    // Dump a graphviz representation of what I will do.
    GraphvizHelpers::dumpDeviceSpec2Graphviz(std::cout, deviceSpecs);
    exit(0);
  } else if (generateDDS) {
    dumpDeviceSpec2DDS(std::cout, deviceSpecs, gDeviceExecutions);
    exit(0);
  } else {
    // By default we simply start the main loop
    startProgram = {DriverState::INIT};
  }
 
  doParent(gDeviceInfos, deviceSpecs, gDeviceControls, gDeviceMetricsInfos,
           gDeviceExecutions,
           startProgram, batch);
  return killChildren(gDeviceInfos);
}
