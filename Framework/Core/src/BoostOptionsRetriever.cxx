// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/BoostOptionsRetriever.h"
#include "Framework/ConfigParamSpec.h"
#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>

using namespace o2::framework;
namespace bpo = boost::program_options;

namespace o2 {
namespace framework {

BoostOptionsRetriever::BoostOptionsRetriever(std::vector<ConfigParamSpec> &specs, bool ignoreUnknown, int &argc, char **&argv)
: mVariables{},
  mDescription{"ALICE O2 Framework - Available options"},
  mIgnoreUnknown{ignoreUnknown}
{
  auto options = mDescription.add_options();
  for (auto & spec : specs) {
    const char *name = spec.name.c_str();
    const char *help = spec.help.c_str();
    // FIXME: propagate default value?
    switch(spec.type) {
      case VariantType::Int:
      case VariantType::Int64:
        options = options(name, bpo::value<int>()->default_value(spec.defaultValue.get<int>()), help);
        break;
      case VariantType::Float:
        options = options(name, bpo::value<float>()->default_value(spec.defaultValue.get<float>()), help);
        break;
      case VariantType::Double:
        options = options(name, bpo::value<double>()->default_value(spec.defaultValue.get<double>()), help);
        break;
      case VariantType::String:
        options = options(name, bpo::value<std::string>()->default_value(spec.defaultValue.get<std::string>()), help);
        break;
      case VariantType::Bool:
        options = options(name, bpo::value<bool>()->default_value(spec.defaultValue.get<bool>()), help);
        break;
      case VariantType::Unknown:
      case VariantType::Empty:
        break;
    };
  }
  parseArgs(argc, argv);
}

void BoostOptionsRetriever::parseArgs(int &argc, char **&argv) {
  if (mIgnoreUnknown == false) {
    auto parsed = bpo::parse_command_line(argc, argv, mDescription);
    bpo::store(parsed, mVariables);
    bpo::notify(mVariables);
    return;
  }
  auto parsed = bpo::command_line_parser(argc, argv).options(mDescription).allow_unregistered().run();
  bpo::store(parsed, mVariables);

  bpo::notify(mVariables);
}

int BoostOptionsRetriever::getInt(const char *key) const {
  return mVariables[key].as<int>();
}

float BoostOptionsRetriever::getFloat(const char *key) const {
  return mVariables[key].as<float>();
}

double BoostOptionsRetriever::getDouble(const char *key) const {
  return mVariables[key].as<double>();
}

bool BoostOptionsRetriever::getBool(const char *key) const {
  return mVariables[key].as<bool>();
}

std::string BoostOptionsRetriever::getString(const char *key) const {
  return mVariables[key].as<std::string>();
}

std::vector<std::string> BoostOptionsRetriever::getVString(const char *key) const {
  return mVariables[key].as<std::vector<std::string>>();
}

} // namespace framework
} // namespace o2
