// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ClustererTask.h
/// \brief Definition of the ITS cluster finder task

#ifndef ALICEO2_ENDCAPSLAYERS_CLUSTERERTASK
#define ALICEO2_ENDCAPSLAYERS_CLUSTERERTASK

#include "ECLayersBase/GeometryTGeo.h"
#include "EndCapsReconstruction/ChipMappingEC0.h"
#include "EndCapsReconstruction/PixelReader.h"
#include "EndCapsReconstruction/RawPixelReader.h"
#include "EndCapsReconstruction/DigitPixelReader.h"
#include "EndCapsReconstruction/Clusterer.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <memory>
#include <limits>

namespace o2
{
class MCCompLabel;
namespace dataformats
{
template <typename T>
class MCTruthContainer;
}

namespace ecl
{

class ClustererTask
{
  using Clusterer = o2::endcaps::Clusterer;
  using Cluster = o2::itsmft::Cluster;
  using CompCluster = o2::itsmft::CompCluster;
  using CompClusterExt = o2::itsmft::CompClusterExt;
  using MCTruth = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

 public:
  ClustererTask(bool useMC = true, bool raw = false);
  ~ClustererTask();

  void Init();
  Clusterer& getClusterer() { return mClusterer; }
  void run(const std::string inpName, const std::string outName);
  o2::endcaps::PixelReader* getReader() const { return (o2::endcaps::PixelReader*)mReader; }

  void loadDictionary(std::string fileName) { mClusterer.loadDictionary(fileName); }

  void writeTree(std::string basename, int i);
  void setMaxROframe(int max) { maxROframe = max; }
  int getMaxROframe() const { return maxROframe; }

 private:
  int maxROframe = std::numeric_limits<int>::max();                                   ///< maximal number of RO frames per a file
  bool mRawDataMode = false;                                                          ///< input from raw data or MC digits
  bool mUseMCTruth = true;                                                            ///< flag to use MCtruth if available
  o2::endcaps::PixelReader* mReader = nullptr;                                         ///< Pointer on the relevant Pixel reader
  std::unique_ptr<o2::endcaps::DigitPixelReader> mReaderMC;                            ///< reader for MC data
  std::unique_ptr<o2::endcaps::RawPixelReader<o2::endcaps::ChipMappingEC0>> mReaderRaw; ///< reader for raw data

  const o2::endcaps::GeometryTGeo* mGeometry = nullptr; ///< ITS OR MFT upgrade geometry
  Clusterer mClusterer;                                ///< Cluster finder

  std::vector<Cluster> mFullClus;               //!< vector of full clusters

  std::vector<CompClusterExt> mCompClus;               //!< vector of compact clusters

  std::vector<o2::itsmft::ROFRecord> mROFRecVec;               //!< vector of ROFRecord references

  MCTruth mClsLabels;               //! MC labels

  std::vector<unsigned char> mPatterns;

  ClassDefNV(ClustererTask, 2);
};
} // namespace ecl
} // namespace o2

#endif /* ALICEO2_ENDCAPSLAYERS_CLUSTERERTASK */
