// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  FastMultEst.h
/// \brief Fast multiplicity estimator for ITS
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_ENDCAPSLAYERS_FASTMULTEST_
#define ALICEO2_ENDCAPSLAYERS_FASTMULTEST_

#include "EndCapsReconstruction/ChipMappingEC0.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "EC0Reconstruction/FastMultEstConfig.h"
#include <gsl/span>
#include <array>

namespace o2
{
namespace ecl
{

struct FastMultEst {

  static constexpr int NLayers = o2::endcaps::ChipMappingEC0::NLayers;

  float mult = 0.;                         /// estimated signal clusters multipliciy at reference (1st?) layer
  float noisePerChip = 0.;                 /// estimated or imposed noise per chip
  float cov[3] = {0.};                     /// covariance matrix of estimation
  float chi2 = 0.;                         /// chi2
  int nLayersUsed = 0;                     /// number of layers actually used
  std::array<int, NLayers> nClPerLayer{0}; // measured N Cl per layer

  void fillNClPerLayer(const gsl::span<const o2::itsmft::CompClusterExt>& clusters);
  float process(const std::array<int, NLayers> ncl)
  {
    return FastMultEstConfig::Instance().imposeNoisePerChip > 0 ? processNoiseImposed(ncl) : processNoiseFree(ncl);
  }
  float processNoiseFree(const std::array<int, NLayers> ncl);
  float processNoiseImposed(const std::array<int, NLayers> ncl);
  float process(const gsl::span<const o2::itsmft::CompClusterExt>& clusters)
  {
    fillNClPerLayer(clusters);
    return process(nClPerLayer);
  }

  ClassDefNV(FastMultEst, 1);
};

} // namespace ecl
} // namespace o2

#endif
