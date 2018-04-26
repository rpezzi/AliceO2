// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HitDriftFiler.h
/// \brief Utility functions to filter hits within a given time window
/// \author sandro.wenzel@cern.ch

#ifndef ALICEO2_TPC_HitDriftFilter_H_
#define ALICEO2_TPC_HitDriftFilter_H_

namespace o2
{
namespace TPC
{

template <typename Collection>
void getHits(TChain& chain, const Collection& eventrecords, std::vector<std::vector<o2::TPC::HitGroup>*>& hitvectors,
             std::vector<o2::TPC::TPCHitGroupID>& hitids, const char* branchname, float tmin /*NS*/, float tmax /*NS*/,
             std::function<float(float, float, float)>&& f)
{
  // f is some function taking event time + z of hit and returns final "digit" time
  LOG(DEBUG) << "BR NAME " << branchname;
  auto br = chain.GetBranch(branchname);
  if (!br) {
    LOG(ERROR) << "No branch found";
    return;
  }

  auto nentries = br->GetEntries();

  // do the filtering
  for (int entry = 0; entry < nentries; ++entry) {
    if (tmin > f(eventrecords[entry].timeNS, 0, 0)) {
      continue;
    }
    if (tmax < f(eventrecords[entry].timeNS, 0, 250)) {
      break;
    }

    // This needs to be done only once for any entry
    if (hitvectors[entry] == nullptr) {
      br->SetAddress(&hitvectors[entry]);
      br->GetEntry(entry);
    }

    int groupid = -1;
    auto groups = hitvectors[entry];
    for (auto& singlegroup : *groups) {
      if (singlegroup.getSize() == 0) {
        // there are not hits in this group .. so continue
        // TODO: figure out why such a group would exist??
        continue;
      }
      const auto& pos = singlegroup.getHit(0).getPos();
      // std::cout << "This Group is in sector " << o2::TPC::Sector::ToSector(pos.X(), pos.Y(), pos.Z()) << "\n";
      groupid++;
      auto zmax = singlegroup.mZAbsMax;
      auto zmin = singlegroup.mZAbsMin;
      // in case of secondaries, the time ordering may be reversed
      if (zmax < zmin) {
        std::swap(zmax, zmin);
      }
      // auto tof = singlegroup.
      float tmaxtrack = f(eventrecords[entry].timeNS, 0., zmin);
      float tmintrack = f(eventrecords[entry].timeNS, 0., zmax);
      if (tmin > tmaxtrack || tmax < tmintrack) {
        // std::cout << "DISCARDING " << groupid << " OF ENTRY " << entry << "\n";
        continue;
      }
      // need to record index of the group
      hitids.emplace_back(entry, groupid);
    }
  }
}

// TPC hit selection lambda
auto calcDriftTime = [](float tNS, float tof, float z) {
  // returns time in NS
  return tNS + o2::TPC::ElectronTransport::getDriftTime(z) * 1000 + tof;
};

} // end namespace TPC
} // end namespace o2

#endif
