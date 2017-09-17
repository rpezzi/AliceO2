// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_FIT_GEOMETRY_H_
#define ALICEO2_FIT_GEOMETRY_H_
////////////////////////////////////////////////
// Full geomrtry  hits classes for detector: FIT    //
////////////////////////////////////////////////

#include <Rtypes.h>
#include <TVector3.h>
namespace o2
{
namespace fit
{
// FIT is not tracking detector, Geometry could be used in future but not now. So just simple settings
class Geometry
{
 public:
  ///
  /// Default constructor.
  /// It must be kept public for root persistency purposes,
  /// but should never be called by the outside world
  Geometry();
  ///
  TVector3 centerMCP(int imcp) { return mMCP[imcp]; }
  static constexpr int NCellsA = 24; // number of radiatiors on A side
  static constexpr int NCellsC = 28; // number of radiatiors on C side

 private:
  TVector3 mMCP[52];

  ClassDefNV(Geometry, 1);
};
}
}
#endif
