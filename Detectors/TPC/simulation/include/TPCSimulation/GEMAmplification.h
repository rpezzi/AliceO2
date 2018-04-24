// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GEMAmplification.h
/// \brief Definition of the GEM amplification
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ALICEO2_TPC_GEMAmplification_H_
#define ALICEO2_TPC_GEMAmplification_H_

#include "TPCBase/RandomRing.h"
#include "TPCBase/ParameterGas.h"
#include "TPCBase/ParameterGEM.h"
#include "TPCBase/CRU.h"
#include "TPCBase/PadPos.h"
#include "TPCBase/CalDet.h"

namespace o2
{
namespace TPC
{

/// \class GEMAmplification
/// This class handles the amplification of electrons in the GEM stack
/// The full amplification in a stack of four GEMs can be conducted, or each of the individual processes (Electrons
/// collection, amplification and extraction) can be conducted individually

class GEMAmplification
{
 public:
  /// Default constructor
  static GEMAmplification& instance()
  {
    static GEMAmplification gemAmplification;
    return gemAmplification;
  }

  /// Destructor
  ~GEMAmplification();

  /// Update the OCDB parameters cached in the class. To be called once per event
  void updateParameters();

  /// Compute the number of electrons after amplification in a full stack of four GEM foils
  /// \param nElectrons Number of electrons arriving at the first amplification stage (GEM1)
  /// \return Number of electrons after amplification in a full stack of four GEM foils
  int getStackAmplification(int nElectrons = 1);

  /// Compute the number of electrons after amplification in a full stack of four GEM foils
  /// taking into account local variations of the electron amplification
  /// \param nElectrons Number of electrons arriving at the first amplification stage (GEM1)
  /// \param cru CRU where the electron arrives
  /// \param pos PadPos where the electron arrives
  /// \return Number of electrons after amplification in a full stack of four GEM foils
  int getStackAmplification(const CRU& cru, const PadPos& pos, int nElectrons = 1);

  /// Compute the number of electrons after amplification in a single GEM foil
  /// taking into account collection and extraction efficiencies and fluctuations of the GEM amplification
  /// \param nElectrons Number of electrons to be amplified
  /// \param GEM Number of the GEM in the stack (1, 2, 3, 4)
  /// \return Number of electrons after amplification in a single GEM foil
  int getSingleGEMAmplification(int nElectrons, int GEM);

  /// Compute the electron losses due to extraction or collection efficiencies
  /// \param nElectrons Input number of electrons
  /// \param probability Collection or extraction efficiency
  /// \return Number of electrons after probable losses
  int getElectronLosses(int nElectrons, float probability);

  /// Compute the number of electrons after amplification in a single GEM foil
  /// taking into account avalanche fluctuations (Polya for <500 electrons and Gaus (central limit theorem) for a
  /// larger number of electrons)
  /// \param nElectrons Input number of electrons
  /// \param GEM Number of the GEM in the stack (1, 2, 3, 4)
  /// \return Number of electrons after amplification in the GEM
  int getGEMMultiplication(int nElectrons, int GEM);

 private:
  GEMAmplification();

  /// Circular random buffer containing random Gaus values for gain fluctuation if the number of electrons is larger
  /// (central limit theorem)
  RandomRing mRandomGaus;
  /// Circular random buffer containing flat random values for the collection/extraction
  RandomRing mRandomFlat;
  /// Container with random Polya distributions, one for each GEM in the stack
  std::array<RandomRing, 4> mGain;

  const ParameterGEM* mGEMParam; ///< Caching of the parameter class to avoid multiple CDB calls
  const ParameterGas* mGasParam; ///< Caching of the parameter class to avoid multiple CDB calls
  const CalPad* mGainMap;        ///< Caching of the parameter class to avoid multiple CDB calls
};
}
}

#endif // ALICEO2_TPC_GEMAmplification_H_
