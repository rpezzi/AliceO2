// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackFitter.cxx
/// \brief Implementation of a class to fit a track to a set of clusters
///
/// \author Philippe Pillot, Subatech; adapted by Rafael Pezzi, UFRGS

#include "MFTTracking/TrackFitter.h"


namespace o2
{
namespace mft
{

//_________________________________________________________________________________________________
void TrackFitter::initField(float l3Current)
{
  /// Set the magnetic field for the MFT
}

//_________________________________________________________________________________________________
void TrackFitter::fit(FitterTrackMFT& track, bool smooth, bool finalize,
                      std::list<TrackParam>::reverse_iterator* itStartingParam)
{
  /// Fit a track to its attached clusters
  /// Smooth the track if requested and the smoother enabled
  /// If finalize = true: copy the smoothed parameters, if any, into the regular ones
  /// Fit the entire track or only the part upstream itStartingParam
  /// Throw an exception in case of failure

  // initialize the starting track parameters and cluster
  auto itParam(track.rbegin());
  if (itStartingParam != nullptr) {
    // use the ones pointed to by itStartingParam
    if (*itStartingParam == track.rend()) {
      throw std::runtime_error("invalid starting parameters");
    }
    itParam = *itStartingParam;
  } else {
    // or start from the last cluster and compute the track parameters from its position
    // and the one of the first previous cluster found on a different layer
    auto itPreviousParam(itParam);
    ++itPreviousParam;
    // TODO: Refine criteria for initTrack: cluster from next MFT layer or next disk?
    (*itParam).setInverseMomentum(momentumFromSagitta(track)); // TODO: Estimate momentum from MFT track sagitta and chord;
    initTrack(*itPreviousParam->getClusterPtr(), *itParam->getClusterPtr(), *itParam);
  }

  // recusively add the upstream clusters and update the track parameters
  TrackParam* startingParam = &*itParam;
  while (++itParam != track.rend()) {
    try {
      addCluster(*startingParam, *itParam->getClusterPtr(), *itParam);
      startingParam = &*itParam;
    } catch (std::exception const&) {
      throw;
    }
  }

  // smooth the track if requested and the smoother enabled
  if (smooth && mSmooth) {
    try {
      smoothTrack(track, finalize);
    } catch (std::exception const&) {
      throw;
    }
  }
}

//_________________________________________________________________________________________________
void TrackFitter::initTrack(const o2::itsmft::Cluster& cl1, const o2::itsmft::Cluster& cl2, TrackParam& param)
{

  /// Compute the initial track parameters at the z position of the last cluster (cl2)
  /// The covariance matrix is computed such that the last cluster is the only constraint
  /// (by assigning an infinite dispersion to the other cluster)
  /// These parameters are the seed for the Kalman filter

  // compute the track parameters at the last cluster
  double dZ = cl1.getZ() - cl2.getZ();
  std::cout << "initTrack: clusters  at  cl1.getZ() = " << cl1.getZ() << " and cl2.getZ() = " << cl2.getZ() << std::endl;

  param.setX(cl2.getX());
  param.setY(cl2.getY());
  param.setZ(cl2.getZ());
  param.setXSlope((cl1.getX() - cl2.getX()) / dZ);
  param.setYSlope((cl1.getY() - cl2.getY()) / dZ);
  double impact = cl2.getY() - cl2.getZ() * param.getYSlope();
  double inverseMomentum = param.getInverseMomentum();
  //param.setInverseMomentum(inverseMomentum);

  // compute the track parameter covariances at the last cluster (as if the other clusters did not exist)
  TMatrixD lastParamCov(5, 5);
  lastParamCov.Zero();
  lastParamCov(0, 0) = cl2.getSigmaZ2(); // cl2.getEx2();
  lastParamCov(0, 1) = -lastParamCov(0, 0) / dZ;
  lastParamCov(1, 0) = lastParamCov(0, 1);
  lastParamCov(1, 1) = (1000. * cl1.getSigmaZ2() + lastParamCov(0, 0)) / dZ / dZ;

  lastParamCov(2, 2) = cl2.getSigmaY2(); // cl2.getEy2();
  lastParamCov(2, 3) = -lastParamCov(2, 2) / dZ;
  lastParamCov(3, 2) = lastParamCov(2, 3);
  lastParamCov(3, 3) = (1000. * lastParamCov(2, 2)) / dZ / dZ;
  lastParamCov(4, 4) = inverseMomentum * inverseMomentum;
  param.setCovariances(lastParamCov);

  // set other parameters
  param.setClusterPtr(&cl2);
  param.setTrackChi2(0.);
}

//_________________________________________________________________________________________________
void TrackFitter::addCluster(const TrackParam& startingParam, const o2::itsmft::Cluster& cl, TrackParam& param)
{
  /// Extrapolate the starting track parameters to the z position of the new cluster
  /// accounting for MCS dispersion in the current layer and the other(s) crossed
  /// Recompute the parameters adding the cluster constraint with the Kalman filter
  /// Throw an exception in case of failure

  if (cl.getZ() <= startingParam.getZ()) {
    throw std::runtime_error("The new cluster must be upstream");
  }
  //std::cout << "addCluster: add cluster at cl.getZ() = " << cl.getZ() << std::endl;
  // copy the current parameters into the new ones
  param.setParameters(startingParam.getParameters());
  param.setZ(startingParam.getZ());
  param.setCovariances(startingParam.getCovariances());
  param.setTrackChi2(startingParam.getTrackChi2());

  // add MCS effect in the current layer
  o2::itsmft::ChipMappingMFT mftChipMapper;
  int currentLayer(mftChipMapper.chip2Layer(startingParam.getClusterPtr()->getSensorID()));
  addMCSEffect(&param, SLayerThicknessInX0[currentLayer], -1.);

  // reset propagator for smoother
  if (mSmooth) {
    param.resetPropagator();
  }

  // add MCS in missing layers if any
  int expectedLayer(currentLayer - 1);
  currentLayer = mftChipMapper.chip2Layer(cl.getSensorID());
  while (currentLayer < expectedLayer) {
    if (!extrapToZCov(&param, SDefaultLayerZ[expectedLayer], mSmooth)) {
      throw std::runtime_error("Track extrapolation failed");
    }
    addMCSEffect(&param, SLayerThicknessInX0[expectedLayer], -1.);
    expectedLayer--;
  }

  // extrapolate to the z position of the new cluster
  if (!extrapToZCov(&param, cl.getZ(), mSmooth)) {
    throw std::runtime_error("Track extrapolation failed");
  }

  // save extrapolated parameters and covariances for smoother
  if (mSmooth) {
    param.setExtrapParameters(param.getParameters());
    param.setExtrapCovariances(param.getCovariances());
  }

  // recompute the parameters
  param.setClusterPtr(&cl);
  try {
    runKalmanFilter(param);
  } catch (std::exception const&) {
    throw;
  }
}

//_________________________________________________________________________________________________
void TrackFitter::smoothTrack(FitterTrackMFT& track, bool finalize)
{
  /// Recompute the track parameters at each cluster using the Smoother
  /// Smoothed parameters are stored in dedicated data members
  /// If finalize, they are copied in the regular parameters in case of success
  /// Throw an exception in case of failure

  auto itCurrentParam(track.begin());
  auto itPreviousParam(itCurrentParam);
  ++itCurrentParam;

  // smoothed parameters and covariances at first cluster = filtered parameters and covariances
  itPreviousParam->setSmoothParameters(itPreviousParam->getParameters());
  itPreviousParam->setSmoothCovariances(itPreviousParam->getCovariances());

  // local chi2 at first cluster = last additional chi2 provided by Kalman
  itPreviousParam->setLocalChi2(itPreviousParam->getTrackChi2() - itCurrentParam->getTrackChi2());

  // recursively smooth the next parameters and covariances
  do {
    try {
      runSmoother(*itPreviousParam, *itCurrentParam);
    } catch (std::exception const&) {
      throw;
    }
    ++itPreviousParam;
  } while (++itCurrentParam != track.end());

  // update the regular parameters and covariances if requested
  if (finalize) {
    for (auto& param : track) {
      param.setParameters(param.getSmoothParameters());
      param.setCovariances(param.getSmoothCovariances());
    }
  }
}

//_________________________________________________________________________________________________
void TrackFitter::runKalmanFilter(TrackParam& trackParam)
{
  /// Compute the new track parameters including the attached cluster with the Kalman filter
  /// The current parameters are supposed to have been extrapolated to the cluster z position
  /// Throw an exception in case of failure

  // get actual track parameters (p)
  TMatrixD param(trackParam.getParameters());

  // get new cluster parameters (m)
  const o2::itsmft::Cluster* cluster = trackParam.getClusterPtr();
  TMatrixD clusterParam(5, 1);
  clusterParam.Zero();
  clusterParam(0, 0) = cluster->getX();
  clusterParam(2, 0) = cluster->getY();

  // compute the actual parameter weight (W)
  TMatrixD paramWeight(trackParam.getCovariances());
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    throw std::runtime_error("1. Determinant = 0");
  }

  // compute the new cluster weight (U)
  TMatrixD clusterWeight(5, 5);
  clusterWeight.Zero();
  clusterWeight(0, 0) = 1. / cluster->getSigmaZ2(); // 1. / cluster->getEx2();
  clusterWeight(2, 2) = 1. / cluster->getSigmaY2(); //  1. / cluster->getEy2();

  // compute the new parameters covariance matrix ((W+U)^-1)
  TMatrixD newParamCov(paramWeight, TMatrixD::kPlus, clusterWeight);
  if (newParamCov.Determinant() != 0) {
    newParamCov.Invert();
  } else {
    throw std::runtime_error("2. Determinant = 0");
  }
  trackParam.setCovariances(newParamCov);

  // compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(clusterParam, TMatrixD::kMinus, param);   // m-p
  TMatrixD tmp2(clusterWeight, TMatrixD::kMult, tmp);    // U(m-p)
  TMatrixD newParam(newParamCov, TMatrixD::kMult, tmp2); // ((W+U)^-1)U(m-p)
  newParam += param;                                     // ((W+U)^-1)U(m-p) + p
  trackParam.setParameters(newParam);

  // compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam;                                                // p'
  tmp -= param;                                                  // (p'-p)
  TMatrixD tmp3(paramWeight, TMatrixD::kMult, tmp);              // W(p'-p)
  TMatrixD addChi2Track(tmp, TMatrixD::kTransposeMult, tmp3);    // ((p'-p)^-1)W(p'-p)
  tmp = newParam;                                                // p'
  tmp -= clusterParam;                                           // (p'-m)
  TMatrixD tmp4(clusterWeight, TMatrixD::kMult, tmp);            // U(p'-m)
  addChi2Track += TMatrixD(tmp, TMatrixD::kTransposeMult, tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  trackParam.setTrackChi2(trackParam.getTrackChi2() + addChi2Track(0, 0));
}

//_________________________________________________________________________________________________
void TrackFitter::runSmoother(const TrackParam& previousParam, TrackParam& param)
{
  /// Recompute the track parameters starting from the previous ones
  /// Throw an exception in case of failure

  // get variables
  const TMatrixD& extrapParameters = previousParam.getExtrapParameters();           // X(k+1 k)
  const TMatrixD& filteredParameters = param.getParameters();                       // X(k k)
  const TMatrixD& previousSmoothParameters = previousParam.getSmoothParameters();   // X(k+1 n)
  const TMatrixD& propagator = previousParam.getPropagator();                       // F(k)
  const TMatrixD& extrapCovariances = previousParam.getExtrapCovariances();         // C(k+1 k)
  const TMatrixD& filteredCovariances = param.getCovariances();                     // C(k k)
  const TMatrixD& previousSmoothCovariances = previousParam.getSmoothCovariances(); // C(k+1 n)

  // compute smoother gain: A(k) = C(kk) * F(k)^t * (C(k+1 k))^-1
  TMatrixD extrapWeight(extrapCovariances);
  if (extrapWeight.Determinant() != 0) {
    extrapWeight.Invert(); // (C(k+1 k))^-1
  } else {
    throw std::runtime_error("3. Determinant = 0");
  }
  TMatrixD smootherGain(filteredCovariances, TMatrixD::kMultTranspose, propagator); // C(kk) * F(k)^t
  smootherGain *= extrapWeight;                                                     // C(kk) * F(k)^t * (C(k+1 k))^-1

  // compute smoothed parameters: X(k n) = X(k k) + A(k) * (X(k+1 n) - X(k+1 k))
  TMatrixD tmpParam(previousSmoothParameters, TMatrixD::kMinus, extrapParameters); // X(k+1 n) - X(k+1 k)
  TMatrixD smoothParameters(smootherGain, TMatrixD::kMult, tmpParam);              // A(k) * (X(k+1 n) - X(k+1 k))
  smoothParameters += filteredParameters;                                          // X(k k) + A(k) * (X(k+1 n) - X(k+1 k))
  param.setSmoothParameters(smoothParameters);

  // compute smoothed covariances: C(k n) = C(k k) + A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
  TMatrixD tmpCov(previousSmoothCovariances, TMatrixD::kMinus, extrapCovariances); // C(k+1 n) - C(k+1 k)
  TMatrixD tmpCov2(tmpCov, TMatrixD::kMultTranspose, smootherGain);                // (C(k+1 n) - C(k+1 k)) * (A(k))^t
  TMatrixD smoothCovariances(smootherGain, TMatrixD::kMult, tmpCov2);              // A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
  smoothCovariances += filteredCovariances;                                        // C(k k) + A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
  param.setSmoothCovariances(smoothCovariances);

  // compute smoothed residual: r(k n) = cluster - X(k n)
  const o2::itsmft::Cluster* cluster = param.getClusterPtr();
  TMatrixD smoothResidual(2, 1);
  smoothResidual.Zero();
  smoothResidual(0, 0) = cluster->getX() - smoothParameters(0, 0);
  smoothResidual(1, 0) = cluster->getY() - smoothParameters(2, 0);

  // compute weight of smoothed residual: W(k n) = (clusterCov - C(k n))^-1
  TMatrixD smoothResidualWeight(2, 2);
  smoothResidualWeight(0, 0) = cluster->getSigmaZ2() - smoothCovariances(0, 0); // cluster->getEx2() - smoothCovariances(0, 0);
  smoothResidualWeight(0, 1) = -smoothCovariances(0, 2);
  smoothResidualWeight(1, 0) = -smoothCovariances(2, 0);
  smoothResidualWeight(1, 1) = cluster->getSigmaY2() - smoothCovariances(2, 2); // cluster->getEy2() - smoothCovariances(2, 2);
  if (smoothResidualWeight.Determinant() != 0) {
    smoothResidualWeight.Invert();
  } else {
    throw std::runtime_error("4. Determinant = 0");
  }

  // compute local chi2 = (r(k n))^t * W(k n) * r(k n)
  TMatrixD tmpChi2(smoothResidual, TMatrixD::kTransposeMult, smoothResidualWeight); // (r(k n))^t * W(k n)
  TMatrixD localChi2(tmpChi2, TMatrixD::kMult, smoothResidual);                     // (r(k n))^t * W(k n) * r(k n)
  param.setLocalChi2(localChi2(0, 0));
}

//__________________________________________________________________________
void linearExtrapToZ(TrackParam* trackParam, double zEnd)
{
  /// Track parameters linearly extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.

  if (trackParam->getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // Compute track parameters
  double dZ = zEnd - trackParam->getZ();
  trackParam->setX(trackParam->getX() + trackParam->getXSlope() * dZ);
  trackParam->setY(trackParam->getY() + trackParam->getYSlope() * dZ);
  trackParam->setZ(zEnd);
}

//__________________________________________________________________________
void linearExtrapToZCov(TrackParam* trackParam, double zEnd, bool updatePropagator)
{
  /// Track parameters and their covariances linearly extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.

  if (trackParam->getZ() == zEnd) {
    return; // nothing to be done if same z
  }

  // No need to propagate the covariance matrix if it does not exist
  if (!trackParam->hasCovariances()) {
    LOG(WARNING) << "Covariance matrix does not exist";
    // Extrapolate linearly track parameters to "zEnd"
    linearExtrapToZ(trackParam, zEnd);
    return;
  }

  // Compute track parameters
  double dZ = zEnd - trackParam->getZ();
  trackParam->setX(trackParam->getX() + trackParam->getXSlope() * dZ);
  trackParam->setY(trackParam->getY() + trackParam->getYSlope() * dZ);
  trackParam->setZ(zEnd);

  // Calculate the jacobian related to the track parameters linear extrapolation to "zEnd"
  TMatrixD jacob(5, 5);
  jacob.UnitMatrix();
  jacob(0, 1) = dZ;
  jacob(2, 3) = dZ;

  // Extrapolate track parameter covariances to "zEnd"
  TMatrixD tmp(trackParam->getCovariances(), TMatrixD::kMultTranspose, jacob);
  TMatrixD tmp2(jacob, TMatrixD::kMult, tmp);
  trackParam->setCovariances(tmp2);

  // Update the propagator if required
  if (updatePropagator) {
    trackParam->updatePropagator(jacob);
  }
}

//__________________________________________________________________________
bool extrapToZ(TrackParam* trackParam, double zEnd, bool isFieldON)
{
  /// Interface to track parameter extrapolation to the plane at "Z".
  /// On return, the track parameters resulting from the extrapolation are updated in trackParam.
  if (!isFieldON) {
    linearExtrapToZ(trackParam, zEnd);
    return true;
  } else {
    linearExtrapToZ(trackParam, zEnd); // FIXME: Add proper helix extrapolation
    return true;
  }
}

//__________________________________________________________________________
bool extrapToZCov(TrackParam* trackParam, double zEnd, bool updatePropagator, bool isFieldON)
{
  /// Track parameters and their covariances extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.

  //++sNCallExtrapToZCov;

  if (trackParam->getZ() == zEnd) {
    return true; // nothing to be done if same z
  }

  if (!isFieldON) { // linear extrapolation if no magnetic field
    linearExtrapToZCov(trackParam, zEnd, updatePropagator);
    return true;
  }

  // No need to propagate the covariance matrix if it does not exist
  if (!trackParam->hasCovariances()) {
    LOG(WARNING) << "Covariance matrix does not exist";
    // Extrapolate track parameters to "zEnd"
    return extrapToZ(trackParam, zEnd);
  }

  // Save the actual track parameters
  TrackParam trackParamSave(*trackParam);
  TMatrixD paramSave(trackParamSave.getParameters());
  double zBegin = trackParamSave.getZ();

  // Get reference to the parameter covariance matrix
  const TMatrixD& kParamCov = trackParam->getCovariances();

  // Extrapolate track parameters to "zEnd"
  // Do not update the covariance matrix if the extrapolation failed
  if (!extrapToZ(trackParam, zEnd)) {
    return false;
  }

  // Get reference to the extrapolated parameters
  const TMatrixD& extrapParam = trackParam->getParameters();

  // Calculate the jacobian related to the track parameters extrapolation to "zEnd"
  TMatrixD jacob(5, 5);
  jacob.Zero();
  TMatrixD dParam(5, 1);
  double direction[5] = {-1., -1., 1., 1., -1.};
  for (int i = 0; i < 5; i++) {
    // Skip jacobian calculation for parameters with no associated error
    if (kParamCov(i, i) <= 0.) {
      continue;
    }

    // Small variation of parameter i only
    for (int j = 0; j < 5; j++) {
      if (j == i) {
        dParam(j, 0) = TMath::Sqrt(kParamCov(i, i));
        dParam(j, 0) *= TMath::Sign(1., direction[j] * paramSave(j, 0)); // variation always in the same direction
      } else {
        dParam(j, 0) = 0.;
      }
    }

    // Set new parameters
    trackParamSave.setParameters(paramSave);
    trackParamSave.addParameters(dParam);
    trackParamSave.setZ(zBegin);

    // Extrapolate new track parameters to "zEnd"
    if (!extrapToZ(&trackParamSave, zEnd)) {
      LOG(WARNING) << "Bad covariance matrix";
      return false;
    }

    // Calculate the jacobian
    TMatrixD jacobji(trackParamSave.getParameters(), TMatrixD::kMinus, extrapParam);
    jacobji *= 1. / dParam(i, 0);
    jacob.SetSub(0, i, jacobji);
  }

  // Extrapolate track parameter covariances to "zEnd"
  TMatrixD tmp(kParamCov, TMatrixD::kMultTranspose, jacob);
  TMatrixD tmp2(jacob, TMatrixD::kMult, tmp);
  trackParam->setCovariances(tmp2);

  // Update the propagator if required
  if (updatePropagator) {
    trackParam->updatePropagator(jacob);
  }

  return true;
}

//__________________________________________________________________________
void addMCSEffect(TrackParam* trackParam, double dZ, double x0, bool isFieldON)
{
  /// Add to the track parameter covariances the effects of multiple Coulomb scattering
  /// through a material of thickness "abs(dZ)" and of radiation length "x0"
  /// assuming linear propagation and using the small angle approximation.
  /// dZ = zOut - zIn (sign is important) and "param" is assumed to be given zOut.
  /// If x0 <= 0., assume dZ = pathLength/x0 and consider the material thickness as negligible.

  double ySlope = trackParam->getYSlope();
  double xSlope = trackParam->getXSlope();
  double inverseMomentum = trackParam->getInverseMomentum();
  double inverseTotalMomentum2 = inverseMomentum * inverseMomentum * (1.0 + ySlope * ySlope) /
                                 (1.0 + ySlope * ySlope + xSlope * xSlope);
  // Path length in the material
  double signedPathLength = dZ * TMath::Sqrt(1.0 + ySlope * ySlope + xSlope * xSlope);
  double pathLengthOverX0 = (x0 > 0.) ? TMath::Abs(signedPathLength) / x0 : TMath::Abs(signedPathLength);
  // relativistic velocity
  double velo = 1.;
  // Angular dispersion square of the track (variance) in a plane perpendicular to the trajectory
  double theta02 = 0.0136 / velo * (1 + 0.038 * TMath::Log(pathLengthOverX0));
  theta02 *= theta02 * inverseTotalMomentum2 * pathLengthOverX0;

  double varCoor = (x0 > 0.) ? signedPathLength * signedPathLength * theta02 / 3. : 0.;
  double varSlop = theta02;
  double covCorrSlope = (x0 > 0.) ? signedPathLength * theta02 / 2. : 0.;

  // Set MCS covariance matrix
  TMatrixD newParamCov(trackParam->getCovariances());
  // Non bending plane
  newParamCov(0, 0) += varCoor;
  newParamCov(0, 1) += covCorrSlope;
  newParamCov(1, 0) += covCorrSlope;
  newParamCov(1, 1) += varSlop;
  // Bending plane
  newParamCov(2, 2) += varCoor;
  newParamCov(2, 3) += covCorrSlope;
  newParamCov(3, 2) += covCorrSlope;
  newParamCov(3, 3) += varSlop;

  // Set momentum related covariances if B!=0
  if (isFieldON) {
    // compute derivative d(q/Pxy) / dSlopeX and d(q/Pxy) / dSlopeY
    double dqPxydSlopeX =
      inverseMomentum * xSlope / (1. + xSlope * xSlope + ySlope * ySlope);
    double dqPxydSlopeY = -inverseMomentum * xSlope * xSlope * ySlope /
                          (1. + ySlope * ySlope) /
                          (1. + xSlope * xSlope + ySlope * ySlope);
    // Inverse bending momentum (due to dependences with bending and non bending slopes)
    newParamCov(4, 0) += dqPxydSlopeX * covCorrSlope;
    newParamCov(0, 4) += dqPxydSlopeX * covCorrSlope;
    newParamCov(4, 1) += dqPxydSlopeX * varSlop;
    newParamCov(1, 4) += dqPxydSlopeX * varSlop;
    newParamCov(4, 2) += dqPxydSlopeY * covCorrSlope;
    newParamCov(2, 4) += dqPxydSlopeY * covCorrSlope;
    newParamCov(4, 3) += dqPxydSlopeY * varSlop;
    newParamCov(3, 4) += dqPxydSlopeY * varSlop;
    newParamCov(4, 4) += (dqPxydSlopeX * dqPxydSlopeX + dqPxydSlopeY * dqPxydSlopeY) * varSlop;
  }

  // Set new covariances
  trackParam->setCovariances(newParamCov);
}

//__________________________________________________________________________
Double_t momentumFromSagitta(FitterTrackMFT& track, Double_t bFieldZ)
{
  std::cout << "   Running momentumFromSagitta for track with nclusters = " << track.getNClusters()  << std::endl;

  auto nPoints = track.getNClusters();
  Double_t* x = new Double_t[nPoints];
  Double_t* y = new Double_t[nPoints];
  TVectorD fitparam(3);
  int n = 0;
  for (auto trackparam = track.begin(); trackparam != track.end(); trackparam++) {
    x[n] = trackparam->getClusterPtr()->getX();
    y[n] = trackparam->getClusterPtr()->getY();
    std::cout << "    adding point to fit for z = " << trackparam->getClusterPtr()->getZ() << " (" << x[n] <<  "," << y[n]  <<  ") "<< std::endl;
    n++;
  }
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << " Fit QuadraticRegression: " << std::endl;
  Double_t q0, q1, q2;
  auto chi2 = QuadraticRegression(nPoints,x ,y, q0, q1,q2);
  std::cout << " Fit Parameters [0] = " << q0 << " [1] =  " << q1 << " [2] = " << q2 << std::endl;
  auto dz = track.last().getZ() - track.first().getZ();
  auto slopeX = (track.last().getX() - track.first().getX()) / dz;
  auto slopeY = (track.last().getY() - track.first().getY()) / dz;
  std::cout << " Radius from QuadraticRegression = " << 0.5/q2 << std::endl;
  auto pt = 0.3*bFieldZ*0.5/q2; // radius = 0.5/q2
  auto p = pt*sqrt(1.0+1.0/(slopeX*slopeX+slopeY*slopeY));
  std::cout << " Seed pt = " << pt << std::endl;
  std::cout << " Seed momentum = " << p << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << " Fit CircleRegression: " << std::endl;
  Double_t radiusCirclereg = CircleRegression(nPoints,x ,y);
  std::cout << " Radius from CircleRegression = " << radiusCirclereg << std::endl;

  std::cout << " Seed pt = " << 0.3*bFieldZ*radiusCirclereg << std::endl;
  std::cout << " Seed momentum = " << radiusCirclereg*sqrt(1.0+1.0/(slopeX*slopeX+slopeY*slopeY)) << std::endl;
  Double_t dist, qq;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << " Fit Sagitta: " << std::endl;

  Sagitta(nPoints,x ,y, dist, qq);
  std::cout << "--------------------------------------------" << std::endl;


  return p;

}

//__________________________________________________________________________
Double_t Sagitta(Int_t nVal, Double_t *xVal, Double_t *yVal, Double_t &distL2, Double_t &q2, double bFieldZ)
{
	/// Calculate sagitta of the track
	/// Return sagitta
	std::cout << "Sagitta" << std::endl;
	Double_t q0,q1,chi2;
	// fit by a polynom of 2nd order
	chi2 = QuadraticRegression(nVal,xVal ,yVal, q0, q1,q2);
	std::cout << Form("q param = %f  %f  %f  chi2 = %e ",q0, q1,q2,chi2) << std::endl;
	std::cout << Form("pt from 2nd order parameter %f ",0.01/q2/2.*0.3*bFieldZ/10.) << std::endl;


	Double_t maxDist = 0.;
	Int_t i1 = -1;
	Int_t i2 = -1;
	Double_t sagitta = 0.;
	Double_t dist ;
	Double_t distTot = 0.;

	for (int i=0; i<nVal-1; i++) {
		distTot += TMath::Sqrt((xVal[i]-xVal[i+1]) * (xVal[i]-xVal[i+1]) + (yVal[i]-yVal[i+1]) * (yVal[i]-yVal[i+1]));

		for (int j = i+1; j<nVal; j++) {
			Double_t dist = (xVal[i]-xVal[j]) * (xVal[i]-xVal[j]) + (yVal[i]-yVal[j]) * (yVal[i]-yVal[j]);
			if(dist>maxDist){
				i1 = i;
				i2 = j;
				maxDist = dist;
			}
		}
	}

	std::cout << Form("distMAx = %f distTot =%f   i1 = %d i2 =%d",TMath::Sqrt(maxDist),distTot,i1,i2) << std::endl;

	distL2 = 1.e-2*TMath::Sqrt((xVal[i1]-xVal[i2]) * (xVal[i1]-xVal[i2]) + (yVal[i1]-yVal[i2]) * (yVal[i1]-yVal[i2]) );

	Double_t p1 = (yVal[i1]-yVal[i2]) / (xVal[i1]-xVal[i2]);
	Double_t p0 = yVal[i2] - xVal[i2] * p1;
	std::cout << Form("p param = %f  %f  ",p0, p1) << std::endl;
	q2 = TMath::Sign(q2,q1*q2*(-(p0 + p1 * (xVal[i1]+xVal[i2])/2.))*(bFieldZ));


	Double_t y12 = q0+q1*xVal[i1]+q2*xVal[i1]*xVal[i1];
	Double_t y22 = q0+q1*xVal[i2]+q2*xVal[i2]*xVal[i2];
	std::cout << Form("x1, y1 %f  %f  ",xVal[i1], y12) << std::endl;
	std::cout << Form("x2, y2 %f  %f  ",xVal[i2], y22) << std::endl;

//	Double_t p12 = (y12-y22) / (xVal[i1]-xVal[i2]);
//	Double_t p02 = y22 - xVal[i2] * p12;
//	std::cout << Form("p2 param = %f  %f  ",p02, p12) << std::endl;

	Double_t maxSagitta = 0.;
	for (int i=0; i<20; i++) {
		Double_t step =TMath::Abs(xVal[i1]-xVal[i2])/20.;
		Double_t xmiddle = xVal[i1] + i*step;

		Double_t y1 = p0 + p1 * xmiddle;
		Double_t p1perp = -1./p1;
		Double_t p0perp = p0 + xmiddle *(p1-p1perp);
		//std::cout << Form("p param perp = %f  %f  ",p0perp, p1perp) << std::endl;

		Double_t xmax = (p1perp - q1 + TMath::Sqrt(p1perp*p1perp - 2*p1perp*q1 + q1*q1 + 4*p0perp*q2 - 4*q0*q2))/(2*q2);
		Double_t xmax2 = -(q1 -p1perp + TMath::Sqrt(p1perp*p1perp - 2*p1perp*q1 + q1*q1 + 4*p0perp*q2 - 4*q0*q2))/(2*q2);
		//		std::cout << Form("xmax = %f   xmax2  %f",xmax,xmax2) << std::endl;

		if (TMath::Abs(xmax2-xmiddle) < TMath::Abs(xmax-xmiddle)) xmax = xmax2;


		Double_t y2 = q0 + q1 * xmax  + q2 * xmax * xmax;

		sagitta = 1e-2*TMath::Sqrt((xmiddle-xmax) * (xmiddle-xmax) + (y1-y2) * (y1-y2) );

		//		std::cout << Form("sagitta = %f   Bz = %f",sagitta,bFieldZ) << std::endl;

		sagitta = TMath::Sign(sagitta,q1*q2*(-y1)*(bFieldZ));

		if(TMath::Abs(sagitta)>TMath::Abs(maxSagitta)) maxSagitta = sagitta;

	}
	std::cout << Form(" Max sagitta = %e => pt = %f",maxSagitta, 0.3*0.5*distL2*distL2/8./maxSagitta) << std::endl;

	return maxSagitta;
}


//__________________________________________________________________________
Double_t QuadraticRegression(Int_t nVal, Double_t *xVal, Double_t *yVal, Double_t &p0, Double_t &p1, Double_t &p2)
{
	/// Perform a Quadratic Regression
	/// Assume same error on all clusters = 1
	/// Return ~ Chi2

	TMatrixD y(nVal,1);
	TMatrixD x(nVal,3);
	TMatrixD xtrans(3,nVal);

	for (int i=0; i<nVal; i++) {
		y(i,0) = yVal[i];
		x(i,0) = 1.;
		x(i,1) = xVal[i];
		x(i,2) = xVal[i]*xVal[i];
		xtrans(0,i) = 1.;
		xtrans(1,i) = xVal[i];
		xtrans(2,i) = xVal[i]*xVal[i];
	}
	TMatrixD tmp(xtrans,TMatrixD::kMult,x);
	tmp.Invert();

	TMatrixD tmp2(xtrans,TMatrixD::kMult,y);
	TMatrixD b(tmp,TMatrixD::kMult,tmp2);

	p0 = b(0,0);
	p1 = b(1,0);
	p2 = b(2,0);

	// chi2 = (y-xb)^t . W . (y-xb)
	TMatrixD tmp3(x,TMatrixD::kMult,b);
	TMatrixD tmp4(y,TMatrixD::kMinus,tmp3);
	TMatrixD chi2(tmp4,TMatrixD::kTransposeMult,tmp4);


	return chi2(0,0);
}

//__________________________________________________________________________
Double_t CircleRegression(Int_t nVal, Double_t *xVal, Double_t *yVal)
{
	/// Perform a Circular Regression
	/// Assume same error on all clusters = 1
	/// Return Radius evaluation
	Double_t sumxi2 =0., sumxi =0., sumxiyi =0., sumyi =0.,sumyi2 =0., sumxi3 =0., sumyi3 =0.;
	Double_t sumxi2yi=0., sumxiyi2=0.;
	Double_t xi,yi, ri;
	for (int i=0; i<nVal; i++) {
		xi = xVal[i]/100.;
		yi = yVal[i]/100.;
		ri = TMath::Sqrt(xi*xi + yi*yi);
//		xi /= ri*ri;
//		yi /= ri*ri;

		sumxi += xi;
		sumyi += yi;
		sumxi2 += xi*xi;
		sumyi2 += yi*yi;
		sumxi3 += xi*xi*xi;
		sumyi3 += yi*yi*yi;
		sumxiyi += xi*yi;
		sumxi2yi += xi*xi*yi;
		sumxiyi2 += xi*yi*yi;
	}

	Double_t A = nVal * sumxi2 - sumxi*sumxi;
	Double_t B = nVal * sumxiyi - sumxi*sumyi;
	Double_t C = nVal * sumyi2 - sumyi*sumyi;
	Double_t D = 0.5*(nVal*sumxiyi2 -sumxi*sumyi2  +nVal*sumxi3 -sumxi*sumxi2);
	Double_t E = 0.5*(nVal*sumxi2yi -sumyi*sumxi2  +nVal*sumyi3 -sumyi*sumyi2);

	Double_t aM = (D*C-B*E) / (A*C-B*B);
	Double_t bM = (A*E-B*D) / (A*C-B*B);

	Double_t rM = 0.;
	Double_t rK = 0.;

	for (int i=0; i<nVal; i++) {
		xi = xVal[i]/100.;
		yi = yVal[i]/100.;

		rM += TMath::Sqrt( (xi-aM)*(xi-aM) + (yi-bM)*(yi-bM) );
		rK +=  ((xi-aM)*(xi-aM) + (yi-bM)*(yi-bM) );
	}
	rM /= nVal;
	rK = TMath::Sqrt(rK/nVal);

	std::cout << Form("aM %f bM %f rM %f rK %f => pt = %f or %f ",aM,bM,rM,rK,rM*0.3*0.5, rK*0.3*0.5) << std::endl;

	return (rM+rK)/2.;
}

} // namespace mft
} // namespace o2
