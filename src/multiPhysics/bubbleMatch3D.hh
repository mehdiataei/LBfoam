/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BUBBLE_MATCH_3D_HH
#define BUBBLE_MATCH_3D_HH

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiPhysics/bubbleMatch3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

template <typename T>
void BubbleMatch3D::execute(MultiScalarField3D<int> &flag,
                            MultiScalarField3D<T> &volumeFraction) {
  bubbleBucketFill(flag);
  pluint numBubbles = countAndTagBubbles();
  bubbleVolume.clear();
  bubbleCenter.clear();
  bubbleDensity.clear();
  bubbleDensity.resize(numBubbles);
  std::fill(bubbleDensity.begin(), bubbleDensity.end(), (T)0.);
  bubbleAnalysis(flag, volumeFraction, numBubbles);
}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void BubbleMatch3D::executeWithGrowth(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag,
    MultiTensorField3D<T, 3> &normal, MultiScalarField3D<T> &rhoBar,
    MultiScalarField3D<T> &mass, MultiTensorField3D<T, 3> &j,
    MultiBlockLattice3D<T, ScalarDescriptor> &adLattice,
    MultiScalarField3D<T> &oldVolumeFraction, bool incompressibleModel, T pi) {
  bubbleBucketFill(flag);
  pluint numBubbles = countAndTagBubbles();
  bubbleVolume.clear();
  bubbleCenter.clear();
  bubbleDensity.clear();
  bubbleDisjoiningPressure.clear();
  bubbleAnalysis(flag, volumeFraction, numBubbles);
  calculateBubbleGrowth(flag, volumeFraction, oldVolumeFraction, adLattice,
                        numBubbles);
  calculateDisjoiningPressure<T, FluidDescriptor>(
      pi, flag, volumeFraction, normal, rhoBar, mass, j, incompressibleModel,
      numBubbles);
  computeBubbleData(numBubbles);
}

template <typename T, template <typename U1> class AD_Descriptor>
void BubbleMatch3D::calculateBubbleGrowth(
    MultiScalarField3D<int> &flag, MultiScalarField3D<T> &volumeFraction,
    MultiScalarField3D<T> &oldVolumeFraction,
    MultiBlockLattice3D<T, AD_Descriptor> &adLattice, pluint numBubbles)

{
  std::vector<MultiBlock3D *> args;
  args.push_back(tagMatrix);
  args.push_back(&flag);
  args.push_back(bubbleAnalysisContainer);
  args.push_back(&volumeFraction);
  args.push_back(&oldVolumeFraction);
  args.push_back(&adLattice);
  applyProcessingFunctional(
      new CalculateBubbleGrowth3D<T, AD_Descriptor>(numBubbles),
      bubbleAnalysisContainer->getBoundingBox(), args);
}

template <typename T, template <typename U> class Descriptor>
void BubbleMatch3D::calculateDisjoiningPressure(
    T pi, MultiScalarField3D<int> &flag, MultiScalarField3D<T> &volumeFraction,
    MultiTensorField3D<T, 3> &normal, MultiScalarField3D<T> &rhoBar,
    MultiScalarField3D<T> &mass, MultiTensorField3D<T, 3> &j,
    bool incompressibleModel, pluint numBubbles) {
  std::vector<MultiBlock3D *> args;
  args.push_back(tagMatrix);
  args.push_back(&flag);
  args.push_back(bubbleAnalysisContainer);
  args.push_back(&volumeFraction);
  args.push_back(&normal);
  args.push_back(&rhoBar);
  args.push_back(&mass);
  args.push_back(&j);

  setToConstant(*disjoiningPressureField,
                disjoiningPressureField->getBoundingBox(), (T)0.);
  args.push_back(disjoiningPressureField);

  applyProcessingFunctional(new CalculateDisjoiningPressure3D<T, Descriptor>(
                                pi, numBubbles, incompressibleModel),
                            bubbleAnalysisContainer->getBoundingBox(), args);
}

template <typename T>
void BubbleMatch3D::bubbleAnalysis(MultiScalarField3D<int> &flag,
                                   MultiScalarField3D<T> &volumeFraction,
                                   pluint numBubbles) {
  std::vector<MultiBlock3D *> args;
  args.push_back(tagMatrix);
  args.push_back(&flag);
  args.push_back(bubbleAnalysisContainer);
  args.push_back(&volumeFraction);
  applyProcessingFunctional(new AnalyzeBubbles3D<T>(numBubbles),
                            bubbleAnalysisContainer->getBoundingBox(), args);
}

/* *************** Class AnalyzeBubbles3D ******************************** */

template <typename T>
AnalyzeBubbles3D<T>::AnalyzeBubbles3D(pluint numBubbles_)
    : numBubbles(numBubbles_) {}

template <typename T>
AnalyzeBubbles3D<T> *AnalyzeBubbles3D<T>::clone() const {
  return new AnalyzeBubbles3D<T>(*this);
}

template <typename T>
void AnalyzeBubbles3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 4);
  ScalarField3D<plint> *pTagMatrix =
      dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField3D<plint> &tagMatrix = *pTagMatrix;

  ScalarField3D<int> *pFlagMatrix =
      dynamic_cast<ScalarField3D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField3D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock3D *pDataBlock =
      dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock3D &dataBlock = *pDataBlock;
  BubbleAnalysisData3D *pData =
      dynamic_cast<BubbleAnalysisData3D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData3D &data = *pData;

  ScalarField3D<T> *pVolumeFraction =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField3D<T> &volumeFraction = *pVolumeFraction;

  Dot3D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot3D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);
  Dot3D absOfs = tagMatrix.getLocation();

  std::vector<double> bubbleVolume(numBubbles);
  std::fill(bubbleVolume.begin(), bubbleVolume.end(), 0.);
  std::vector<Array<double, 3>> bubbleCenter(numBubbles);
  std::fill(bubbleCenter.begin(), bubbleCenter.end(),
            Array<double, 3>(0., 0., 0.));

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
        plint tag = tagMatrix.get(iX, iY, iZ);
        if (tag >= 0) {
          if (flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y,
                             iZ + flagOffset.z) == freeSurfaceFlag3D::empty) {
            PLB_ASSERT(tag < (plint)bubbleVolume.size());
            bubbleVolume[tag] += 1.0;
            bubbleCenter[tag] +=
                Array<double, 3>((double)iX + absOfs.x, (double)iY + absOfs.y,
                                 (double)iZ + absOfs.z);
          } else if (flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y,
                                    iZ + flagOffset.z) ==
                     freeSurfaceFlag3D::interface) {
            PLB_ASSERT(tag < (plint)bubbleVolume.size());
            double vf =
                1.0 - (double)volumeFraction.get(
                          iX + vfOffset.x, iY + vfOffset.y, iZ + vfOffset.z);
            bubbleVolume[tag] += vf;
            bubbleCenter[tag] += vf * Array<double, 3>((double)iX + absOfs.x,
                                                       (double)iY + absOfs.y,
                                                       (double)iZ + absOfs.z);
          } else {
            PLB_ASSERT(false);
          }
        }
      }
    }
  }

  data.bubbleVolume = bubbleVolume;
  data.bubbleCenter = bubbleCenter;
}

/* *************** Class CalculateBubbleGrowth3D
 * ********************************* */
template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth3D<T, AD_Descriptor>::CalculateBubbleGrowth3D(
    pluint numBubbles_)
    : numBubbles(numBubbles_)

{}

template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth3D<T, AD_Descriptor>::CalculateBubbleGrowth3D(
    CalculateBubbleGrowth3D<T, AD_Descriptor> const &rhs)
    : numBubbles(rhs.numBubbles)

{}

template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth3D<T, AD_Descriptor>
    *CalculateBubbleGrowth3D<T, AD_Descriptor>::clone() const {
  return new CalculateBubbleGrowth3D<T, AD_Descriptor>(*this);
}

template <typename T, template <typename U1> class AD_Descriptor>
void CalculateBubbleGrowth3D<T, AD_Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::staticVariables;  // tags.
  modified[1] = modif::nothing;          // flags.
  modified[2] = modif::nothing;          // data.
  modified[3] = modif::nothing;          // volume fraction.
  modified[4] = modif::nothing;          // old volume fraction.
  modified[5] = modif::nothing;          // gas lattice.
}

template <typename T, template <typename U1> class AD_Descriptor>
void CalculateBubbleGrowth3D<T, AD_Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) {
  using namespace freeSurfaceFlag3D;

  PLB_ASSERT(atomicBlocks.size() == 6);
  ScalarField3D<plint> *pTagMatrix =
      dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField3D<plint> &tagMatrix = *pTagMatrix;

  ScalarField3D<int> *pFlagMatrix =
      dynamic_cast<ScalarField3D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField3D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock3D *pDataBlock =
      dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock3D &dataBlock = *pDataBlock;
  BubbleAnalysisData3D *pData =
      dynamic_cast<BubbleAnalysisData3D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData3D &data = *pData;

  ScalarField3D<T> *pVolumeFraction =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField3D<T> &volumeFraction = *pVolumeFraction;

  ScalarField3D<T> *pOldVolumeFraction =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[4]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField3D<T> &oldVolumeFraction = *pOldVolumeFraction;

  BlockLattice3D<T, AD_Descriptor> *pGasLattice =
      dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(atomicBlocks[5]);
  PLB_ASSERT(pGasLattice);
  BlockLattice3D<T, AD_Descriptor> &gasLattice = *pGasLattice;

  std::vector<double> bubbleDensity(numBubbles);
  std::fill(bubbleDensity.begin(), bubbleDensity.end(), 0.);

  Dot3D gasOffset = computeRelativeDisplacement(tagMatrix, gasLattice);
  Dot3D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot3D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY)
      for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
        plint innerFlag = flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y,
                                         iZ + flagOffset.z);
        if (innerFlag == interface) {
          Cell<T, AD_Descriptor> const &gasCell1 = gasLattice.get(
              iX + gasOffset.x, iY + gasOffset.y, iZ + gasOffset.z);

          T incrim = (T)0.;
          for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop)

          {
            plint nextX = iX + AD_Descriptor<T>::c[iPop][0];
            plint nextY = iY + AD_Descriptor<T>::c[iPop][1];
            plint nextZ = iZ + AD_Descriptor<T>::c[iPop][2];

            plint nextFlag =
                flagMatrix.get(nextX + flagOffset.x, nextY + flagOffset.y,
                               nextZ + flagOffset.z);
            if (nextFlag == fluid) {
              Cell<T, AD_Descriptor> const &gasCell2 =
                  gasLattice.get(nextX + gasOffset.x, nextY + gasOffset.y,
                                 nextZ + gasOffset.z);

              plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);

              incrim += gasCell1[opp] - gasCell2[iPop];
            }
          }

          plint tag = tagMatrix.get(iX, iY, iZ);

          PLB_ASSERT(tag >= 0)

          T delta_vof = volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y,
                                           iZ + vfOffset.z) -
                        oldVolumeFraction.get(iX + vfOffset.x, iY + vfOffset.y,
                                              iZ + vfOffset.z);
          bubbleDensity[tag] += incrim - gasCell1.computeDensity() * delta_vof;
        }
      }
  }

  data.bubbleDensity = bubbleDensity;
}

/* *************** Class CalculateDisjoiningPressure3D *
 * ******************************** */
template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure3D<T, Descriptor>::CalculateDisjoiningPressure3D(
    T pi_, pluint numBubbles_, bool incompressibleModel_)
    : pi(pi_),
      numBubbles(numBubbles_),
      incompressibleModel(incompressibleModel_)

{}

template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure3D<T, Descriptor>::CalculateDisjoiningPressure3D(
    CalculateDisjoiningPressure3D<T, Descriptor> const &rhs)
    : pi(rhs.pi),
      numBubbles(rhs.numBubbles),
      incompressibleModel(rhs.incompressibleModel)

{}

template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure3D<T, Descriptor>
    *CalculateDisjoiningPressure3D<T, Descriptor>::clone() const {
  return new CalculateDisjoiningPressure3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CalculateDisjoiningPressure3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::staticVariables;  // tags.
  modified[1] = modif::nothing;          // flags.
  modified[2] = modif::nothing;          // data.
  modified[3] = modif::staticVariables;  // volume fraction.
  modified[4] = modif::nothing;          // normals
  modified[5] = modif::staticVariables;  // rhoBar.
  modified[6] = modif::nothing;          // Mass.
  if (incompressibleModel) {
    modified[7] = modif::nothing;  // j.
  } else {
    modified[7] = modif::staticVariables;  // j.
  }
  modified[8] = modif::staticVariables;  // disjoining pressure.
}

template <typename T, template <typename U> class Descriptor>
void CalculateDisjoiningPressure3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks) {
  using namespace freeSurfaceFlag3D;

  PLB_ASSERT(atomicBlocks.size() == 9);
  ScalarField3D<plint> *pTagMatrix =
      dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField3D<plint> &tagMatrix = *pTagMatrix;

  ScalarField3D<int> *pFlagMatrix =
      dynamic_cast<ScalarField3D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField3D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock3D *pDataBlock =
      dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock3D &dataBlock = *pDataBlock;
  BubbleAnalysisData3D *pData =
      dynamic_cast<BubbleAnalysisData3D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData3D &data = *pData;

  ScalarField3D<T> *pVolumeFraction =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField3D<T> &volumeFraction = *pVolumeFraction;

  TensorField3D<T, 3> *pNormal =
      dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[4]);
  PLB_ASSERT(pNormal);
  TensorField3D<T, 3> &normal = *pNormal;

  ScalarField3D<T> *pDensity =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[5]);
  PLB_ASSERT(pDensity);
  ScalarField3D<T> &density = *pDensity;

  ScalarField3D<T> *pMass = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[6]);
  PLB_ASSERT(pMass);
  ScalarField3D<T> &mass = *pMass;

  TensorField3D<T, 3> *pJ =
      dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[7]);
  PLB_ASSERT(pJ);
  TensorField3D<T, 3> &j = *pJ;

  ScalarField3D<T> *pDisjoiningPressureField =
      dynamic_cast<ScalarField3D<T> *>(atomicBlocks[8]);
  ScalarField3D<T> &disjoiningPressureField = *pDisjoiningPressureField;
  PLB_ASSERT(pDisjoiningPressureField);

  std::vector<double> bubbleDisjoiningPressure(numBubbles);
  std::fill(bubbleDisjoiningPressure.begin(), bubbleDisjoiningPressure.end(),
            0.);

  Dot3D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot3D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);
  Dot3D normalOffset = computeRelativeDisplacement(tagMatrix, normal);
  Dot3D massOffset = computeRelativeDisplacement(tagMatrix, mass);
  Dot3D jOffset = computeRelativeDisplacement(tagMatrix, j);
  Dot3D densityOffset = computeRelativeDisplacement(tagMatrix, density);
  Dot3D disjoiningOffset =
      computeRelativeDisplacement(tagMatrix, disjoiningPressureField);

  typedef typename InterfaceLists3D<T, Descriptor>::Node Node;
  typedef Descriptor<T> D;

  std::map<Node, T> densitiesMap;

  std::vector<std::vector<T>> vertp;

  Array<T, 3> iNormal;
  // outward normal (away from the bubble)
  Array<T, 3> outwardiNormal;

  T halfLength = (T)0.5;
  T fullLength = (T)1.;

  T tMaxX;
  T tMaxY;
  T tMaxZ;
  T tDeltaX, tDeltaY, tDeltaZ;

  plint stepX, stepY, stepZ;

  // Disjoining pressure effective distance in terms of lattice units.
  int maxLim = 4;
  plint nextX, nextY, nextZ;

  plint iFlag, jFlag;

  plint iBubbleTag;

  T iVF;
  //    T jVF;
  T iDelta = 0.;
  //    T jDelta = 0.;

  T c, dx;
  dx = 1.;

  int no_of_rows = 120;
  int no_of_cols = 3;
  T initial_value = 0.;

  T distance = T();
  T densityChange = T();

  Array<int, 3> integrationVector;
  Array<T, 3> dAndDeltaStar;

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY)
      for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
        iFlag = flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y,
                               iZ + flagOffset.z);
        if (iFlag == interface) {
          iBubbleTag = tagMatrix.get(iX, iY, iZ);

          // outward normal (away from the bubble)
          iNormal = normal.get(iX + normalOffset.x, iY + normalOffset.y,
                               iZ + normalOffset.z);
          outwardiNormal = -1. * iNormal;

          tMaxX = std::abs(halfLength / outwardiNormal[0]);
          tMaxY = std::abs(halfLength / outwardiNormal[1]);
          tMaxZ = std::abs(halfLength / outwardiNormal[2]);
          tDeltaX = std::abs(fullLength / outwardiNormal[0]);
          tDeltaY = std::abs(fullLength / outwardiNormal[1]);
          tDeltaZ = std::abs(fullLength / outwardiNormal[2]);

          stepX = outwardiNormal[0] >= 0.0 ? 1 : -1;
          stepY = outwardiNormal[1] >= 0.0 ? 1 : -1;
          stepZ = outwardiNormal[2] >= 0.0 ? 1 : -1;

          // Disjoining pressure effective distance in terms of lattice units.
          nextX = iX;
          nextY = iY;
          nextZ = iZ;

          // Traversal algorithm
          for (int d = 1; d <= maxLim; d++) {
            if ((tMaxX < tMaxY || std::isinf(tMaxY)) &&
                tMaxX != -std::numeric_limits<T>::infinity()) {
              if ((tMaxX < tMaxZ || std::isinf(tMaxZ)) &&
                  tMaxX != -std::numeric_limits<T>::infinity()) {
                tMaxX += tDeltaX;
                nextX += stepX;
              } else {
                nextZ += stepZ;
                tMaxZ += tDeltaZ;
              }
            } else {
              if ((tMaxY < tMaxZ || std::isinf(tMaxZ)) &&
                  tMaxY != -std::numeric_limits<T>::infinity()) {
                nextY += stepY;
                tMaxY += tDeltaY;
              } else {
                nextZ += stepZ;
                tMaxZ = tMaxZ + tDeltaZ;
              }
            }

            // Skip iBubble interface cells
            if (tagMatrix.get(nextX, nextY, nextZ) == iBubbleTag) {
              continue;
            }

            jFlag = flagMatrix.get(nextX + flagOffset.x, nextY + flagOffset.y,
                                   nextZ + flagOffset.z);

            if (jFlag == wall || jFlag == slipWall) {
              break;
            } else if (jFlag == interface) {
              T abs0 = std::fabs(iNormal[0]);
              T abs1 = std::fabs(iNormal[1]);
              T abs2 = std::fabs(iNormal[2]);

              int integrationDirection = 2;

              if (abs0 > abs1) {
                if (abs0 > abs2) {
                  integrationDirection = 0;
                } else {
                  integrationDirection = 2;
                }
              } else {
                if (abs1 > abs2) {
                  integrationDirection = 1;
                }
                // abs0 <= abs1 && abs1 <= abs2
                else {
                  integrationDirection = 2;
                }
              }

              // Compute the vector parallel to the integration direction.
              integrationVector[0] = integrationDirection == 0 ? 1 : 0;
              integrationVector[1] = integrationDirection == 1 ? 1 : 0;
              integrationVector[2] = integrationDirection == 2 ? 1 : 0;

              iVF = volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y,
                                       iZ + vfOffset.z);

              // jVF = volumeFraction.get(nextX + vfOffset.x, nextY +
              // vfOffset.y, nextZ + vfOffset.z);

              // jNormal = normal.get(nextX + normalOffset.x, nextY +
              // normalOffset.y, nextZ + normalOffset.z);
              //   outwardjNormal = -1. * jNormal;

              // Calculation of delta

              vertp.clear();
              vertp.resize(no_of_rows,
                           std::vector<T>(no_of_cols, initial_value));

              cubicmesh(vertp);
              enforv3dsz(c, dx, dx, dx, iVF, vertp, outwardiNormal);

              iDelta =
                  (outwardiNormal[0] + outwardiNormal[1] + outwardiNormal[2]) *
                      0.5 +
                  c;

              dAndDeltaStar[0] =
                  (std::abs(nextX - iX) + 0.5 -
                   (1. - volumeFraction.get(nextX + vfOffset.x,
                                            nextY + vfOffset.y,
                                            nextZ + vfOffset.z))) *
                  integrationVector[0];
              dAndDeltaStar[1] =
                  (std::abs(nextY - iY) + 0.5 -
                   (1. - volumeFraction.get(nextX + vfOffset.x,
                                            nextY + vfOffset.y,
                                            nextZ + vfOffset.z))) *
                  integrationVector[1];
              dAndDeltaStar[2] =
                  (std::abs(nextZ - iZ) + 0.5 -
                   (1. - volumeFraction.get(nextX + vfOffset.x,
                                            nextY + vfOffset.y,
                                            nextZ + vfOffset.z))) *
                  integrationVector[2];

              T dAndDelta = std::abs(dot(dAndDeltaStar, outwardiNormal));

              distance = dAndDelta - iDelta;

              distance = distance < maxLim ? distance : maxLim;

              bubbleDisjoiningPressure[iBubbleTag] += pi * (maxLim - distance);

              disjoiningPressureField.get(iX + disjoiningOffset.x,
                                          iY + disjoiningOffset.y,
                                          iZ + disjoiningOffset.z) =
                  pi * (maxLim - distance) * D::invCs2;

              densityChange = pi * (maxLim - distance) * D::invCs2;

              // densitiesMap.insert(std::pair<Node,T> (Node(iX, iY, iZ),
              // densityChange));
              densitiesMap[Node(iX, iY, iZ)] += densityChange;
            }
          }
        }
      }
  }

  typename std::map<Node, T>::const_iterator nodes = densitiesMap.begin();

  Array<T, 3> iniJ, newJ;

  for (; nodes != densitiesMap.end(); ++nodes) {
    Node node = nodes->first;
    plint iX = node[0];
    plint iY = node[1];
    plint iZ = node[2];
    T densityChange = nodes->second;
    T iniDensity = Descriptor<T>::fullRho(density.get(
        iX + densityOffset.x, iY + densityOffset.y, iZ + densityOffset.z));
    T newDensity = iniDensity - densityChange;
    density.get(iX + densityOffset.x, iY + densityOffset.y,
                iZ + densityOffset.z) = Descriptor<T>::rhoBar(newDensity);
    volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y, iZ + vfOffset.z) =
        mass.get(iX + massOffset.x, iY + massOffset.y, iZ + massOffset.z) /
        newDensity;
    // On interface cells, adjust the pressure to incorporate disjoining
    // pressure.
    if (!incompressibleModel) {
      iniJ = j.get(iX + jOffset.x, iY + jOffset.y, iZ + jOffset.z);
      newJ = iniJ * newDensity / iniDensity;
      j.get(iX + jOffset.x, iY + jOffset.y, iZ + jOffset.z) = newJ;
    }
  }

  data.bubbleDisjoiningPressure = bubbleDisjoiningPressure;
}

template <typename T, template <typename U> class Descriptor>
void BubbleMatch3D::stabilizeLattice(
    FreeSurfaceFields3D<T, Descriptor> &fields) {
  applyProcessingFunctional(new FreeSurfaceStabilize3D<T, Descriptor>(),
                            fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(fields.rhoDefault),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
          fields.dynamics->clone(), fields.force),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(
          fields.rhoDefault),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(new FreeSurfaceComputeStatistics3D<T, Descriptor>,
                            fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);
}

}  // namespace plb

#endif  // BUBBLE_MATCH_3D_HH
