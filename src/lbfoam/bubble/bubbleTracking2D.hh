/*
 *  LBfoam: An open-source software package for the simulation of foaming
 *  using the Lattice Boltzmann Method
 *  Copyright (C) 2020 Mohammadmehdi Ataei
 *  m.ataei@mail.utoronto.ca
 *  This file is part of LBfoam.
 *
 *  LBfoam is free software: you can redistribute it and/or modify it under
 *  the terms of the GNU Affero General Public License as published by the
 *  Free Software Foundation version 3.
 *
 *  LBfoam is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this Program. If not, see <http://www.gnu.org/licenses/>.
 *
 *  #############################################################################
 *
 *  Author:         Mohammadmehdi Ataei, 2020
 *
 *  #############################################################################
 *
 *  Parts of the LBfoam code that originate from Palabos are distributed
 *  under the terms of the AGPL 3.0 license with the following copyright
 *  notice:
 *
 *  This file is part of the Palabos library.
 *
 *  Copyright (C) 2011-2017 FlowKit Sarl
 *  Route d'Oron 2
 *  1010 Lausanne, Switzerland
 *  E-mail contact: contact@flowkit.com
 *
 *  The most recent release of Palabos can be downloaded at
 *  <http://www.palabos.org/>
 *
 *  The library Palabos is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  The library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef BUBBLE_TRACKING_2D_HH
#define BUBBLE_TRACKING_2D_HH

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "lbfoam/bubble/bubbleTracking2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"

namespace plb {
namespace lbfoam {

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void BubbleTracking2D::execute(
    MultiScalarField2D<T> &volumeFraction, MultiScalarField2D<int> &flag,
    MultiTensorField2D<T, 2> &normal, MultiScalarField2D<T> &rhoBar,
    MultiScalarField2D<T> &mass, MultiTensorField2D<T, 2> &j,
    MultiBlockLattice2D<T, ScalarDescriptor> &adLattice,
    MultiScalarField2D<T> &oldVolumeFraction, bool incompressibleModel, T pi) {
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
void BubbleTracking2D::calculateBubbleGrowth(
    MultiScalarField2D<int> &flag, MultiScalarField2D<T> &volumeFraction,
    MultiScalarField2D<T> &oldVolumeFraction,
    MultiBlockLattice2D<T, AD_Descriptor> &adLattice, pluint numBubbles)

{
  std::vector<MultiBlock2D *> args;
  args.push_back(tagMatrix);
  args.push_back(&flag);
  args.push_back(bubbleAnalysisContainer);
  args.push_back(&volumeFraction);
  args.push_back(&oldVolumeFraction);
  args.push_back(&adLattice);
  applyProcessingFunctional(
      new CalculateBubbleGrowth2D<T, AD_Descriptor>(numBubbles),
      bubbleAnalysisContainer->getBoundingBox(), args);
}

template <typename T, template <typename U> class Descriptor>
void BubbleTracking2D::calculateDisjoiningPressure(
    T pi, MultiScalarField2D<int> &flag, MultiScalarField2D<T> &volumeFraction,
    MultiTensorField2D<T, 2> &normal, MultiScalarField2D<T> &rhoBar,
    MultiScalarField2D<T> &mass, MultiTensorField2D<T, 2> &j,
    bool incompressibleModel, pluint numBubbles) {
  std::vector<MultiBlock2D *> args;

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

  applyProcessingFunctional(new CalculateDisjoiningPressure2D<T, Descriptor>(
                                pi, numBubbles, incompressibleModel),
                            bubbleAnalysisContainer->getBoundingBox(), args);
}

template <typename T>
void BubbleTracking2D::bubbleAnalysis(MultiScalarField2D<int> &flag,
                                      MultiScalarField2D<T> &volumeFraction,
                                      pluint numBubbles) {
  std::vector<MultiBlock2D *> args;
  args.push_back(tagMatrix);
  args.push_back(&flag);
  args.push_back(bubbleAnalysisContainer);
  args.push_back(&volumeFraction);
  applyProcessingFunctional(new AnalyzeBubbles2D<T>(numBubbles),
                            bubbleAnalysisContainer->getBoundingBox(), args);
}

/* *************** Class AnalyzeBubbles2D ******************************** */

template <typename T>
AnalyzeBubbles2D<T>::AnalyzeBubbles2D(pluint numBubbles_)
    : numBubbles(numBubbles_) {}

template <typename T>
AnalyzeBubbles2D<T> *AnalyzeBubbles2D<T>::clone() const {
  return new AnalyzeBubbles2D<T>(*this);
}

template <typename T>
void AnalyzeBubbles2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 4);
  ScalarField2D<plint> *pTagMatrix =
      dynamic_cast<ScalarField2D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField2D<plint> &tagMatrix = *pTagMatrix;

  ScalarField2D<int> *pFlagMatrix =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField2D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock2D *pDataBlock =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock2D &dataBlock = *pDataBlock;
  BubbleAnalysisData2D *pData =
      dynamic_cast<BubbleAnalysisData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData2D &data = *pData;

  ScalarField2D<T> *pVolumeFraction =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField2D<T> &volumeFraction = *pVolumeFraction;

  Dot2D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot2D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);
  Dot2D absOfs = tagMatrix.getLocation();

  std::vector<double> bubbleVolume(numBubbles);
  std::fill(bubbleVolume.begin(), bubbleVolume.end(), 0.);
  std::vector<Array<double, 2>> bubbleCenter(numBubbles);
  std::fill(bubbleCenter.begin(), bubbleCenter.end(), Array<double, 2>(0., 0.));

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint tag = tagMatrix.get(iX, iY);
      if (tag >= 0) {
        if (flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y) ==
            freeSurfaceFlag2D::empty) {
          PLB_ASSERT(tag < (plint)bubbleVolume.size());
          bubbleVolume[tag] += 1.0;
          bubbleCenter[tag] +=
              Array<double, 2>((double)iX + absOfs.x, (double)iY + absOfs.y);
        } else if (flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y) ==
                   freeSurfaceFlag2D::interface) {
          PLB_ASSERT(tag < (plint)bubbleVolume.size());
          double vf = 1.0 - (double)volumeFraction.get(iX + vfOffset.x,
                                                       iY + vfOffset.y);
          bubbleVolume[tag] += vf;
          bubbleCenter[tag] += vf * Array<double, 2>((double)iX + absOfs.x,
                                                     (double)iY + absOfs.y);
        } else {
          PLB_ASSERT(false);
        }
      }
    }
  }
  data.bubbleVolume = bubbleVolume;
  data.bubbleCenter = bubbleCenter;
}

/* *************** Class CalculateBubbleGrowth2D
 * ********************************* */
template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth2D<T, AD_Descriptor>::CalculateBubbleGrowth2D(
    pluint numBubbles_)
    : numBubbles(numBubbles_)

{}

template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth2D<T, AD_Descriptor>::CalculateBubbleGrowth2D(
    CalculateBubbleGrowth2D<T, AD_Descriptor> const &rhs)
    : numBubbles(rhs.numBubbles)

{}

template <typename T, template <typename U1> class AD_Descriptor>
CalculateBubbleGrowth2D<T, AD_Descriptor>
    *CalculateBubbleGrowth2D<T, AD_Descriptor>::clone() const {
  return new CalculateBubbleGrowth2D<T, AD_Descriptor>(*this);
}

template <typename T, template <typename U1> class AD_Descriptor>
void CalculateBubbleGrowth2D<T, AD_Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::staticVariables;  // tags.
  modified[1] = modif::nothing;          // flags.
  modified[2] = modif::nothing;          // data.
  modified[3] = modif::nothing;          // volume fraction.
  modified[4] = modif::nothing;          // old volume fraction.
  modified[5] = modif::nothing;          // gas lattice.
}

template <typename T, template <typename U1> class AD_Descriptor>
void CalculateBubbleGrowth2D<T, AD_Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  PLB_ASSERT(atomicBlocks.size() == 6);
  ScalarField2D<plint> *pTagMatrix =
      dynamic_cast<ScalarField2D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField2D<plint> &tagMatrix = *pTagMatrix;

  ScalarField2D<int> *pFlagMatrix =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField2D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock2D *pDataBlock =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock2D &dataBlock = *pDataBlock;
  BubbleAnalysisData2D *pData =
      dynamic_cast<BubbleAnalysisData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData2D &data = *pData;

  ScalarField2D<T> *pVolumeFraction =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField2D<T> &volumeFraction = *pVolumeFraction;

  ScalarField2D<T> *pOldVolumeFraction =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[4]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField2D<T> &oldVolumeFraction = *pOldVolumeFraction;

  BlockLattice2D<T, AD_Descriptor> *pGasLattice =
      dynamic_cast<BlockLattice2D<T, AD_Descriptor> *>(atomicBlocks[5]);
  PLB_ASSERT(pGasLattice);
  BlockLattice2D<T, AD_Descriptor> &gasLattice = *pGasLattice;

  std::vector<double> bubbleDensity(numBubbles);
  std::fill(bubbleDensity.begin(), bubbleDensity.end(), 0.);

  Dot2D gasOffset = computeRelativeDisplacement(tagMatrix, gasLattice);
  Dot2D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot2D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint innerFlag = flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y);
      if (innerFlag == interface) {
        Cell<T, AD_Descriptor> const &gasCell1 =
            gasLattice.get(iX + gasOffset.x, iY + gasOffset.y);

        T incrim = (T)0.;
        for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop)

        {
          plint nextX = iX + AD_Descriptor<T>::c[iPop][0];
          plint nextY = iY + AD_Descriptor<T>::c[iPop][1];

          plint nextFlag =
              flagMatrix.get(nextX + flagOffset.x, nextY + flagOffset.y);
          if (nextFlag == fluid) {
            Cell<T, AD_Descriptor> const &gasCell2 =
                gasLattice.get(nextX + gasOffset.x, nextY + gasOffset.y);

            plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);

            incrim += gasCell1[opp] - gasCell2[iPop];

            //                        pcout << "  cellinfluid: " <<
            //                        gasCell2.computeDensity() << "
            //                        cellininterface:   " <<
            //                        gasCell1.computeDensity()   << std::endl;
          }
        }

        plint tag = tagMatrix.get(iX, iY);

        PLB_ASSERT(tag >= 0)

        T delta_vof = volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y) -
                      oldVolumeFraction.get(iX + vfOffset.x, iY + vfOffset.y);
        //                bubbleDensity[tag] += incrim -
        //                gasCell1.computeDensity() * delta_vof;
        bubbleDensity[tag] += (incrim - gasCell1.computeDensity() * delta_vof);

        //                // testing here
        //                bubbleDensity[tag] = 0;

        //                pcout << "  incrim: " << incrim << "  increase: " <<
        //                incrim - gasCell1.computeDensity() * delta_vof   <<
        //                std::endl;
      }
    }
  }

  data.bubbleDensity = bubbleDensity;
}

/* *************** Class CalculateDisjoiningPressure2D *
 * ******************************** */
template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure2D<T, Descriptor>::CalculateDisjoiningPressure2D(
    T pi_, pluint numBubbles_, bool incompressibleModel_)
    : pi(pi_),
      numBubbles(numBubbles_),
      incompressibleModel(incompressibleModel_)

{}

template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure2D<T, Descriptor>::CalculateDisjoiningPressure2D(
    CalculateDisjoiningPressure2D<T, Descriptor> const &rhs)
    : pi(rhs.pi),
      numBubbles(rhs.numBubbles),
      incompressibleModel(rhs.incompressibleModel)

{}

template <typename T, template <typename U> class Descriptor>
CalculateDisjoiningPressure2D<T, Descriptor>
    *CalculateDisjoiningPressure2D<T, Descriptor>::clone() const {
  return new CalculateDisjoiningPressure2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CalculateDisjoiningPressure2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // tags.
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
void CalculateDisjoiningPressure2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  PLB_ASSERT(atomicBlocks.size() == 9);
  ScalarField2D<plint> *pTagMatrix =
      dynamic_cast<ScalarField2D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField2D<plint> &tagMatrix = *pTagMatrix;

  ScalarField2D<int> *pFlagMatrix =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(pFlagMatrix);
  ScalarField2D<int> &flagMatrix = *pFlagMatrix;

  AtomicContainerBlock2D *pDataBlock =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[2]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock2D &dataBlock = *pDataBlock;
  BubbleAnalysisData2D *pData =
      dynamic_cast<BubbleAnalysisData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleAnalysisData2D &data = *pData;

  ScalarField2D<T> *pVolumeFraction =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
  PLB_ASSERT(pVolumeFraction);
  ScalarField2D<T> &volumeFraction = *pVolumeFraction;

  TensorField2D<T, 2> *pNormal =
      dynamic_cast<TensorField2D<T, 2> *>(atomicBlocks[4]);
  PLB_ASSERT(pNormal);
  TensorField2D<T, 2> &normal = *pNormal;

  ScalarField2D<T> *pDensity =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[5]);
  PLB_ASSERT(pDensity);
  ScalarField2D<T> &density = *pDensity;

  ScalarField2D<T> *pMass = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[6]);
  PLB_ASSERT(pMass);
  ScalarField2D<T> &mass = *pMass;

  TensorField2D<T, 2> *pJ =
      dynamic_cast<TensorField2D<T, 2> *>(atomicBlocks[7]);
  PLB_ASSERT(pJ);
  TensorField2D<T, 2> &j = *pJ;

  ScalarField2D<T> *pDisjoiningPressureField =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[8]);
  ScalarField2D<T> &disjoiningPressureField = *pDisjoiningPressureField;
  PLB_ASSERT(pDisjoiningPressureField);

  std::vector<double> bubbleDisjoiningPressure(numBubbles);
  std::fill(bubbleDisjoiningPressure.begin(), bubbleDisjoiningPressure.end(),
            0.);

  Dot2D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  Dot2D vfOffset = computeRelativeDisplacement(tagMatrix, volumeFraction);
  Dot2D normalOffset = computeRelativeDisplacement(tagMatrix, normal);
  Dot2D massOffset = computeRelativeDisplacement(tagMatrix, mass);
  Dot2D jOffset = computeRelativeDisplacement(tagMatrix, j);
  Dot2D densityOffset = computeRelativeDisplacement(tagMatrix, density);
  Dot2D disjoiningOffset =
      computeRelativeDisplacement(tagMatrix, disjoiningPressureField);

  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  typedef Descriptor<T> D;

  std::map<Node, T> densitiesMap;

  std::vector<std::vector<T>> vertp;

  Array<T, 2> iNormal;
  // outward normal (away from the bubble)
  Array<T, 2> outwardiNormal;

  // Disjoining pressure effective distance in terms of lattice units.
  int maxLim = 4;
  plint nextX, nextY;

  plint iFlag, jFlag;

  plint iBubbleTag;

  T iVF;
  T iDelta = 0.;
  // T jDelta = 0.;

  T c, dx;
  dx = 1.;
  T dAndDelta;

  int no_of_rows = 120;
  int no_of_cols = 3;
  T initial_value = 0.;
  T distance = T();
  T densityChange = T();
  Array<int, 2> marchingVector;

  Array<T, 2> iniJ;
  Array<T, 2> newJ;

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      iFlag = flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y);
      if (iFlag == interface) {
        iBubbleTag = tagMatrix.get(iX, iY);

        // outward normal
        iNormal = normal.get(iX + normalOffset.x, iY + normalOffset.y);
        outwardiNormal = -1. * iNormal;

        rayTracer2D<T> rayTracer = rayTracer2D<T>(iX, iY, outwardiNormal);
        for (int d = 1; d <= maxLim; d++) {
          Array<T, 2> nextCell = rayTracer.findNextCell();
          nextX = nextCell[0];
          nextY = nextCell[1];

          // Skip iBubble interface cells
          if (tagMatrix.get(nextX, nextY) == iBubbleTag) {
            continue;
          }

          jFlag = flagMatrix.get(nextX + flagOffset.x, nextY + flagOffset.y);

          if (jFlag == wall || jFlag == slipWall) {
            break;
          } else if (jFlag == interface) {
            T abs0 = std::fabs(iNormal[0]);
            T abs1 = std::fabs(iNormal[1]);

            iVF = volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y);

            int marchingDirection = 1;

            if (abs0 > abs1) {
              marchingDirection = 0;
            }

            // Compute the vector parallel to the integration direction.
            marchingVector[0] = marchingDirection == 0 ? 1 : 0;
            marchingVector[1] = marchingDirection == 1 ? 1 : 0;

            // Calculation of delta

            vertp.clear();
            vertp.resize(no_of_rows, std::vector<T>(no_of_cols, initial_value));

            squareMesh(vertp);
            enforv2dsz(c, dx, dx, iVF, vertp, outwardiNormal);

            iDelta = (outwardiNormal[0] + outwardiNormal[1]) * 0.5 + c;

            T dAndDeltaStar[2];
            dAndDeltaStar[0] = (std::abs(nextX - iX) + 0.5 -
                                (1. - volumeFraction.get(nextX + vfOffset.x,
                                                         nextY + vfOffset.y))) *
                               marchingVector[0];
            dAndDeltaStar[1] = (std::abs(nextY - iY) + 0.5 -
                                (1. - volumeFraction.get(nextX + vfOffset.x,
                                                         nextY + vfOffset.y))) *
                               marchingVector[1];

            dAndDelta = std::abs(dot(dAndDeltaStar, outwardiNormal));

            distance = dAndDelta - iDelta;

            distance = distance < maxLim ? std::abs(distance) : maxLim;

            bubbleDisjoiningPressure[iBubbleTag] += pi * (maxLim - distance);

            disjoiningPressureField.get(iX + disjoiningOffset.x,
                                        iY + disjoiningOffset.y) =
                pi * std::abs(maxLim - distance) * D::invCs2;

            densityChange = pi * (maxLim - distance) * D::invCs2;

            PLB_ASSERT(std::isfinite(densityChange))

            densitiesMap[Node(iX, iY)] += densityChange;
          }
        }
      }
    }
  }

  typename std::map<Node, T>::const_iterator nodes = densitiesMap.begin();
  for (; nodes != densitiesMap.end(); ++nodes) {
    Node node = nodes->first;
    plint iX = node[0];
    plint iY = node[1];
    T densityChange = nodes->second;
    T iniDensity = Descriptor<T>::fullRho(
        density.get(iX + densityOffset.x, iY + densityOffset.y));
    T newDensity = iniDensity - densityChange;
    density.get(iX + densityOffset.x, iY + densityOffset.y) =
        Descriptor<T>::rhoBar(newDensity);
    volumeFraction.get(iX + vfOffset.x, iY + vfOffset.y) =
        mass.get(iX + massOffset.x, iY + massOffset.y) / newDensity;
    // On interface cells, adjust the pressure to incorporate disjoining
    // pressure.
    if (!incompressibleModel) {
      iniJ = j.get(iX + jOffset.x, iY + jOffset.y);
      newJ = iniJ * newDensity / iniDensity;
      j.get(iX + jOffset.x, iY + jOffset.y) = newJ;
    }
  }

  data.bubbleDisjoiningPressure = bubbleDisjoiningPressure;
}

}  // namespace lbfoam
}  // namespace plb

#endif  // BUBBLE_TRACKING_2D_HH
