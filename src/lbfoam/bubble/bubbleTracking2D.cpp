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

#include "lbfoam/bubble/bubbleTracking2D.h"

#include <limits>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField2D.hh"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional2D.hh"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "dataProcessors/dataInitializerFunctional2D.hh"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "dataProcessors/dataInitializerWrapper2D.hh"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiBlockGenerator2D.hh"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataField2D.hh"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.hh"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/nonLocalTransfer2D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/serialMultiDataField2D.h"
#include "multiBlock/serialMultiDataField2D.hh"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"
#include "parallelism/parallelMultiDataField2D.h"
#include "parallelism/parallelMultiDataField2D.hh"

namespace plb {

namespace lbfoam {

/* ************** class BubbleTracking2D ********************************** */

BubbleTracking2D::BubbleTracking2D(MultiBlock2D &templ)
    : bubbleContainer(
          createContainerBlock(templ, new BubbleCounterData2D(maxNumBubbles))),
      bubbleAnalysisContainer(
          createContainerBlock(templ, new BubbleAnalysisData2D())),
      bubbleRemapContainer(
          createContainerBlock(templ, new BubbleRemapData2D(maxNumBubbles))),
      mpiData(*bubbleContainer),
      tagMatrix(new MultiScalarField2D<plint>(*bubbleContainer)),
      disjoiningPressureField(
          new MultiScalarField2D<double>(*bubbleContainer)) {
  setToConstant(*tagMatrix, tagMatrix->getBoundingBox(), (plint)-1);
  tagMatrix->setMultiBlockManagement().changeEnvelopeWidth((plint)4);
}

BubbleTracking2D::~BubbleTracking2D() {
  delete bubbleContainer;
  delete bubbleAnalysisContainer;
  delete bubbleRemapContainer;
  delete tagMatrix;
  delete disjoiningPressureField;
}

pluint BubbleTracking2D::countAndTagBubbles() {
  std::vector<MultiBlock2D *> args;
  args.push_back(tagMatrix);
  args.push_back(bubbleRemapContainer);
  applyProcessingFunctional(new CollectBubbleTags2D(),
                            bubbleRemapContainer->getBoundingBox(), args);
  plint numBubbles = globalBubbleIds();
  applyProcessingFunctional(new ApplyTagRemap2D(),
                            bubbleRemapContainer->getBoundingBox(), args);
  return numBubbles;
}

void BubbleTracking2D::computeBubbleData(pluint numBubbles) {
  std::vector<double> bubbleCenterX(numBubbles), bubbleCenterY(numBubbles);
  bubbleVolume.resize(numBubbles);
  bubbleCenter.resize(numBubbles);
  bubbleDensity.resize(numBubbles);
  bubbleDisjoiningPressure.resize(numBubbles);

  std::vector<plint> const &localIds = mpiData.getLocalIds();
  for (pluint i = 0; i < localIds.size(); ++i) {
    plint id = localIds[i];
    AtomicContainerBlock2D &atomicDataContainer =
        bubbleAnalysisContainer->getComponent(id);
    BubbleAnalysisData2D *pData =
        dynamic_cast<BubbleAnalysisData2D *>(atomicDataContainer.getData());
    PLB_ASSERT(pData);

    BubbleAnalysisData2D &data = *pData;
    std::vector<double> const &nextVolume = data.bubbleVolume;
    std::vector<double> const &nextDensity = data.bubbleDensity;
    std::vector<double> const &nextDisjoiningPressure =
        data.bubbleDisjoiningPressure;
    std::vector<Array<double, 2>> const &nextCenter = data.bubbleCenter;
    PLB_ASSERT(nextVolume.size() == numBubbles);
    PLB_ASSERT(nextCenter.size() == numBubbles);

    for (pluint i = 0; i < numBubbles; ++i) {
      bubbleVolume[i] += nextVolume[i];
      bubbleDensity[i] += nextDensity[i];
      bubbleDisjoiningPressure[i] += nextDisjoiningPressure[i];
      bubbleCenterX[i] += nextCenter[i][0];
      bubbleCenterY[i] += nextCenter[i][1];
    }
  }
#ifdef PLB_MPI_PARALLEL
  global::mpi().allReduceVect(bubbleVolume, MPI_SUM);
  global::mpi().allReduceVect(bubbleDensity, MPI_SUM);
  global::mpi().allReduceVect(bubbleDisjoiningPressure, MPI_SUM);
  global::mpi().allReduceVect(bubbleCenterX, MPI_SUM);
  global::mpi().allReduceVect(bubbleCenterY, MPI_SUM);
#endif
  static const double epsilon = std::numeric_limits<double>::epsilon() * 1.e4;
  for (pluint i = 0; i < numBubbles; ++i) {
    bubbleCenter[i] = Array<double, 2>(bubbleCenterX[i], bubbleCenterY[i]);
    double volume = bubbleVolume[i];
    if (volume > epsilon) {
      bubbleCenter[i] /= volume;
    }
  }
}

plint BubbleTracking2D::globalBubbleIds() {
  plint localNumUniqueBubbles = 0;
  std::vector<plint> const &localIds = mpiData.getLocalIds();
  for (pluint i = 0; i < localIds.size(); ++i) {
    plint id = localIds[i];
    AtomicContainerBlock2D &atomicDataContainer =
        bubbleRemapContainer->getComponent(id);
    BubbleRemapData2D *pData =
        dynamic_cast<BubbleRemapData2D *>(atomicDataContainer.getData());
    PLB_ASSERT(pData);
    BubbleRemapData2D &data = *pData;
    localNumUniqueBubbles += data.getUniqueTags().size();
  }

  std::vector<plint> allNumBubbles(global::mpi().getSize());
  allNumBubbles[global::mpi().getRank()] = localNumUniqueBubbles;
#ifdef PLB_MPI_PARALLEL
  global::mpi().allReduceVect(allNumBubbles, MPI_SUM);
#endif

  std::vector<plint> cumNumBubbles(global::mpi().getSize());
  PLB_ASSERT(cumNumBubbles.size() > 0);
  cumNumBubbles[0] = allNumBubbles[0];
  for (pluint i = 1; i < cumNumBubbles.size(); ++i) {
    cumNumBubbles[i] = allNumBubbles[i] + cumNumBubbles[i - 1];
  }
  plint totNumBubbles = cumNumBubbles.back();

  std::vector<plint> bubbleIds(totNumBubbles);
  pluint offset = global::mpi().getRank() == 0
                      ? 0
                      : cumNumBubbles[global::mpi().getRank() - 1];

  for (pluint i = 0; i < localIds.size(); ++i) {
    plint id = localIds[i];
    AtomicContainerBlock2D &atomicDataContainer =
        bubbleRemapContainer->getComponent(id);
    BubbleRemapData2D *pData =
        dynamic_cast<BubbleRemapData2D *>(atomicDataContainer.getData());
    PLB_ASSERT(pData);
    BubbleRemapData2D &data = *pData;
    std::vector<plint> uniqueTags = data.getUniqueTags();
    for (pluint i = 0; i < uniqueTags.size(); ++i, ++offset) {
      bubbleIds[offset] = uniqueTags[i];
    }
  }
#ifdef PLB_MPI_PARALLEL
  global::mpi().allReduceVect(bubbleIds, MPI_SUM);
#endif
  std::map<plint, plint> tagRemap;
  for (pluint i = 0; i < bubbleIds.size(); ++i) {
    tagRemap[bubbleIds[i]] = i;
  }

  for (pluint i = 0; i < localIds.size(); ++i) {
    plint id = localIds[i];
    AtomicContainerBlock2D &atomicDataContainer =
        bubbleRemapContainer->getComponent(id);
    BubbleRemapData2D *pData =
        dynamic_cast<BubbleRemapData2D *>(atomicDataContainer.getData());
    PLB_ASSERT(pData);
    BubbleRemapData2D &data = *pData;
    data.getTagRemap() = tagRemap;
  }

  return totNumBubbles;
}

void BubbleTracking2D::resetBubbleContainer() {
  MultiBlockManagement2D const &management =
      bubbleContainer->getMultiBlockManagement();
  ThreadAttribution const &threadAttribution =
      management.getThreadAttribution();
  SparseBlockStructure2D const &sparseBlock =
      management.getSparseBlockStructure();
  std::map<plint, Box2D> const &domains = sparseBlock.getBulks();

  std::map<plint, Box2D>::const_iterator it = domains.begin();
  plint pos = 0;
  for (; it != domains.end(); ++it, ++pos) {
    plint id = it->first;
    if (threadAttribution.isLocal(id)) {
      AtomicContainerBlock2D &atomicDataContainer =
          bubbleContainer->getComponent(id);
      dynamic_cast<BubbleCounterData2D *>(atomicDataContainer.getData())
          ->reset();
    }
  }
}

void BubbleTracking2D::bubbleBucketFill(MultiScalarField2D<int> &flag) {
  setToConstant(*tagMatrix, tagMatrix->getBoundingBox(), (plint)-1);
  resetBubbleContainer();
  plint numIter = 2;
  while (numIter > 0) {
    CountBubbleIteration2D functional;
    std::vector<MultiBlock2D *> args;
    args.push_back(tagMatrix);
    args.push_back(&flag);
    args.push_back(bubbleContainer);
    applyProcessingFunctional(functional, bubbleContainer->getBoundingBox(),
                              args);
    plint numConflicts = functional.getNumConflicts();
    if (numConflicts > 0) {
      numIter = 2;
    } else {
      --numIter;
    }
  }
}

/* ************** class BubbleMPIdata ********************************** */

BubbleMPIdata2D::BubbleMPIdata2D(MultiBlock2D &block) {
  computeLocalIds(block);
}

std::vector<plint> const &BubbleMPIdata2D::getLocalIds() const {
  return localIds;
}

void BubbleMPIdata2D::computeLocalIds(MultiBlock2D &block) {
  MultiBlockManagement2D const &management = block.getMultiBlockManagement();
  ThreadAttribution const &threadAttribution =
      management.getThreadAttribution();
  SparseBlockStructure2D const &sparseBlock =
      management.getSparseBlockStructure();
  std::map<plint, Box2D> const &domains = sparseBlock.getBulks();

  std::map<plint, Box2D>::const_iterator it = domains.begin();
  plint pos = 0;
  for (; it != domains.end(); ++it, ++pos) {
    plint id = it->first;
    if (threadAttribution.isLocal(id)) {
      localIds.push_back(id);
    }
  }
}

/* ************** class BubbleCounterData2D **********************************
 */

BubbleCounterData2D::BubbleCounterData2D(plint maxNumBubbles_)
    : maxNumBubbles(maxNumBubbles_) {}

BubbleCounterData2D *BubbleCounterData2D::clone() const {
  return new BubbleCounterData2D(*this);
}

bool BubbleCounterData2D::convertCell(plint &tag0, plint tag1, plint tag2,
                                      plint tag3, plint tag4,

                                      plint tag1_, plint tag2_, plint tag3_,
                                      plint tag4_) {
  bool hasConflict =
      processNeighbor(tag0, tag1) || processNeighbor(tag0, tag2) ||
      processNeighbor(tag0, tag3) || processNeighbor(tag0, tag4) ||
      processNeighbor(tag0, tag1_) || processNeighbor(tag0, tag2_) ||
      processNeighbor(tag0, tag3_) || processNeighbor(tag0, tag4_);
  if (tag0 == -1) {
    tag0 = getNextTag();
  }
  return hasConflict;
}

bool BubbleCounterData2D::processNeighbor(plint &tag0, plint tag1) {
  plint myTag = convertTag(tag0);
  plint otherTag = convertTag(tag1);
  tag0 = myTag;  // re-assign tag in case it got re-assigned in "convertTag()".
  bool hasConflict = false;
  if (otherTag == -1) {
    if (myTag == -1) {
      tag0 = getNextTag();
    }
  } else {
    if (myTag < otherTag) {
      tag0 = otherTag;
      registerConflict(myTag, otherTag);
      hasConflict = true;
    }
  }
  return hasConflict;
}

plint BubbleCounterData2D::getNextTag() {
  plint nextTag = getUniqueID() * maxNumBubbles + nextCellId;
  ++nextCellId;
  return nextTag;
}

void BubbleCounterData2D::reset() {
  nextCellId = 0;
  retagging.clear();
  // uniqueTags.clear();
  // tagRemap.clear();
}

plint BubbleCounterData2D::convertTag(plint tag) const {
  if (tag == -1) return tag;
  std::map<plint, plint>::const_iterator it = retagging.find(tag);
  if (it == retagging.end()) {
    return tag;
  } else {
    return it->second;
  }
}

void BubbleCounterData2D::registerConflict(plint oldTag, plint newTag) {
  retagging.insert(std::pair<plint, plint>(oldTag, newTag));
  std::map<plint, plint>::iterator it = retagging.begin();
  // If some of the map's items are still pointing to the old
  // tag, redirect them to the new one.
  for (; it != retagging.end(); ++it) {
    if (it->second == oldTag) {
      it->second = newTag;
    }
  }
}

/* *************** Class BubbleRemapData2D ******************************** */

BubbleRemapData2D *BubbleRemapData2D::clone() const {
  return new BubbleRemapData2D(*this);
}

bool BubbleRemapData2D::isMyTag(plint tag) {
  return tag / maxNumBubbles == getUniqueID();
}

/* *************** Class CountBubbleIteration2D ********************************
 */

CountBubbleIteration2D::CountBubbleIteration2D()
    : numConflictsId(this->getStatistics().subscribeIntSum()) {}

CountBubbleIteration2D *CountBubbleIteration2D::clone() const {
  return new CountBubbleIteration2D(*this);
}

plint CountBubbleIteration2D::getNumConflicts() const {
  return this->getStatistics().getIntSum(numConflictsId);
}

void CountBubbleIteration2D::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 3);
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
  BubbleCounterData2D *pData =
      dynamic_cast<BubbleCounterData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleCounterData2D &data = *pData;

  Dot2D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      int currentFlag = flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y);
      if (currentFlag == freeSurfaceFlag2D::empty ||
          currentFlag == freeSurfaceFlag2D::interface) {
        bool isConflicting = data.convertCell(
            tagMatrix.get(iX, iY), tagMatrix.get(iX - 1, iY),
            tagMatrix.get(iX, iY - 1), tagMatrix.get(iX - 1, iY - 1),
            tagMatrix.get(iX - 1, iY + 1),

            tagMatrix.get(iX + 1, iY), tagMatrix.get(iX, iY + 1),
            tagMatrix.get(iX + 1, iY + 1), tagMatrix.get(iX + 1, iY - 1)

        );
        if (isConflicting) {
          statistics.gatherIntSum(numConflictsId, 1);
        }
      }
    }
  }
}

/* *************** Class CollectBubbleTags2D ******************************** */

CollectBubbleTags2D *CollectBubbleTags2D::clone() const {
  return new CollectBubbleTags2D(*this);
}

void CollectBubbleTags2D::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 2);
  ScalarField2D<plint> *pTagMatrix =
      dynamic_cast<ScalarField2D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField2D<plint> &tagMatrix = *pTagMatrix;

  AtomicContainerBlock2D *pDataBlock =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[1]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock2D &dataBlock = *pDataBlock;
  BubbleRemapData2D *pData =
      dynamic_cast<BubbleRemapData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleRemapData2D &data = *pData;

  std::set<plint> uniqueTags;
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint tag = tagMatrix.get(iX, iY);
      if (tag >= 0) {
        uniqueTags.insert(tag);
      }
    }
  }

  data.getUniqueTags().clear();
  std::set<plint>::const_iterator it = uniqueTags.begin();
  for (; it != uniqueTags.end(); ++it) {
    plint tag = *it;
    if (data.isMyTag(tag)) {
      data.getUniqueTags().push_back(tag);
    }
  }
}

/* *************** Class ApplyTagRemap2D ******************************** */

ApplyTagRemap2D *ApplyTagRemap2D::clone() const {
  return new ApplyTagRemap2D(*this);
}

void ApplyTagRemap2D::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 2);
  ScalarField2D<plint> *pTagMatrix =
      dynamic_cast<ScalarField2D<plint> *>(atomicBlocks[0]);
  PLB_ASSERT(pTagMatrix);
  ScalarField2D<plint> &tagMatrix = *pTagMatrix;

  AtomicContainerBlock2D *pDataBlock =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[1]);
  PLB_ASSERT(pDataBlock);
  AtomicContainerBlock2D &dataBlock = *pDataBlock;
  BubbleRemapData2D *pData =
      dynamic_cast<BubbleRemapData2D *>(dataBlock.getData());
  PLB_ASSERT(pData);
  BubbleRemapData2D &data = *pData;

  std::map<plint, plint> const &tagRemap = data.getTagRemap();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint tag = tagMatrix.get(iX, iY);
      if (tag >= 0) {
        std::map<plint, plint>::const_iterator it = tagRemap.find(tag);
        PLB_ASSERT(it != tagRemap.end());
        tagMatrix.get(iX, iY) = it->second;
      }
    }
  }
}

}  // namespace lbfoam
}  // namespace plb
