/*
"#############################################################################"
"                                                                             "
"  LBfoam: An open-source software package for the simulation of foaming      "
"  using the Lattice Boltzmann Method               	                        "
"  Copyright (C) 2020 Mohammadmehdi Ataei                                     "
"  m.ataei@mail.utoronto.ca                                                   "
"                                                                             "
"  This file is part of LBfoam.                                               "
"                                                                             "
"  LBfoam is free software: you can redistribute it and/or modify it under    "
"  the terms of the GNU Affero General Public License as published by the     "
"  Free Software Foundation version 3.                                        "
"                                                                             "
"  LBfoam is distributed in the hope that it will be useful, but WITHOUT ANY  "
"  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  "
"  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for   "
"  more details.                                                              "
"                                                                             "
"  You should have received a copy of the GNU Affero General Public License   "
"  along with this Program. If not, see <http://www.gnu.org/licenses/>.       "
"                                                                             "
"#############################################################################"
"                                                                             "
"  Author:         Mohammadmehdi Ataei, 2020                                  "
"                                                                             "
"#############################################################################"
"                                                                             "
"  Parts of the LBfoam code that originate from Palabos are distributed       "
"  under the terms of the AGPL 3.0 license with the following copyright       "
"  notice:                                                                    "
"                                                                             "
"  This file is part of the Palabos library.                                  "
"                                                                             "
"  Copyright (C) 2011-2017 FlowKit Sarl                                       "
"  Route d'Oron 2                                                             "
"  1010 Lausanne, Switzerland                                                 "
"  E-mail contact: contact@flowkit.com                                        "
"                                                                             "
"  The most recent release of Palabos can be downloaded at                    "
"  <http://www.palabos.org/>                                                  "
"                                                                             "
"  The library Palabos is free software: you can redistribute it and/or       "
"  modify it under the terms of the GNU Affero General Public License as      "
"  published by the Free Software Foundation, either version 3 of the         "
"  License, or (at your option) any later version.                            "
"                                                                             "
"  The library is distributed in the hope that it will be useful,             "
"  but WITHOUT ANY WARRANTY; without even the implied warranty of             "
"  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               "
"  GNU Affero General Public License for more details.                        "
"                                                                             "
"  You should have received a copy of the GNU Affero General Public License   "
"  along with this program.  If not, see <http://www.gnu.org/licenses/>.      "
"                                                                             "
"                                                                             "
"                                                                             "
"#############################################################################"
*/

#ifndef BUBBLE_TRACKING_2D_H
#define BUBBLE_TRACKING_2D_H

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"
#include "lbfoam/models/freeSurfaceModel2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"

namespace plb

{

namespace lbfoam {

class BubbleMPIdata2D {
 public:
  BubbleMPIdata2D(MultiBlock2D &block);
  std::vector<plint> const &getLocalIds() const;

 private:
  void computeLocalIds(MultiBlock2D &block);

 private:
  std::vector<plint> localIds;
};

class BubbleTracking2D {
 public:
  BubbleTracking2D(MultiBlock2D &templ);
  ~BubbleTracking2D();
  template <typename T, template <typename U1> class FluidDescriptor,
            template <typename U2> class ScalarDescriptor>
  void execute(MultiScalarField2D<T> &volumeFraction,
               MultiScalarField2D<int> &flag, MultiTensorField2D<T, 2> &normal,
               MultiScalarField2D<T> &rhoBar, MultiScalarField2D<T> &mass,
               MultiTensorField2D<T, 2> &j,
               MultiBlockLattice2D<T, ScalarDescriptor> &adLattice,
               MultiScalarField2D<T> &oldVolumeFraction,
               bool incompressibleModel, T pi);

  MultiScalarField2D<plint> *getTagMatrix() { return tagMatrix; }

  MultiScalarField2D<double> *getDisjoiningPressureField() {
    return disjoiningPressureField;
  }

  void setTagMatrix(MultiScalarField2D<plint> *newTagMatrix) {
    tagMatrix = newTagMatrix;
  }
  void setDisjoiningPressureField(MultiScalarField2D<double> *newDisjoining) {
    disjoiningPressureField = newDisjoining;
  }
  std::vector<double> const &getBubbleVolume() { return bubbleVolume; }
  std::vector<double> const &getBubbleDensity() { return bubbleDensity; }
  std::vector<double> const &getBubbleDisjoiningPressure() {
    return bubbleDisjoiningPressure;
  }
  std::vector<Array<double, 2>> const &getBubbleCenter() {
    return bubbleCenter;
  }
  pluint numBubbles() const { return bubbleVolume.size(); }

 private:
  // Re-assign a continuously numbered ID to the detected bubbles.
  pluint countAndTagBubbles();
  // Computes the volumes and centers of all new bubbles.
  template <typename T>
  void bubbleAnalysis(MultiScalarField2D<int> &flag,
                      MultiScalarField2D<T> &volumeFraction, pluint numBubbles);
  template <typename T, template <typename U1> class AD_Descriptor>
  void calculateBubbleGrowth(MultiScalarField2D<int> &flag,
                             MultiScalarField2D<T> &volumeFraction,
                             MultiScalarField2D<T> &oldVolumeFraction,
                             MultiBlockLattice2D<T, AD_Descriptor> &adLattice,
                             pluint numBubbles);
  // calculates DisjoiningPressure
  template <typename T, template <typename U> class Descriptor>
  void calculateDisjoiningPressure(T pi, MultiScalarField2D<int> &flag,
                                   MultiScalarField2D<T> &volumeFraction,
                                   MultiTensorField2D<T, 2> &normal,
                                   MultiScalarField2D<T> &rhoBar,
                                   MultiScalarField2D<T> &mass,
                                   MultiTensorField2D<T, 2> &j,
                                   bool incompressibleModel, pluint numBubbles);

  // template<typename T, template<typename U> class Descriptor>
  // void stabilizeLattice(FreeSurfaceFields2D<T,Descriptor>& fields);
  // Implements all required MPI operations needed to compute the bubble volume
  // and centers,
  // after calling AnalyzeBubbles2D.
  void computeBubbleData(pluint numBubbles);
  // Implements all required MPI operations needed to compute the global IDs of
  // the current
  // bubbbles, after calling CollectBubbleTags2D.
  plint globalBubbleIds();
  // Prepare the bubble container for the next time iteration.
  void resetBubbleContainer();
  // Parallel bucket-fill algorithm to assign a unique ID to every contiguous
  // region.
  void bubbleBucketFill(MultiScalarField2D<int> &flag);

 private:
  BubbleTracking2D(BubbleTracking2D const &rhs) : mpiData(rhs.mpiData) {
    PLB_ASSERT(false);
  }
  BubbleTracking2D &operator=(BubbleTracking2D const &rhs) {
    PLB_ASSERT(false);
    return *this;
  }

 private:
  MultiContainerBlock2D *bubbleContainer, *bubbleAnalysisContainer,
      *bubbleRemapContainer;
  BubbleMPIdata2D mpiData;
  MultiScalarField2D<plint> *tagMatrix;
  MultiScalarField2D<double> *disjoiningPressureField;
  std::vector<double> bubbleVolume;
  std::vector<double> bubbleDensity;
  std::vector<double> bubbleDisjoiningPressure;
  std::vector<Array<double, 2>> bubbleCenter;
  static const plint maxNumBubbles = 100000;
};

/**
 * Data for the bubble counter, associated to one block.
 * It holds information that changes during time:
 *  nextCellId: next available ID, if a new cell-type is found.
 *  retagging: a map that relates equivalent IDs, when two domains
 *             are found to be contiguous.
 *  maxNumBubbles: an upper bound for the number of bubbles per block,
 *                 so every block can create a globally unique bubble ID.
 **/
class BubbleCounterData2D : public ContainerBlockData {
 public:
  BubbleCounterData2D(plint maxNumBubbles_);
  virtual BubbleCounterData2D *clone() const;
  // The assumption here is that cell 0 is a cell that needs to
  // be tagged ("a bubble cell"). Depending on the tag of neighboring
  // cells, either one of the neighbor's tag or a new tag is assigned.
  // There is a conflict if cell 0 has been previously converted in
  // a way that is incompatible with the neighbor's tags.
  bool convertCell(plint &tag0, plint tag1, plint tag2, plint tag3, plint tag4,

                   plint tag1_, plint tag2_, plint tag3_, plint tag4_);
  // The assumption here is that cell 0 is a cell that needs to
  // be tagged ("a bubble cell"). Depending on the tag of neighboring
  // cell 1, either cell 1's tag or a new tag is assigned.
  // There is a conflict if cell 0 has been previously converted in
  // a way that is incompatible with cell 1's tag.
  bool processNeighbor(plint &tag0, plint tag1);
  plint getNextTag();
  void reset();
  // Get the "real" value (after remapping) of a given tag.
  plint convertTag(plint tag) const;
  // Important: it is assumed that oldTag and newTag are not
  // themselves mapped. This means that for arbitrary tags tag1 and
  // tag2 you should not call registerConflict(tag1,tag2), but
  // registerConflict(convertTag(tag1), convertTag(tag2).

  void registerConflict(plint oldTag, plint newTag);

 private:
  plint nextCellId;
  std::map<plint, plint> retagging;

 private:
  plint maxNumBubbles;
};

class BubbleRemapData2D : public ContainerBlockData {
 public:
  BubbleRemapData2D(plint maxNumBubbles_ = 0) : maxNumBubbles(maxNumBubbles_) {}
  virtual BubbleRemapData2D *clone() const;
  std::vector<plint> &getUniqueTags() { return uniqueTags; }
  std::vector<plint> const &getUniqueTags() const { return uniqueTags; }
  std::map<plint, plint> &getTagRemap() { return tagRemap; }
  bool isMyTag(plint tag);

 private:
  plint maxNumBubbles;
  std::vector<plint> uniqueTags;
  std::map<plint, plint> tagRemap;
};

struct BubbleAnalysisData2D : public ContainerBlockData {
  virtual BubbleAnalysisData2D *clone() const {
    return new BubbleAnalysisData2D(*this);
  }
  std::vector<double> bubbleVolume;
  std::vector<double> bubbleDensity;
  std::vector<double> bubbleDisjoiningPressure;
  std::vector<Array<double, 2>> bubbleCenter;
};

class CountBubbleIteration2D : public PlainReductiveBoxProcessingFunctional2D {
 public:
  CountBubbleIteration2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual CountBubbleIteration2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // tags.
    modified[1] = modif::nothing;          // flags.
    modified[2] = modif::nothing;          // data.
  }
  plint getNumConflicts() const;

 private:
  plint numConflictsId;
};

template <typename T>
class AnalyzeBubbles2D : public BoxProcessingFunctional2D {
 public:
  AnalyzeBubbles2D(pluint numBubbles_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual AnalyzeBubbles2D<T> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // tags.
    modified[1] = modif::nothing;          // flags.
    modified[2] = modif::nothing;          // data.
    modified[3] = modif::nothing;          // volume fraction.
  }

 private:
  pluint numBubbles;
};

// Converts the information about the overall available bubble tags available
// to all processors.
class CollectBubbleTags2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual CollectBubbleTags2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // tags.
    modified[1] = modif::nothing;  // data.
  }
};

// Assign a new tag to all bubble cells (they must have been uniquely tagged
// previously).
// The only field in the BubbleCounterData2D which is used here is tagRemap.
class ApplyTagRemap2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual ApplyTagRemap2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // tags.
    modified[1] = modif::nothing;          // data.
  }
};

template <typename T, template <typename U1> class AD_Descriptor>
class CalculateBubbleGrowth2D : public BoxProcessingFunctional2D {
 public:
  CalculateBubbleGrowth2D(pluint numBubbles_);
  //~CalculateBubbleGrowth2D();
  CalculateBubbleGrowth2D(CalculateBubbleGrowth2D<T, AD_Descriptor> const &rhs);
  CalculateBubbleGrowth2D<T, AD_Descriptor> &operator=(
      CalculateBubbleGrowth2D<T, AD_Descriptor> const &rhs);
  void swap(CalculateBubbleGrowth2D<T, AD_Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual CalculateBubbleGrowth2D<T, AD_Descriptor> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  pluint numBubbles;
};

template <typename T, template <typename U> class Descriptor>
class CalculateDisjoiningPressure2D : public BoxProcessingFunctional2D {
 public:
  CalculateDisjoiningPressure2D(T pi_, pluint numBubbles_,
                                bool incompressibleModel);
  CalculateDisjoiningPressure2D(
      CalculateDisjoiningPressure2D<T, Descriptor> const &rhs);
  CalculateDisjoiningPressure2D<T, Descriptor> &operator=(
      CalculateDisjoiningPressure2D<T, Descriptor> const &rhs);
  void swap(CalculateDisjoiningPressure2D<T, Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual CalculateDisjoiningPressure2D<T, Descriptor> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  T pi;
  pluint numBubbles;
  bool incompressibleModel;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // BUBBLE_TRACKING_2D_H
