/*
"#############################################################################"
"                                                                             "
"  LBfoam: An open-source software package for the simulation of foaming      "
"  using the Lattice Boltzmann Method               	                      "
"  Copyright (C) 2020 Mohammadmehdi Ataei                                     "
"  m.ataei@mail.utoronto.ca                                                  "
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

#ifndef BUBBLE_GROWTH_2D_H
#define BUBBLE_GROWTH_2D_H

#include <sstream>
#include <string>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"
#include "lbfoam/bubble/bubbleTracking2D.h"
#include "lbfoam/models/freeSurfaceModel2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"

namespace plb {

namespace lbfoam {

struct BubbleTransition2D;
class BubbleInfo2D;
struct FullBubbleRecord2D;

template <typename T>
class BubbleGrowth2D {
 public:
  BubbleGrowth2D(MultiBlock2D &templ);
  ~BubbleGrowth2D();
  void transition(BubbleTracking2D &bubbleTracking, plint iterationStep,
                  T temperature, T R_s, T p_sat, T dx, T fluidRho,
                  T newBubbleVolumeCorrection = (T)1.,
                  bool entrapBubbles = false, plint numRememberedVolumes = 1);
  // Based on the bubble information (original volume, current volume), update
  // the pressure field for all bubbles.
  void updateBubblePressure(MultiScalarField2D<T> &outsideDensity, T rhoEmpty,
                            T alpha = (T)-1, T beta = (T)-1, T gamma = (T)1);
  void updateBubbleGrowth(MultiScalarField2D<T> &outsideDensity, T rhoEmpty,
                          T alpha = (T)-1, T beta = (T)-1, T gamma = (T)1);
  void scaleAndIncreaseDensity(T temperature, T R_s, T p_sat, T dx, T fluidRho);
  MultiScalarField2D<plint> *getOldTagMatrix() { return oldTagMatrix; }
  void freeze();
  void freezeLargestBubble();
  void timeHistoryLog(std::string fName);
  void fullBubbleLog(std::string fName);
  std::map<plint, BubbleInfo2D> const &getBubbles() const { return bubbles; }

 private:
  // Implement the time evolution of the bubbles:
  // - Reassign an ID to the new bubbles, compatible with the old ID.
  // - Create a data structure for the new bubbles.
  // - Print a message if bubbles are created/deleted.
  void matchAndRemapBubbles(BubbleTracking2D &bubbleTracking,
                            MultiScalarField2D<plint> &tagMatrix1,
                            MultiScalarField2D<plint> &tagMatrix2,
                            T newBubbleVolumeCorrection, plint iterationStep,
                            bool entrapBubbles, plint numRememberedVolumes);
  // Create the map "newToAllOldMap" which maps the temporary new IDs
  // (continuously numbered) to the "old" IDs (either existing ones or newly
  // attributed).
  void correlateBubbleIds(MultiScalarField2D<plint> &tagMatrix1,
                          MultiScalarField2D<plint> &tagMatrix2,
                          std::vector<std::vector<plint>> &newToAllOldMap,
                          pluint numBubbles);
  // Updates the (non-parallel) data structure which holds the overview of the
  // currently available bubbles.
  void updateBubbleInformation(
      std::vector<BubbleTransition2D> &bubbleTransitions,
      std::vector<double> const &bubbleVolume,
      std::vector<double> const &bubbleDensity,
      std::vector<double> const &bubbleDisjoiningPressure,
      T newBubbleVolumeCorrection, plint iterationStep, bool entrapBubbles,
      plint numRememberedVolumes);
  void computeNewBubbles(std::set<plint> &oldIDs, std::set<plint> &newIDs,
                         std::vector<double> const &bubbleVolume,
                         std::vector<double> const &bubbleDensity,
                         std::vector<double> const &bubbleDisjoiningPressure,
                         T newBubbleVolumeCorrection,
                         std::map<plint, BubbleInfo2D> &newBubbles,
                         std::map<plint, plint> &newToFinal, bool entrapBubbles,
                         plint numRememberedVolumes);
  void updateBubbleLog(BubbleTransition2D &bubbleTransition,
                       std::vector<double> const &bubbleVolume,
                       plint iterationStep,
                       std::map<plint, BubbleInfo2D> &newBubbles,
                       std::map<plint, plint> &newToFinal);
  static void computeBubbleTransitions(
      std::vector<std::vector<plint>> const &newToAllOldMap,
      std::map<plint, std::vector<plint>> const &oldToAllNewMap,
      std::vector<BubbleTransition2D> &bubbleTransitions);
  static T getCutoffVolume() { return 2.0; }

 private:
  BubbleGrowth2D(BubbleGrowth2D<T> const &rhs) { PLB_ASSERT(false); }
  BubbleGrowth2D<T> &operator=(BubbleGrowth2D<T> const &rhs) {
    PLB_ASSERT(false);
    return *this;
  }

 private:
  MultiContainerBlock2D *bubbleAnalysisContainer, *bubbleCorrelationContainer,
      *bubbleRemapContainer;
  BubbleMPIdata2D mpiData;
  MultiScalarField2D<plint> *oldTagMatrix;
  std::map<plint, BubbleInfo2D> bubbles;
  plint nextBubbleID;
  //    <Iteration>   <IDs of created bubbles>,<IDs of vanished bubbles>
  std::map<plint, std::pair<std::vector<plint>, std::vector<plint>>>
      timeHistory;
  std::vector<FullBubbleRecord2D> FullBubbleRecord;
  /// IMPORTANT ///
  // Checkpointing needs to be implemented for this code. To checkpoint, you
  // need to save
  // - bubbles
  // - FullBubbleRecord2D
  // - timeHistory (optional)
  // - oldTagMatrix
  // - nextBubbleID
};

struct BubbleCorrelationData2D : public ContainerBlockData {
  virtual BubbleCorrelationData2D *clone() const {
    return new BubbleCorrelationData2D(*this);
  }
  // Maps new bubble tags to old bubble tags at successive iterations.
  // This is needed because a new tag may map to several old tags, due
  // to a bubble merging process.
  std::vector<plint> newToOldMap0, newToOldMap1;
};

class BubbleInfo2D {
 public:
  BubbleInfo2D()
      : referenceVolume(0.),
        currentVolume(0.),
        frozen(false),
        currentDensity((double)1.),
        incrim(0.),
        disjoiningPressure(0.)

  {
    numRememberedVolumes = 0;
    rememberedVolumes.clear();
    invNumRememberedVolumes = 0.0;
    numRegisteredVolumes = 0;
    meanVolume = 0.0;
  }
  BubbleInfo2D(double volume, plint numRememberedVolumes_, double density)
      : referenceVolume(volume),
        currentVolume(volume),
        frozen(false),
        currentDensity(density),
        incrim(0.),
        disjoiningPressure(0.),
        numRememberedVolumes(numRememberedVolumes_)

  {
    PLB_ASSERT(numRememberedVolumes > 0);

    if (numRememberedVolumes != 1) {
      rememberedVolumes.resize(numRememberedVolumes, referenceVolume);
      invNumRememberedVolumes = 1.0 / numRememberedVolumes;
      numRegisteredVolumes = 0;
      meanVolume = referenceVolume;
    } else {
      rememberedVolumes.clear();
      invNumRememberedVolumes = 0.0;
      numRegisteredVolumes = 0;
      meanVolume = referenceVolume;
    }
  }
  void freeze() { frozen = true; }

  void unfreeze() { frozen = false; }

  void setVolume(double newVolume) {
    if (!frozen) {
      currentVolume = newVolume;

      if (numRememberedVolumes != 1) {
        plint i = numRegisteredVolumes % numRememberedVolumes;
        numRegisteredVolumes++;
        meanVolume +=
            invNumRememberedVolumes * (currentVolume - rememberedVolumes[i]);
        rememberedVolumes[i] = currentVolume;
      } else {
        meanVolume = currentVolume;
      }
    }
  }
  void increaseDensity(double incrim) {
    if (!frozen) {
      currentDensity += incrim;
    }
  }

  void setDensity(double density) { currentDensity = density; }

  void setIncrim(double inc)

  {
    incrim = inc;
  }

  void setDisjoiningPressure(double dp) { disjoiningPressure = dp; }

  double getVolumeRatio() const {
    static const double epsilon = std::numeric_limits<double>::epsilon() * 1.e4;
    if (std::fabs(meanVolume) > epsilon) {
      return referenceVolume / meanVolume;
    } else {
      return 1.0;
    }
  }
  double getReferenceVolume() const { return referenceVolume; }
  void setReferenceVolume(double vol) { referenceVolume = vol; }
  double getCurrentDensity() const { return currentDensity; }
  double getVolume() const { return currentVolume; }
  double getIncrim() const { return incrim; }
  double getDisjoiningPressure() const { return disjoiningPressure; }
  bool isFrozen() const { return frozen; }

 private:
  double referenceVolume;
  double currentVolume;
  bool frozen;
  double currentDensity;
  double incrim;
  double disjoiningPressure;

  // The following variables relate to the hysteresis mechanism.
  plint numRememberedVolumes;
  std::vector<double> rememberedVolumes;
  double invNumRememberedVolumes;
  plint numRegisteredVolumes;
  double meanVolume;
};

struct BubbleTransition2D {
  bool empty() const { return oldIDs.empty() && newIDs.empty(); }
  std::string description() const {
    std::stringstream sstream;
    if (oldIDs.empty() && newIDs.empty()) {
    } else if (oldIDs.empty()) {
      PLB_ASSERT(newIDs.size() == 1);
      sstream << "Bubble " << *newIDs.begin() << " created";
    } else if (newIDs.empty()) {
      PLB_ASSERT(oldIDs.size() == 1);
      sstream << "Bubble " << *oldIDs.begin() << " vanished";
    } else if (oldIDs.size() == 1) {
      if (newIDs.size() == 1) {
        sstream << "Straight transition from ID " << *oldIDs.begin()
                << " to ID " << *newIDs.begin();
      } else {
        sstream << "Splitting from ID " << *oldIDs.begin() << " to IDs";
        std::set<plint>::const_iterator it = newIDs.begin();
        for (; it != newIDs.end(); ++it) {
          sstream << " " << *it;
        }
      }
    } else if (newIDs.size() == 1) {
      sstream << "Merging IDs";
      std::set<plint>::const_iterator it = oldIDs.begin();
      for (; it != oldIDs.end(); ++it) {
        sstream << " " << *it;
      }
      sstream << " into ID " << *newIDs.begin();
    } else {
      sstream << "Transition from IDs";
      std::set<plint>::const_iterator it1 = oldIDs.begin();
      for (; it1 != oldIDs.end(); ++it1) {
        sstream << " " << *it1;
      }
      sstream << " into IDs ";
      std::set<plint>::const_iterator it2 = newIDs.begin();
      for (; it2 != newIDs.end(); ++it2) {
        sstream << " " << *it2;
      }
    }
    return sstream.str();
  }
  std::set<plint> oldIDs, newIDs;
};

struct FullBubbleRecord2D {
  FullBubbleRecord2D(double initialVolume_, plint beginIteration_)
      : initialVolume(initialVolume_),
        finalVolume(0.),
        beginIteration(beginIteration_),
        endIteration(beginIteration_),
        frozen(false) {}

  std::string description(plint ID) const {
    std::stringstream sstream;
    sstream << "Bubble " << ID << ". At it. " << beginIteration << ", Vol. "
            << initialVolume << ": " << beginTransition.description() << ".";
    if (endIteration > beginIteration) {
      sstream << " Removed at it. " << endIteration << ", Vol. " << finalVolume
              << ": " << endTransition.description() << ".";
    }
    sstream << " Frozen bubble: " << (frozen ? "Yes." : "No.");
    return sstream.str();
  }

  double initialVolume, finalVolume;
  plint beginIteration, endIteration;
  BubbleTransition2D beginTransition, endTransition;
  bool frozen;
};

template <typename T>
class CorrelateBubbleIds2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual CorrelateBubbleIds2D<T> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // tags before.
    modified[1] = modif::nothing;  // tags after.
    modified[2] = modif::nothing;  // data.
  }
};

template <typename T>
class UpdateBubblePressure2D : public BoxProcessingFunctional2D_SS<plint, T> {
 public:
  UpdateBubblePressure2D(std::map<plint, BubbleInfo2D> const &bubbles_, T rho0_,
                         T alpha_, T beta_, T gamma_);
  virtual void process(Box2D domain, ScalarField2D<plint> &tags,
                       ScalarField2D<T> &density);
  virtual UpdateBubblePressure2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;          // tags.
    modified[1] = modif::staticVariables;  // density.
  }

 private:
  static T getCutoffDensityRatio() { return 2.0; }
  // static T getCutoffDensityRatio() { return 1.01; }
 private:
  std::map<plint, BubbleInfo2D> bubbles;
  T rho0;
  T alpha;
  T beta;
  T gamma;
};

/* Mehdi*/
template <typename T>
class UpdateBubbleGrowth2D : public BoxProcessingFunctional2D_SS<plint, T> {
 public:
  UpdateBubbleGrowth2D(std::map<plint, BubbleInfo2D> const &bubbles_, T rho0_,
                       T alpha_, T beta_, T gamma_);
  virtual void process(Box2D domain, ScalarField2D<plint> &tags,
                       ScalarField2D<T> &density);
  virtual UpdateBubbleGrowth2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;          // tags.
    modified[1] = modif::staticVariables;  // density.
  }

 private:
  static T getCutoffDensityRatio() { return 2.0; }

 private:
  std::map<plint, BubbleInfo2D> bubbles;
  T rho0;
  T alpha;
  T beta;
  T gamma;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // BUBBLE_GROWTH_2D_H
