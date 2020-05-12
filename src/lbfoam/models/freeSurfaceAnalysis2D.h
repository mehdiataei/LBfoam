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

#ifndef FREE_SURFACE_ANALYSIS_2D_H
#define FREE_SURFACE_ANALYSIS_2D_H

#include <memory>

#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock2D *> freeSurfaceArgs,
                         Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageMassFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_AverageMassFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_AverageMassFunctional2D<T, Descriptor> *clone() const;
  T getAverageMass() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint averageMassId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock2D *> freeSurfaceArgs,
                       Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_TotalMassFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_TotalMassFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_TotalMassFunctional2D<T, Descriptor> *clone() const;
  T getTotalMass() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint totalMassId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock2D *> freeSurfaceArgs,
                            Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageDensityFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_AverageDensityFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_AverageDensityFunctional2D<T, Descriptor> *clone() const;
  T getAverageDensity() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint averageDensityId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock2D *> freeSurfaceArgs,
                                   Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageVolumeFractionFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_AverageVolumeFractionFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_AverageVolumeFractionFunctional2D<T, Descriptor> *clone() const;
  T getAverageVolumeFraction() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint averageVfId;
};

template <typename T, template <typename U> class Descriptor>
plint countFreeSurfaceElements(std::vector<MultiBlock2D *> freeSurfaceArgs,
                               plint flagToLookFor, Box2D domain);

template <typename T, template <typename U> class Descriptor>
class CountFreeSurfaceElementsFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  CountFreeSurfaceElementsFunctional2D(plint flagToLookFor_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual CountFreeSurfaceElementsFunctional2D<T, Descriptor> *clone() const;
  plint getNumElements() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint numCellsId;
  plint flagToLookFor;
};

template <typename T, template <typename U> class Descriptor>
Array<T, 2> freeSurfaceAverageMomentum(
    std::vector<MultiBlock2D *> freeSurfaceArgs, Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageMomentumFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_AverageMomentumFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_AverageMomentumFunctional2D<T, Descriptor> *clone() const;
  Array<T, 2> getAverageMomentum() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  Array<plint, 2> averageMomentumId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock2D *> freeSurfaceArgs,
                           Box2D domain);

template <typename T, template <typename U> class Descriptor>
class FS_AverageHeightFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FS_AverageHeightFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FS_AverageHeightFunctional2D<T, Descriptor> *clone() const;
  T getAverageHeight() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint averageHeightId;
};

template <typename T, template <typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock2D *> freeSurfaceArgs, plint N,
                       Box2D domain);

template <typename T, template <typename U> class Descriptor>
class GetWaterLevelAtxyFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  GetWaterLevelAtxyFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual GetWaterLevelAtxyFunctional2D<T, Descriptor> *clone() const;
  plint getNumFluidCellsAtXY() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    for (pluint i = 0; i < modified.size(); ++i) {
      modified[i] = modif::nothing;
    }
  }

 private:
  plint numFluidOccupiedCellId;
};

// Here we count both the no-slip and the free-slip wall cells.
class FreeSurfaceCountWallCells2D
    : public ReductiveBoxProcessingFunctional2D_S<int> {
 public:
  FreeSurfaceCountWallCells2D();
  virtual void process(Box2D domain, ScalarField2D<int> &flag);
  virtual FreeSurfaceCountWallCells2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // Flags.
  }
  plint getNumWallCells() const;

 private:
  plint sumScalarId;
};

plint freeSurfaceCountWallCells(MultiScalarField2D<int> &flag, Box2D domain);

template <typename T>
T freeSurfaceComputePorosity(MultiScalarField2D<int> &flag, Box2D domain);

template <typename T>
class FreeSurfaceComputeFluidVolume2D
    : public ReductiveBoxProcessingFunctional2D_SS<T, int> {
 public:
  FreeSurfaceComputeFluidVolume2D();
  virtual void process(Box2D domain, ScalarField2D<T> &volumeFraction,
                       ScalarField2D<int> &flag);
  virtual FreeSurfaceComputeFluidVolume2D<T> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // Volume fraction.
    modified[1] = modif::nothing;  // Flags.
  }
  T getFluidVolume() const;

 private:
  plint sumScalarId;
};

template <typename T>
class MaskedFreeSurfaceComputeFluidVolume2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  MaskedFreeSurfaceComputeFluidVolume2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual MaskedFreeSurfaceComputeFluidVolume2D<T> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // Volume fraction.
    modified[1] = modif::nothing;  // Flags.
    modified[2] = modif::nothing;  // Mask.
  }
  T getFluidVolume() const;

 private:
  plint sumScalarId;
};

template <typename T>
T freeSurfaceComputeFluidVolume(MultiScalarField2D<T> &volumeFraction,
                                MultiScalarField2D<int> &flag, Box2D domain);

template <typename T>
T freeSurfaceComputeFluidVolume(MultiScalarField2D<T> &volumeFraction,
                                MultiScalarField2D<int> &flag,
                                MultiScalarField2D<int> &mask);

// The porosity must have been computed inside the exact same domain as the
// argument of this function. Usually the porosity is constant in time, this is
// why we pass it as an argument, and do not recompute it inside this function.
template <typename T>
T freeSurfaceComputeSaturation(T porosity,
                               MultiScalarField2D<T> &volumeFraction,
                               MultiScalarField2D<int> &flag, Box2D domain);

// Here we find the smallest Box2D that includes all fluid and interface cells
class FreeSurfaceComputeFluidBoundingBox2D
    : public ReductiveBoxProcessingFunctional2D_S<int> {
 public:
  FreeSurfaceComputeFluidBoundingBox2D();
  virtual void process(Box2D domain, ScalarField2D<int> &flag);
  virtual FreeSurfaceComputeFluidBoundingBox2D *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;  // Flags.
  }
  Box2D getFluidBoundingBox() const;

 private:
  plint minIdX, minIdY;
  plint maxIdX, maxIdY;
};

Box2D freeSurfaceComputeFluidBoundingBox(MultiScalarField2D<int> &flag,
                                         Box2D domain);
Box2D freeSurfaceComputeFluidBoundingBox(MultiScalarField2D<int> &flag);

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceBoxSumForcedEnergyFunctional2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  FreeSurfaceBoxSumForcedEnergyFunctional2D();
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor> *clone()
      const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;
  T getAverageEnergy() const;

 private:
  plint sumEnergyId;
  plint sumCellsId;
};

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain);

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, MultiScalarField2D<T> &velocityNorm,
    Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, MultiScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain, plint iComponent);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, plint iComponent);

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag,
    MultiTensorField2D<T, Descriptor<T>::d> &velocity, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T, Descriptor<T>::d>>
freeSurfaceComputeForcedVelocity(MultiBlockLattice2D<T, Descriptor> &lattice,
                                 MultiTensorField2D<T, Descriptor<T>::d> &force,
                                 MultiScalarField2D<int> &flag, Box2D domain);

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T, Descriptor<T>::d>>
freeSurfaceComputeForcedVelocity(MultiBlockLattice2D<T, Descriptor> &lattice,
                                 MultiTensorField2D<T, Descriptor<T>::d> &force,
                                 MultiScalarField2D<int> &flag);

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_2D_H
