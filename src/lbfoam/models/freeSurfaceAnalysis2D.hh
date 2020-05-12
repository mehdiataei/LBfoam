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

#ifndef FREE_SURFACE_ANALYSIS_2D_HH
#define FREE_SURFACE_ANALYSIS_2D_HH

#include "lbfoam/models/freeSurfaceAnalysis2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
FS_AverageMassFunctional2D<T, Descriptor>::FS_AverageMassFunctional2D()
    : averageMassId(this->getStatistics().subscribeAverage()) {}

template <typename T, template <typename U> class Descriptor>
void FS_AverageMassFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      statistics.gatherAverage(averageMassId, param.mass(iX, iY));
      statistics.incrementStats();
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMassFunctional2D<T, Descriptor>
    *FS_AverageMassFunctional2D<T, Descriptor>::clone() const {
  return new FS_AverageMassFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageMassFunctional2D<T, Descriptor>::getAverageMass() const {
  return this->getStatistics().getAverage(averageMassId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock2D *> freeSurfaceArgs,
                         Box2D domain) {
  FS_AverageMassFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getAverageMass();
}

template <typename T, template <typename U> class Descriptor>
FS_TotalMassFunctional2D<T, Descriptor>::FS_TotalMassFunctional2D()
    : totalMassId(this->getStatistics().subscribeSum()) {}

template <typename T, template <typename U> class Descriptor>
void FS_TotalMassFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      statistics.gatherSum(totalMassId, param.mass(iX, iY));
      statistics.incrementStats();
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_TotalMassFunctional2D<T, Descriptor>
    *FS_TotalMassFunctional2D<T, Descriptor>::clone() const {
  return new FS_TotalMassFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_TotalMassFunctional2D<T, Descriptor>::getTotalMass() const {
  return this->getStatistics().getSum(totalMassId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock2D *> freeSurfaceArgs,
                       Box2D domain) {
  FS_TotalMassFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getTotalMass();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageDensityFunctional2D<T, Descriptor>::FS_AverageDensityFunctional2D()
    : averageDensityId(this->getStatistics().subscribeAverage()) {}

template <typename T, template <typename U> class Descriptor>
void FS_AverageDensityFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.cell(iX, iY).getDynamics().hasMoments()) {
        statistics.gatherAverage(averageDensityId, param.getDensity(iX, iY));
        statistics.incrementStats();
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageDensityFunctional2D<T, Descriptor>
    *FS_AverageDensityFunctional2D<T, Descriptor>::clone() const {
  return new FS_AverageDensityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageDensityFunctional2D<T, Descriptor>::getAverageDensity() const {
  return this->getStatistics().getAverage(averageDensityId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock2D *> freeSurfaceArgs,
                            Box2D domain) {
  FS_AverageDensityFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getAverageDensity();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageVolumeFractionFunctional2D<
    T, Descriptor>::FS_AverageVolumeFractionFunctional2D()
    : averageVfId(this->getStatistics().subscribeAverage()) {}

template <typename T, template <typename U> class Descriptor>
void FS_AverageVolumeFractionFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.cell(iX, iY).getDynamics().hasMoments()) {
        statistics.gatherAverage(averageVfId, param.volumeFraction(iX, iY));
        statistics.incrementStats();
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageVolumeFractionFunctional2D<T, Descriptor>
    *FS_AverageVolumeFractionFunctional2D<T, Descriptor>::clone() const {
  return new FS_AverageVolumeFractionFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageVolumeFractionFunctional2D<
    T, Descriptor>::getAverageVolumeFraction() const {
  return this->getStatistics().getAverage(averageVfId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock2D *> freeSurfaceArgs,
                                   Box2D domain) {
  FS_AverageVolumeFractionFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getAverageVolumeFraction();
}

template <typename T, template <typename U> class Descriptor>
CountFreeSurfaceElementsFunctional2D<
    T, Descriptor>::CountFreeSurfaceElementsFunctional2D(plint flagToLookFor_)
    : numCellsId(this->getStatistics().subscribeIntSum()),
      flagToLookFor(flagToLookFor_) {}

template <typename T, template <typename U> class Descriptor>
void CountFreeSurfaceElementsFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      int materialIndex = param.flag(iX, iY);
      if (materialIndex == flagToLookFor) {  // Fluid Cell
        statistics.gatherIntSum(numCellsId, 1);
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
CountFreeSurfaceElementsFunctional2D<T, Descriptor>
    *CountFreeSurfaceElementsFunctional2D<T, Descriptor>::clone() const {
  return new CountFreeSurfaceElementsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint CountFreeSurfaceElementsFunctional2D<T, Descriptor>::getNumElements()
    const {
  return this->getStatistics().getIntSum(numCellsId);
}

template <typename T, template <typename U> class Descriptor>
plint countFreeSurfaceElements(std::vector<MultiBlock2D *> freeSurfaceArgs,
                               plint flagToLookFor, Box2D domain) {
  CountFreeSurfaceElementsFunctional2D<T, Descriptor> functional(flagToLookFor);
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getNumElements();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMomentumFunctional2D<T, Descriptor>::FS_AverageMomentumFunctional2D()
    : averageMomentumId(this->getStatistics().subscribeAverage(),
                        this->getStatistics().subscribeAverage(),
                        this->getStatistics().subscribeAverage()) {}

template <typename T, template <typename U> class Descriptor>
void FS_AverageMomentumFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.cell(iX, iY).getDynamics().hasMoments()) {
        Array<T, Descriptor<T>::d> j = param.getMomentum(iX, iY);
        statistics.gatherAverage(averageMomentumId[0], j[0]);
        statistics.gatherAverage(averageMomentumId[1], j[1]);
        statistics.incrementStats();
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageMomentumFunctional2D<T, Descriptor>
    *FS_AverageMomentumFunctional2D<T, Descriptor>::clone() const {
  return new FS_AverageMomentumFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
Array<T, 2> FS_AverageMomentumFunctional2D<T, Descriptor>::getAverageMomentum()
    const {
  return Array<T, 2>(this->getStatistics().getAverage(averageMomentumId[0]),
                     this->getStatistics().getAverage(averageMomentumId[1]));
}

template <typename T, template <typename U> class Descriptor>
Array<T, 2> freeSurfaceAverageMomentum(
    std::vector<MultiBlock2D *> freeSurfaceArgs, Box2D domain) {
  FS_AverageMomentumFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getAverageMomentum();
}

template <typename T, template <typename U> class Descriptor>
FS_AverageHeightFunctional2D<T, Descriptor>::FS_AverageHeightFunctional2D()
    : averageHeightId(this->getStatistics().subscribeAverage()) {}

template <typename T, template <typename U> class Descriptor>
void FS_AverageHeightFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();
  Dot2D absOffset = param.absOffset();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.cell(iX, iY).getDynamics().hasMoments()) {
        T localHeight = T(0);
        if (util::greaterEqual(param.volumeFraction(iX, iY), (T)1)) {
          localHeight = T(absOffset.y);  // Mehdi: z was returned, not sure if
                                         // returning y is correct?
        }
        statistics.gatherAverage(averageHeightId, localHeight);
        statistics.incrementStats();
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FS_AverageHeightFunctional2D<T, Descriptor>
    *FS_AverageHeightFunctional2D<T, Descriptor>::clone() const {
  return new FS_AverageHeightFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
T FS_AverageHeightFunctional2D<T, Descriptor>::getAverageHeight() const {
  return this->getStatistics().getAverage(averageHeightId);
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock2D *> freeSurfaceArgs,
                           Box2D domain) {
  FS_AverageHeightFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  return functional.getAverageHeight();
}

template <typename T, template <typename U> class Descriptor>
GetWaterLevelAtxyFunctional2D<T, Descriptor>::GetWaterLevelAtxyFunctional2D()
    : numFluidOccupiedCellId(this->getStatistics().subscribeIntSum()) {}

template <typename T, template <typename U> class Descriptor>
void GetWaterLevelAtxyFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.cell(iX, iY).getDynamics().hasMoments()) {
        if (param.volumeFraction(iX, iY) >= 0.5) {
          statistics.gatherIntSum(numFluidOccupiedCellId, 1);
        }
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
GetWaterLevelAtxyFunctional2D<T, Descriptor>
    *GetWaterLevelAtxyFunctional2D<T, Descriptor>::clone() const {
  return new GetWaterLevelAtxyFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint GetWaterLevelAtxyFunctional2D<T, Descriptor>::getNumFluidCellsAtXY()
    const {
  return this->getStatistics().getIntSum(numFluidOccupiedCellId);
}

template <typename T, template <typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock2D *> freeSurfaceArgs, plint N,
                       Box2D domain) {
  GetWaterLevelAtxyFunctional2D<T, Descriptor> functional;
  applyProcessingFunctional(functional, domain, freeSurfaceArgs);
  plint length_domain =
      domain.x1 - domain.x0;  // number of cell along y direction
  if (length_domain == 0) length_domain = 1;
  plint width_domain =
      domain.y1 - domain.y0;  // number of cell along y direction
  if (width_domain == 0) width_domain = 1;
  T heightAtXY =
      functional.getNumFluidCellsAtXY() / (T(N) * length_domain * width_domain);
  return heightAtXY;
}

template <typename T>
T freeSurfaceComputePorosity(MultiScalarField2D<int> &flag, Box2D domain) {
  plint numWallCells = freeSurfaceCountWallCells(flag, domain);
  T porosity = (T)1 - (T)numWallCells / (T)domain.nCells();
  return porosity;
}

template <typename T>
FreeSurfaceComputeFluidVolume2D<T>::FreeSurfaceComputeFluidVolume2D()
    : sumScalarId(this->getStatistics().subscribeSum()) {}

template <typename T>
void FreeSurfaceComputeFluidVolume2D<T>::process(
    Box2D domain, ScalarField2D<T> &volumeFraction, ScalarField2D<int> &flag) {
  Dot2D offset = computeRelativeDisplacement(volumeFraction, flag);
  BlockStatistics &statistics = this->getStatistics();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      int localFlag = flag.get(iX + offset.x, iY + offset.y);
      if (freeSurfaceFlag2D::isFullWet(localFlag)) {
        statistics.gatherSum(sumScalarId, (double)1.0);
      } else if (localFlag == freeSurfaceFlag2D::interface) {
        statistics.gatherSum(sumScalarId, (double)volumeFraction.get(iX, iY));
      }
    }
  }
}

template <typename T>
FreeSurfaceComputeFluidVolume2D<T> *FreeSurfaceComputeFluidVolume2D<T>::clone()
    const {
  return new FreeSurfaceComputeFluidVolume2D<T>(*this);
}

template <typename T>
T FreeSurfaceComputeFluidVolume2D<T>::getFluidVolume() const {
  double fluidVolume = this->getStatistics().getSum(sumScalarId);
  // The sum is internally computed on floating-point values. If T is
  //   integer, the value must be rounded at the end.
  if (std::numeric_limits<T>::is_integer) {
    return (T)util::roundToInt(fluidVolume);
  }
  return (T)fluidVolume;
}

template <typename T>
T freeSurfaceComputeFluidVolume(MultiScalarField2D<T> &volumeFraction,
                                MultiScalarField2D<int> &flag, Box2D domain) {
  FreeSurfaceComputeFluidVolume2D<T> functional;
  applyProcessingFunctional(functional, domain, volumeFraction, flag);
  return functional.getFluidVolume();
}

template <typename T>
MaskedFreeSurfaceComputeFluidVolume2D<
    T>::MaskedFreeSurfaceComputeFluidVolume2D()
    : sumScalarId(this->getStatistics().subscribeSum()) {}

template <typename T>
void MaskedFreeSurfaceComputeFluidVolume2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> fields) {
  ScalarField2D<T> &volumeFraction =
      *dynamic_cast<ScalarField2D<T> *>(fields[0]);
  ScalarField2D<int> &flag = *dynamic_cast<ScalarField2D<int> *>(fields[1]);
  ScalarField2D<int> &mask = *dynamic_cast<ScalarField2D<int> *>(fields[2]);
  Dot2D offsetFlag = computeRelativeDisplacement(volumeFraction, flag);
  Dot2D offsetMask = computeRelativeDisplacement(volumeFraction, mask);
  BlockStatistics &statistics = this->getStatistics();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (mask.get(iX + offsetMask.x, iY + offsetMask.y)) {
        int localFlag = flag.get(iX + offsetFlag.x, iY + offsetFlag.y);
        if (freeSurfaceFlag2D::isFullWet(localFlag)) {
          statistics.gatherSum(sumScalarId, (double)1.0);
        } else if (localFlag == freeSurfaceFlag2D::interface) {
          statistics.gatherSum(sumScalarId, (double)volumeFraction.get(iX, iY));
        }
      }
    }
  }
}

template <typename T>
MaskedFreeSurfaceComputeFluidVolume2D<T>
    *MaskedFreeSurfaceComputeFluidVolume2D<T>::clone() const {
  return new MaskedFreeSurfaceComputeFluidVolume2D<T>(*this);
}

template <typename T>
T MaskedFreeSurfaceComputeFluidVolume2D<T>::getFluidVolume() const {
  double fluidVolume = this->getStatistics().getSum(sumScalarId);
  // The sum is internally computed on floating-point values. If T is
  //   integer, the value must be rounded at the end.
  if (std::numeric_limits<T>::is_integer) {
    return (T)util::roundToInt(fluidVolume);
  }
  return (T)fluidVolume;
}

template <typename T>
T freeSurfaceComputeFluidVolume(MultiScalarField2D<T> &volumeFraction,
                                MultiScalarField2D<int> &flag,
                                MultiScalarField2D<int> &mask) {
  MaskedFreeSurfaceComputeFluidVolume2D<T> functional;
  std::vector<MultiBlock2D *> args;
  args.push_back(&volumeFraction);
  args.push_back(&flag);
  args.push_back(&mask);
  applyProcessingFunctional(functional, volumeFraction.getBoundingBox(), args);
  return functional.getFluidVolume();
}

template <typename T>
T freeSurfaceComputeSaturation(T porosity,
                               MultiScalarField2D<T> &volumeFraction,
                               MultiScalarField2D<int> &flag, Box2D domain) {
  T fluidVolume = freeSurfaceComputeFluidVolume(volumeFraction, flag, domain);
  T totalVolume = domain.nCells();
  T saturation = fluidVolume / (porosity * totalVolume);
  return saturation;
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceBoxSumForcedEnergyFunctional2D<
    T, Descriptor>::FreeSurfaceBoxSumForcedEnergyFunctional2D()
    : sumEnergyId(this->getStatistics().subscribeSum()),
      sumCellsId(this->getStatistics().subscribeSum()) {}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceBoxSumForcedEnergyFunctional2D<
    T, Descriptor>::processGenericBlocks(Box2D domain,
                                         std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 3);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  TensorField2D<T, Descriptor<T>::d> *force =
      dynamic_cast<TensorField2D<T, Descriptor<T>::d> *>(blocks[1]);
  PLB_ASSERT(force);
  ScalarField2D<int> *flag = dynamic_cast<ScalarField2D<int> *>(blocks[2]);
  PLB_ASSERT(flag);

  Dot2D ofsForce = computeRelativeDisplacement(*lattice, *force);
  Dot2D ofsFlag = computeRelativeDisplacement(*lattice, *flag);

  BlockStatistics &statistics = this->getStatistics();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (freeSurfaceFlag2D::isWet(flag->get(iX + ofsFlag.x, iY + ofsFlag.y))) {
        Array<T, Descriptor<T>::d> velocity;
        lattice->get(iX, iY).computeVelocity(velocity);
        Array<T, Descriptor<T>::d> const &f =
            force->get(iX + ofsForce.x, iY + ofsForce.y);
        velocity[0] += (T)0.5 * f[0];
        velocity[1] += (T)0.5 * f[1];
        T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
        statistics.gatherSum(sumEnergyId, uNormSqr);
        statistics.gatherSum(sumCellsId, (T)1);
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor>
    *FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor>::clone() const {
  return new FreeSurfaceBoxSumForcedEnergyFunctional2D(*this);
}

template <typename T, template <typename U> class Descriptor>
T FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor>::getAverageEnergy()
    const {
  T sumEnergy = this->getStatistics().getSum(sumEnergyId) / (T)2;
  T sumCells = this->getStatistics().getSum(sumCellsId);
  T averageEnergy = (T)0;
  if (!util::isZero(sumCells)) {
    averageEnergy = sumEnergy / sumCells;
  }
  return averageEnergy;
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;  // Lattice.
  modified[1] = modif::nothing;  // Force.
  modified[2] = modif::nothing;  // Flag.
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain) {
  FreeSurfaceBoxSumForcedEnergyFunctional2D<T, Descriptor> functional;
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&force);
  args.push_back(&flag);
  applyProcessingFunctional(functional, domain, args);
  return functional.getAverageEnergy();
}

template <typename T, template <typename U> class Descriptor>
T freeSurfaceComputeAverageForcedEnergy(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag) {
  return freeSurfaceComputeAverageForcedEnergy(lattice, force, flag,
                                               lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, MultiScalarField2D<T> &velocityNorm,
    Box2D domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&force);
  args.push_back(&velocityNorm);
  applyProcessingFunctional(
      new BoxForcedVelocityNormFunctional2D<T, Descriptor>, domain, args);

  setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag2D::empty, domain,
                   (T)0);
  setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag2D::protectEmpty,
                   domain, (T)0);
  setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag2D::wall, domain,
                   (T)0);
  setToConstant<T>(velocityNorm, flag, (int)freeSurfaceFlag2D::slipWall, domain,
                   (T)0);
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain) {
  std::auto_ptr<MultiScalarField2D<T>> velocityNorm =
      generateMultiScalarField<T>(lattice, domain);

  // The domain needs to be extended to the outer envelopes, for the following
  // reason. Imagine that the domain is smaller than the bounding-box. Given
  // that the BoxForcedVelocityNormFunctional2D() acts on both bulk and
  // envelope, you would expect the envelope layer around the domain, on the
  // velocityNorm multi-block, to be assigned some proper values too. By
  // default, this is however not what happens, because the physical space
  // occupied by these envelopes does not intersect with the domain "domain". We
  // work around this issue by extending the domain. There's no problem if the
  // enlarged domain gets beyond the actual extent of the lattice, because
  // Palabos handles these situations properly.

  freeSurfaceComputeForcedVelocityNorm(
      lattice, force, flag, *velocityNorm,
      domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
  return velocityNorm;
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityNorm(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag) {
  return freeSurfaceComputeForcedVelocityNorm(lattice, force, flag,
                                              lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, MultiScalarField2D<T> &velocityComponent,
    Box2D domain, plint iComponent) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&force);
  args.push_back(&velocityComponent);
  applyProcessingFunctional(
      new BoxForcedVelocityComponentFunctional2D<T, Descriptor>(iComponent),
      domain, args);

  setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag2D::empty,
                   domain, (T)0);
  setToConstant<T>(velocityComponent, flag,
                   (int)freeSurfaceFlag2D::protectEmpty, domain, (T)0);
  setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag2D::wall,
                   domain, (T)0);
  setToConstant<T>(velocityComponent, flag, (int)freeSurfaceFlag2D::slipWall,
                   domain, (T)0);
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, Box2D domain, plint iComponent) {
  std::auto_ptr<MultiScalarField2D<T>> velocityComponent =
      generateMultiScalarField<T>(lattice, domain);

  // The domain needs to be extended to the outer envelopes, for the following
  // reason. Imagine that the domain is smaller than the bounding-box. Given
  // that the BoxForcedVelocityComponentFunctional2D() acts on both bulk and
  // envelope, you would expect the envelope layer around the domain, on the
  // velocityComponent multi-block, to be assigned some proper values too. By
  // default, this is however not what happens, because the physical space
  // occupied by these envelopes does not intersect with the domain "domain". We
  // work around this issue by extending the domain. There's no problem if the
  // enlarged domain gets beyond the actual extent of the lattice, because
  // Palabos handles these situations properly.

  freeSurfaceComputeForcedVelocityComponent(
      lattice, force, flag, *velocityComponent,
      domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()),
      iComponent);
  return velocityComponent;
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T>> freeSurfaceComputeForcedVelocityComponent(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag, plint iComponent) {
  return freeSurfaceComputeForcedVelocityComponent(
      lattice, force, flag, lattice.getBoundingBox(), iComponent);
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeForcedVelocity(
    MultiBlockLattice2D<T, Descriptor> &lattice,
    MultiTensorField2D<T, Descriptor<T>::d> &force,
    MultiScalarField2D<int> &flag,
    MultiTensorField2D<T, Descriptor<T>::d> &velocity, Box2D domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&force);
  args.push_back(&velocity);
  applyProcessingFunctional(new BoxForcedVelocityFunctional2D<T, Descriptor>,
                            domain, args);

  setToConstant<T, Descriptor<T>::d>(velocity, flag,
                                     (int)freeSurfaceFlag2D::empty, domain,
                                     Array<T, Descriptor<T>::d>::zero());
  setToConstant<T, Descriptor<T>::d>(
      velocity, flag, (int)freeSurfaceFlag2D::protectEmpty, domain,
      Array<T, Descriptor<T>::d>::zero());
  setToConstant<T, Descriptor<T>::d>(velocity, flag,
                                     (int)freeSurfaceFlag2D::wall, domain,
                                     Array<T, Descriptor<T>::d>::zero());
  setToConstant<T, Descriptor<T>::d>(velocity, flag,
                                     (int)freeSurfaceFlag2D::slipWall, domain,
                                     Array<T, Descriptor<T>::d>::zero());
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T, Descriptor<T>::d>>
freeSurfaceComputeForcedVelocity(MultiBlockLattice2D<T, Descriptor> &lattice,
                                 MultiTensorField2D<T, Descriptor<T>::d> &force,
                                 MultiScalarField2D<int> &flag, Box2D domain) {
  std::auto_ptr<MultiTensorField2D<T, Descriptor<T>::d>> velocity =
      generateMultiTensorField<T, Descriptor<T>::d>(lattice, domain);

  // The domain needs to be extended to the outer envelopes, for the following
  // reason. Imagine that the domain is smaller than the bounding-box. Given
  // that the BoxForcedVelocityFunctional2D() acts on both bulk and envelope,
  // you would expect the envelope layer around the domain, on the velocity
  // multi-block, to be assigned some proper values too. By default, this is
  // however not what happens, because the physical space occupied by these
  // envelopes does not intersect with the domain "domain". We work around this
  // issue by extending the domain. There's no problem if the enlarged domain
  // gets beyond the actual extent of the lattice, because Palabos handles these
  // situations properly.

  freeSurfaceComputeForcedVelocity(
      lattice, force, flag, *velocity,
      domain.enlarge(lattice.getMultiBlockManagement().getEnvelopeWidth()));
  return velocity;
}

template <typename T, template <typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T, Descriptor<T>::d>>
freeSurfaceComputeForcedVelocity(MultiBlockLattice2D<T, Descriptor> &lattice,
                                 MultiTensorField2D<T, Descriptor<T>::d> &force,
                                 MultiScalarField2D<int> &flag) {
  return freeSurfaceComputeForcedVelocity(lattice, force, flag,
                                          lattice.getBoundingBox());
}

}  // namespace lbfoam

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_2D_HH
