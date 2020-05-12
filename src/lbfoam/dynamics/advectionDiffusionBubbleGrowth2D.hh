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


#ifndef ADVECTION_DIFFUSION_BUBBLE_GROWTH_2D_HH
#define ADVECTION_DIFFUSION_BUBBLE_GROWTH_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "core/dynamics.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/dynamics/advectionDiffusionBubbleGrowth2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"

namespace plb {
namespace lbfoam {

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
LatticeToPassiveAdvDiff2D<
    T, FluidDescriptor, ScalarDescriptor>::LatticeToPassiveAdvDiff2D(T scaling_)
    : scaling(scaling_) {}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>::process(
    Box2D domain, BlockLattice2D<T, FluidDescriptor> &fluid,
    BlockLattice2D<T, ScalarDescriptor> &scalar) {
  Dot2D offset = computeRelativeDisplacement(fluid, scalar);
  const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T *u = scalar.get(iX + offset.x, iY + offset.y).getExternal(velOffset);
      Array<T, 2> velocity;
      fluid.get(iX, iY).computeVelocity(velocity);
      velocity *= scaling;
      velocity.to_cArray(u);
    }
  }
}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor> *
LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>::clone() const {
  return new LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>(
      *this);
}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT LatticeToPassiveAdvDiff2D<
    T, FluidDescriptor, ScalarDescriptor>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;
  modified[1] = modif::staticVariables;
}

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void latticeToPassiveAdvDiff(MultiBlockLattice2D<T, FluidDescriptor> &fluid,
                             MultiBlockLattice2D<T, ScalarDescriptor> &scalar,
                             Box2D domain) {
  applyProcessingFunctional(
      new LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>(),
      domain, fluid, scalar);
}

/* Mehdi: Coupling between advection-diffusion and free-surface. */

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    AdvectionDiffusionFreeSurfaceCoupling2D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t)
    : conductiveDynamics(conductiveDynamics_),
      adiabaticDynamics(adiabaticDynamics_),
      iniVal(iniVal_) {
  PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
  PLB_ASSERT(!util::isZero(Pr_t));
  C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    ~AdvectionDiffusionFreeSurfaceCoupling2D() {
  delete conductiveDynamics;
  delete adiabaticDynamics;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    AdvectionDiffusionFreeSurfaceCoupling2D(
        AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                FS_Descriptor> const &rhs)
    : conductiveDynamics(rhs.conductiveDynamics->clone()),
      adiabaticDynamics(rhs.adiabaticDynamics->clone()),
      iniVal(rhs.iniVal),
      C(rhs.C) {}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor> &
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
operator=(AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                  FS_Descriptor> const &rhs) {
  AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>(rhs)
      .swap(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    swap(AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                 FS_Descriptor> &rhs) {
  std::swap(conductiveDynamics, rhs.conductiveDynamics);
  std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
  std::swap(iniVal, rhs.iniVal);
  std::swap(C, rhs.C);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  using namespace freeSurfaceFlag2D;

  PLB_ASSERT(blocks.size() == 3);
  BlockLattice2D<T, AD_Descriptor> *tempLattice =
      dynamic_cast<BlockLattice2D<T, AD_Descriptor> *>(blocks[0]);
  BlockLattice2D<T, FS_Descriptor> *fluidLattice =
      dynamic_cast<BlockLattice2D<T, FS_Descriptor> *>(blocks[1]);
  ScalarField2D<int> *flags = dynamic_cast<ScalarField2D<int> *>(blocks[2]);
  PLB_ASSERT(tempLattice);
  PLB_ASSERT(fluidLattice);
  PLB_ASSERT(flags);

  plint conductiveId = conductiveDynamics->getId();
  plint adiabaticId = adiabaticDynamics->getId();
  BlockLattice2D<T, AD_Descriptor> tempBak(*tempLattice);

  Dot2D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
  Dot2D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);

  const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
  Array<T, 2> zero;
  zero.resetToZero();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y);
      bool tempIsAdiabatic =
          tempBak.get(iX, iY).getDynamics().getId() == adiabaticId;

      if (flag == fluid && tempIsAdiabatic) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          cell[iPop] = T();
        }
        plint numNeighbors = 0;
        plint d = 1;
        for (plint dx = -d; dx <= d; dx++) {
          plint i = iX + dx;
          for (plint dy = -d; dy <= d; dy++) {
            plint j = iY + dy;
            Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j);
            // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y);
            // if(!(dx==0 && dy==0) &&
            // (nextCell.getDynamics().getId()==conductiveId) &&
            //        (nextFlag!=wall)) {
            if (!(dx == 0 && dy == 0) &&
                (nextCell.getDynamics().getId() == conductiveId)) {
              for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                cell[iPop] += nextCell[iPop];
              }
              ++numNeighbors;
            }
          }
        }
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          if (numNeighbors == 0) {
            cell[iPop] = AD_Descriptor<T>::t[iPop];
          } else {
            cell[iPop] /= numNeighbors;
          }
        }
      }

      if (flag == interface) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        iniCellAtEquilibrium(cell, iniVal, zero);
        tempLattice->attributeDynamics(iX, iY, adiabaticDynamics->clone());
      }

      // Set the advection velocity.

      Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY);
      Cell<T, FS_Descriptor> &fluidCell =
          fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y);

      T *u = tempCell.getExternal(velOffset);
      if (isFullWet(flag)) {
        Array<T, 2> velocity;
        fluidCell.computeVelocity(velocity);
        velocity.to_cArray(u);
      } else {
        zero.to_cArray(u);
      }

      // Set the relaxation parameter for the advection-diffusion.

      if (tempCell.getDynamics().getId() != adiabaticId) {
        T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, fluidCell);
        if (!util::isZero(fluidOmega)) {
          T fluidOmega0 = fluidCell.getDynamics().getOmega();
          T tempOmega0 = tempCell.getDynamics().getOmega();
          T tempOmega = (T)1 / ((T)1 / tempOmega0 +
                                C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
          tempCell.getDynamics().setOmega(tempOmega);
        }
      }
    }
  }
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>
    *AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                             FS_Descriptor>::clone() const {
  return new AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                     FS_Descriptor>(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::dataStructure;  // Temperature
  modified[1] = modif::nothing;        // FreeSurface Fluid
  modified[2] = modif::nothing;        // FreeSurface Flags
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t,
        int maskValue_)
    : conductiveDynamics(conductiveDynamics_),
      adiabaticDynamics(adiabaticDynamics_),
      iniVal(iniVal_),
      maskValue(maskValue_) {
  PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
  PLB_ASSERT(!util::isZero(Pr_t));
  C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    ~MaskedAdvectionDiffusionFreeSurfaceCoupling2D() {
  delete conductiveDynamics;
  delete adiabaticDynamics;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
    MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
        MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                      FS_Descriptor> const &rhs)
    : conductiveDynamics(rhs.conductiveDynamics->clone()),
      adiabaticDynamics(rhs.adiabaticDynamics->clone()),
      iniVal(rhs.iniVal),
      C(rhs.C),
      maskValue(rhs.maskValue) {}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor> &
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>::
operator=(MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
          T, AD_Descriptor, FS_Descriptor> const &rhs) {
  MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                FS_Descriptor>(rhs)
      .swap(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                   FS_Descriptor>::
    swap(MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                       FS_Descriptor> &rhs) {
  std::swap(conductiveDynamics, rhs.conductiveDynamics);
  std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
  std::swap(iniVal, rhs.iniVal);
  std::swap(C, rhs.C);
  std::swap(maskValue, rhs.maskValue);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
    T, AD_Descriptor,
    FS_Descriptor>::processGenericBlocks(Box2D domain,
                                         std::vector<AtomicBlock2D *> blocks) {
  using namespace freeSurfaceFlag2D;

  PLB_ASSERT(blocks.size() == 4);
  BlockLattice2D<T, AD_Descriptor> *tempLattice =
      dynamic_cast<BlockLattice2D<T, AD_Descriptor> *>(blocks[0]);
  BlockLattice2D<T, FS_Descriptor> *fluidLattice =
      dynamic_cast<BlockLattice2D<T, FS_Descriptor> *>(blocks[1]);
  ScalarField2D<int> *flags = dynamic_cast<ScalarField2D<int> *>(blocks[2]);
  ScalarField2D<int> *mask = dynamic_cast<ScalarField2D<int> *>(blocks[3]);
  PLB_ASSERT(tempLattice);
  PLB_ASSERT(fluidLattice);
  PLB_ASSERT(flags);
  PLB_ASSERT(mask);

  plint conductiveId = conductiveDynamics->getId();
  plint adiabaticId = adiabaticDynamics->getId();
  BlockLattice2D<T, AD_Descriptor> tempBak(*tempLattice);

  Dot2D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
  Dot2D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);
  Dot2D ofsMask = computeRelativeDisplacement(*tempLattice, *mask);

  const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
  Array<T, 2> zero;
  zero.resetToZero();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (mask->get(iX + ofsMask.x, iY + ofsMask.y) != maskValue) {
        continue;
      }

      plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y);
      bool tempIsAdiabatic =
          tempBak.get(iX, iY).getDynamics().getId() == adiabaticId;

      if (flag == fluid && tempIsAdiabatic) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          cell[iPop] = T();
        }
        plint numNeighbors = 0;
        plint d = 1;
        for (plint dx = -d; dx <= d; dx++) {
          plint i = iX + dx;
          for (plint dy = -d; dy <= d; dy++) {
            plint j = iY + dy;
            Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j);
            // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y);
            // if(!(dx==0 && dy==0) &&
            // (nextCell.getDynamics().getId()==conductiveId) &&
            //        (nextFlag!=wall)) {
            // TODO: Should the mask values be checked for the neighbors?
            if (!(dx == 0 && dy == 0) &&
                (nextCell.getDynamics().getId() == conductiveId)) {
              for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                cell[iPop] += nextCell[iPop];
              }
              ++numNeighbors;
            }
          }
        }
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          if (numNeighbors == 0) {
            cell[iPop] = AD_Descriptor<T>::t[iPop];
          } else {
            cell[iPop] /= numNeighbors;
          }
        }
      }

      if (flag == interface) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        iniCellAtEquilibrium(cell, iniVal, zero);
        tempLattice->attributeDynamics(iX, iY, adiabaticDynamics->clone());
      }

      // Set the advection velocity.

      Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY);
      Cell<T, FS_Descriptor> &fluidCell =
          fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y);

      T *u = tempCell.getExternal(velOffset);
      if (isFullWet(flag)) {
        Array<T, 2> velocity;
        fluidCell.computeVelocity(velocity);
        velocity.to_cArray(u);
      } else {
        zero.to_cArray(u);
      }

      // Set the relaxation parameter for the advection-diffusion.

      if (tempCell.getDynamics().getId() != adiabaticId) {
        T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, fluidCell);
        if (!util::isZero(fluidOmega)) {
          T fluidOmega0 = fluidCell.getDynamics().getOmega();
          T tempOmega0 = tempCell.getDynamics().getOmega();
          T tempOmega = (T)1 / ((T)1 / tempOmega0 +
                                C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
          tempCell.getDynamics().setOmega(tempOmega);
        }
      }
    }
  }
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor> *
MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                              FS_Descriptor>::clone() const {
  return new MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                           FS_Descriptor>(
      *this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                   FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::dataStructure;  // Temperature
  modified[1] = modif::nothing;        // FreeSurface Fluid
  modified[2] = modif::nothing;        // FreeSurface Flags
  modified[3] = modif::nothing;        // Mask
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                FS_Descriptor>::
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t,
        int maskValue_)
    : conductiveDynamics(conductiveDynamics_),
      adiabaticDynamics(adiabaticDynamics_),
      iniVal(iniVal_),
      maskValue(maskValue_) {
  PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
  PLB_ASSERT(!util::isZero(Pr_t));
  C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
    T, AD_Descriptor,
    FS_Descriptor>::~N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D() {
  delete conductiveDynamics;
  delete adiabaticDynamics;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                FS_Descriptor>::
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
            T, AD_Descriptor, FS_Descriptor> const &rhs)
    : conductiveDynamics(rhs.conductiveDynamics->clone()),
      adiabaticDynamics(rhs.adiabaticDynamics->clone()),
      iniVal(rhs.iniVal),
      C(rhs.C),
      maskValue(rhs.maskValue) {}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>
    &N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                     FS_Descriptor>::
    operator=(N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
              T, AD_Descriptor, FS_Descriptor> const &rhs) {
  N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                  FS_Descriptor>(rhs)
      .swap(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                     FS_Descriptor>::
    swap(N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                         FS_Descriptor> &rhs) {
  std::swap(conductiveDynamics, rhs.conductiveDynamics);
  std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
  std::swap(iniVal, rhs.iniVal);
  std::swap(C, rhs.C);
  std::swap(maskValue, rhs.maskValue);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
    T, AD_Descriptor,
    FS_Descriptor>::processGenericBlocks(Box2D domain,
                                         std::vector<AtomicBlock2D *> blocks) {
  using namespace freeSurfaceFlag2D;

  PLB_ASSERT(blocks.size() == 4);
  BlockLattice2D<T, AD_Descriptor> *tempLattice =
      dynamic_cast<BlockLattice2D<T, AD_Descriptor> *>(blocks[0]);
  BlockLattice2D<T, FS_Descriptor> *fluidLattice =
      dynamic_cast<BlockLattice2D<T, FS_Descriptor> *>(blocks[1]);
  ScalarField2D<int> *flags = dynamic_cast<ScalarField2D<int> *>(blocks[2]);
  NTensorField2D<int> *mask = dynamic_cast<NTensorField2D<int> *>(blocks[3]);
  PLB_ASSERT(tempLattice);
  PLB_ASSERT(fluidLattice);
  PLB_ASSERT(flags);
  PLB_ASSERT(mask);

  plint conductiveId = conductiveDynamics->getId();
  plint adiabaticId = adiabaticDynamics->getId();
  BlockLattice2D<T, AD_Descriptor> tempBak(*tempLattice);

  Dot2D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
  Dot2D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);
  Dot2D ofsMask = computeRelativeDisplacement(*tempLattice, *mask);

  const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
  Array<T, 2> zero;
  zero.resetToZero();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (*mask->get(iX + ofsMask.x, iY + ofsMask.y) != maskValue) {
        continue;
      }

      plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y);
      bool tempIsAdiabatic =
          tempBak.get(iX, iY).getDynamics().getId() == adiabaticId;

      if (flag == fluid && tempIsAdiabatic) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          cell[iPop] = T();
        }
        plint numNeighbors = 0;
        plint d = 1;
        for (plint dx = -d; dx <= d; dx++) {
          plint i = iX + dx;
          for (plint dy = -d; dy <= d; dy++) {
            plint j = iY + dy;

            Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j);
            // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y);
            // if(!(dx==0 && dy==0) &&
            // (nextCell.getDynamics().getId()==conductiveId) &&
            //        (nextFlag!=wall)) {
            // TODO: Should the mask values be checked for the neighbors?
            if (!(dx == 0 && dy == 0) &&
                (nextCell.getDynamics().getId() == conductiveId)) {
              for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                cell[iPop] += nextCell[iPop];
              }
              ++numNeighbors;
            }
          }
        }
        for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
          if (numNeighbors == 0) {
            cell[iPop] = AD_Descriptor<T>::t[iPop];
          } else {
            cell[iPop] /= numNeighbors;
          }
        }
      }

      if (flag == interface) {
        tempLattice->attributeDynamics(iX, iY, conductiveDynamics->clone());
        Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY);
        iniCellAtEquilibrium(cell, iniVal, zero);
        tempLattice->attributeDynamics(iX, iY, adiabaticDynamics->clone());
      }

      // Set the advection velocity.

      Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY);
      Cell<T, FS_Descriptor> &fluidCell =
          fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y);

      T *u = tempCell.getExternal(velOffset);
      if (isFullWet(flag)) {
        Array<T, 2> velocity;
        fluidCell.computeVelocity(velocity);
        velocity.to_cArray(u);
      } else {
        zero.to_cArray(u);
      }

      // Set the relaxation parameter for the advection-diffusion.

      if (tempCell.getDynamics().getId() != adiabaticId) {
        T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, fluidCell);
        if (!util::isZero(fluidOmega)) {
          T fluidOmega0 = fluidCell.getDynamics().getOmega();
          T tempOmega0 = tempCell.getDynamics().getOmega();
          T tempOmega = (T)1 / ((T)1 / tempOmega0 +
                                C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
          tempCell.getDynamics().setOmega(tempOmega);
        }
      }
    }
  }
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>
    *N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
        T, AD_Descriptor, FS_Descriptor>::clone() const {
  return new N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                             FS_Descriptor>(
      *this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                     FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::dataStructure;  // Temperature
  modified[1] = modif::nothing;        // FreeSurface Fluid
  modified[2] = modif::nothing;        // FreeSurface Flags
  modified[3] = modif::nothing;        // Mask
}

/* Mehdi: Coupling between advection-diffusion for growth and free-surface. */

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::GrowthCoupling2D(
    Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics_,
    Dynamics<T, AD_Descriptor> *emptyDynamics_, T Kh_,
    std::map<plint, BubbleInfo2D> const &bubbles_, bool surfaceDiffusion_)
    : advectionDiffusionDynamics(advectionDiffusionDynamics_),
      emptyDynamics(emptyDynamics_),
      Kh(Kh_),
      bubbles(bubbles_),
      surfaceDiffusion(surfaceDiffusion_)
{}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::~GrowthCoupling2D() {
  delete advectionDiffusionDynamics;
  delete emptyDynamics;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  using namespace freeSurfaceFlag2D;
  typedef typename InterfaceLists2D<T, AD_Descriptor>::Node Node;
  PLB_ASSERT(blocks.size() == 6);

  BlockLattice2D<T, AD_Descriptor> *gasLattice =
      dynamic_cast<BlockLattice2D<T, AD_Descriptor> *>(blocks[0]);
  BlockLattice2D<T, FS_Descriptor> *fluidLattice =
      dynamic_cast<BlockLattice2D<T, FS_Descriptor> *>(blocks[1]);
  ScalarField2D<int> *flags = dynamic_cast<ScalarField2D<int> *>(blocks[2]);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[3]);
  ScalarField2D<T> *outsideDensity =
      dynamic_cast<ScalarField2D<T> *>(blocks[4]);
  ScalarField2D<plint> *tags = dynamic_cast<ScalarField2D<plint> *>(blocks[5]);

  // plint EmptyId = emptyDynamics->getId();
  plint adID = advectionDiffusionDynamics->getId();

  BlockLattice2D<T, AD_Descriptor> gasBak(*gasLattice);

  PLB_ASSERT(gasLattice);
  PLB_ASSERT(fluidLattice);
  PLB_ASSERT(flags);
  PLB_ASSERT(j);
  PLB_ASSERT(tags);

  std::map<Node, Array<T, AD_Descriptor<T>::q>> neighborOppositePop;
  Dot2D ofsTags = computeRelativeDisplacement(*gasLattice, *tags);
  Dot2D ofsFluid = computeRelativeDisplacement(*gasLattice, *fluidLattice);
  Dot2D ofsFlag = computeRelativeDisplacement(*gasLattice, *flags);
  Dot2D relativeOffsetJ = computeRelativeDisplacement(*gasLattice, *j);
  Dot2D ofsOutsideDensity =
      computeRelativeDisplacement(*gasLattice, *outsideDensity);

  const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
  Array<T, 2> zero;
  zero.resetToZero();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      Cell<T, AD_Descriptor> &gasCell = gasLattice->get(iX, iY);
      Cell<T, FS_Descriptor> &fluidCell =
          fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y);
      plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y);
      T *u = gasCell.getExternal(velOffset);
      if (isFullWet(flag)) {
        gasLattice->attributeDynamics(iX, iY,
                                      advectionDiffusionDynamics->clone());

        Array<T, 2> velocity;
        fluidCell.computeVelocity(velocity);
        velocity.to_cArray(u);
      } else {
        gasLattice->attributeDynamics(iX, iY, emptyDynamics->clone());
        zero.to_cArray(u);
      }

      if (flag == interface) {
        bool needsModification = false;
        Array<T, AD_Descriptor<T>::q> savedPop;
        savedPop[0] = T() - 2.;

        for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop) {
          plint prevX = iX - AD_Descriptor<T>::c[iPop][0];
          plint prevY = iY - AD_Descriptor<T>::c[iPop][1];

          plint nextX = iX + AD_Descriptor<T>::c[iPop][0];
          plint nextY = iY + AD_Descriptor<T>::c[iPop][1];

          plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);

          plint prevFlag = flags->get(prevX + ofsFlag.x, prevY + ofsFlag.y);
          plint nextFlag = flags->get(nextX + ofsFlag.x, nextY + ofsFlag.y);

          if (isEmpty(prevFlag) || isAnyWall(prevFlag))
          // If the f_i[iPop] would be streamed from an empty cell
          {
            if (isEmpty(nextFlag) || isAnyWall(nextFlag)) {
              savedPop[iPop] = (T) + 1000.;
            } else {
              savedPop[iPop] = gasLattice->get(
                  iX, iY)[opp];  // So that g_\bar{i}^in=g(x,t+dt)
              needsModification = true;
            }
          } else {
            savedPop[iPop] = (T)-1.;
          }
        }

        if (surfaceDiffusion == false) {
          plint tag = tags->get(iX + ofsTags.x, iY + ofsTags.y);

          std::map<plint, BubbleInfo2D>::const_iterator it = bubbles.find(tag);

          if (it->second.isFrozen()) {
            needsModification = false;

            gasLattice->attributeDynamics(iX, iY,
                                          advectionDiffusionDynamics->clone());
            Cell<T, AD_Descriptor> &cell = gasLattice->get(iX, iY);
            for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
              cell[iPop] = T();
            }
            plint numNeighbors = 0;
            plint d = 1;
            for (plint dx = -d; dx <= d; dx++) {
              plint i = iX + dx;
              for (plint dy = -d; dy <= d; dy++) {
                plint j = iY + dy;
                Cell<T, AD_Descriptor> const &nextCell = gasBak.get(i, j);
                // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y);
                // if(!(dx==0 && dy==0) &&
                // (nextCell.getDynamics().getId()==conductiveId) &&
                //        (nextFlag!=wall)) {
                if (!(dx == 0 && dy == 0) &&
                    (nextCell.getDynamics().getId() == adID)) {
                  for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                    cell[iPop] += nextCell[iPop];
                  }
                  ++numNeighbors;
                }
              }
            }
            for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
              if (numNeighbors == 0) {
                cell[iPop] = AD_Descriptor<T>::t[iPop];
              } else {
                cell[iPop] /= numNeighbors;
              }
            }
          }
        }

        if (needsModification) {
          neighborOppositePop.insert(
              std::pair<Node, Array<T, AD_Descriptor<T>::q>>(Node(iX, iY),
                                                             savedPop));
        }

        gasLattice->attributeDynamics(iX, iY,
                                      advectionDiffusionDynamics->clone());
      }
    }
  }

  typename std::map<Node, Array<T, AD_Descriptor<T>::q>>::const_iterator nodes =
      neighborOppositePop.begin();
  for (; nodes != neighborOppositePop.end(); ++nodes) {
    Node node = nodes->first;
    plint iX = node[0];
    plint iY = node[1];
    Cell<T, AD_Descriptor> &gasCell = gasLattice->get(iX, iY);
    Array<T, AD_Descriptor<T>::q> neighborOppPop = nodes->second;
    for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop) {
      if (neighborOppPop[iPop] > (T)-1.) {
        // Velocity is simply taken from the previous time step.
        Array<T, 2> j_FS =
            j->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y);
        T jSqr = VectorTemplate<T, FS_Descriptor>::normSqr(j_FS);
        // Remember: the value of pressure on an interface node has been set in
        // FreeSurfaceMacroscopic2D, and is equal to the ambient pressure for a
        // single free-surface fluid, or in the case of a binary pressure, an
        // averaged value.
        T rho_FS_adj = outsideDensity->get(iX + ofsOutsideDensity.x,
                                           iY + ofsOutsideDensity.y) *
                       Kh * AD_Descriptor<T>::cs2;  // Apply Henry's law

        //                    pcout << "  rhoadj: " << rho_FS_adj << "
        //                    outsidedensity: " <<
        //                    outsideDensity->get(iX+ofsOutsideDensity.x,iY+ofsOutsideDensity.y)
        //                          << "  Kh: "  << Kh   << std::endl;

        T rhoBar_FS_adj = AD_Descriptor<T>::rhoBar(rho_FS_adj);
        T gEq_i = gasCell.computeEquilibrium(iPop, rhoBar_FS_adj, j_FS, jSqr);
        plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);
        T gEq_opp_i =
            gasCell.computeEquilibrium(opp, rhoBar_FS_adj, j_FS, jSqr);

        // Condition that the opposite iPops are unkown
        if (util::fpequal(neighborOppPop[iPop], (T) + 1000.) ||
            neighborOppPop[opp] > (T)1.) {
          gasCell[iPop] = gEq_i;
        } else {
          gasCell[iPop] = gEq_i + gEq_opp_i - neighborOppPop[iPop];
        }

        // pcout << "Compdensity:    " << gasCell.computeDensity() << std::endl;
      }
    }
  }
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::GrowthCoupling2D(
    GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> const &rhs)
    : advectionDiffusionDynamics(rhs.advectionDiffusionDynamics->clone()),
      emptyDynamics(rhs.emptyDynamics->clone()),
      Kh(rhs.Kh),
      bubbles(rhs.bubbles),
      surfaceDiffusion(rhs.surfaceDiffusion)
{}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>
    *GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::clone() const {
  return new GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::dataStructure;  // Gas
  modified[1] = modif::nothing;        // FreeSurface Fluid
  modified[2] = modif::nothing;        // FreeSurface Flags
  modified[3] = modif::nothing;        // FreeSurface j
  modified[4] = modif::nothing;        // OutsideDensity
  modified[5] = modif::nothing;        // tags
}

}  // namespace lbfoam
}  // namespace plb

#endif  // ADVECTION_DIFFUSION_2D_HH
