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

#ifndef ADVECTION_DIFFUSION_BUBBLE_GROWTH_2D_H
#define ADVECTION_DIFFUSION_BUBBLE_GROWTH_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"
#include "lbfoam/bubble/bubbleGrowth2D.h"

namespace plb {
namespace lbfoam {

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations: the fluid velocity is copied
 * to the advection-diffusion field, which is advected passively.
 */
template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
class LatticeToPassiveAdvDiff2D
    : public BoxProcessingFunctional2D_LL<T, FluidDescriptor, T,
                                          ScalarDescriptor> {
 public:
  LatticeToPassiveAdvDiff2D(T scaling_ = 1.);
  virtual void process(Box2D domain, BlockLattice2D<T, FluidDescriptor> &fluid,
                       BlockLattice2D<T, ScalarDescriptor> &scalar);
  virtual LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor>
      *clone() const;
  virtual BlockDomain::DomainT appliesTo() const;
  void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

 private:
  T scaling;
};

template <typename T, template <typename U1> class FluidDescriptor,
          template <typename U2> class ScalarDescriptor>
void latticeToPassiveAdvDiff(MultiBlockLattice2D<T, FluidDescriptor> &fluid,
                             MultiBlockLattice2D<T, ScalarDescriptor> &scalar,
                             Box2D domain);

/* Mehdi: Coupling between the 2D free-surface model and an advection-diffusion
 * equation. The "empty space" of the free-surface simulation is approximated by
 * an adiabatic medium in the advection-diffusion simulation.
 */
template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
class AdvectionDiffusionFreeSurfaceCoupling2D
    : public BoxProcessingFunctional2D {
 public:
  // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic
  // viscosity and the eddy diffusivity).
  AdvectionDiffusionFreeSurfaceCoupling2D(
      Dynamics<T, AD_Descriptor> *conductiveDynamics_,
      Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t);
  ~AdvectionDiffusionFreeSurfaceCoupling2D();
  AdvectionDiffusionFreeSurfaceCoupling2D(
      AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                              FS_Descriptor> const &rhs);
  AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor> &
  operator=(AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                    FS_Descriptor> const &rhs);
  void swap(AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                    FS_Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual AdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                  FS_Descriptor>
      *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  Dynamics<T, AD_Descriptor> *conductiveDynamics;
  Dynamics<T, AD_Descriptor> *adiabaticDynamics;
  T iniVal;
  T C;
};

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
class MaskedAdvectionDiffusionFreeSurfaceCoupling2D
    : public BoxProcessingFunctional2D {
 public:
  // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic
  // viscosity and the eddy diffusivity).
  MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
      Dynamics<T, AD_Descriptor> *conductiveDynamics_,
      Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t,
      int maskValue_);
  ~MaskedAdvectionDiffusionFreeSurfaceCoupling2D();
  MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
      MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                    FS_Descriptor> const &rhs);
  MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor, FS_Descriptor>
      &operator=(MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
                 T, AD_Descriptor, FS_Descriptor> const &rhs);
  void swap(MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                          FS_Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                        FS_Descriptor>
      *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  Dynamics<T, AD_Descriptor> *conductiveDynamics;
  Dynamics<T, AD_Descriptor> *adiabaticDynamics;
  T iniVal;
  T C;
  int maskValue;
};

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
class N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D
    : public BoxProcessingFunctional2D {
 public:
  // Pr_t is the turbulent Prandtl number (the ratio between the eddy kinematic
  // viscosity and the eddy diffusivity).
  N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
      Dynamics<T, AD_Descriptor> *conductiveDynamics_,
      Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t,
      int maskValue_);
  ~N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D();
  N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D(
      N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
          T, AD_Descriptor, FS_Descriptor> const &rhs);
  N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                  FS_Descriptor>
      &operator=(N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
                 T, AD_Descriptor, FS_Descriptor> const &rhs);
  void swap(N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<
            T, AD_Descriptor, FS_Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual N_MaskedAdvectionDiffusionFreeSurfaceCoupling2D<T, AD_Descriptor,
                                                          FS_Descriptor>
      *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  Dynamics<T, AD_Descriptor> *conductiveDynamics;
  Dynamics<T, AD_Descriptor> *adiabaticDynamics;
  T iniVal;
  T C;
  int maskValue;
};

/* Mehdi: Coupling between the free-surface model and an advection-diffusion
 * equation for growth. The "empty space" of the free-surface simulation is
 * approximated by an empty medium with empty dynamics.
 */
template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
class GrowthCoupling2D : public BoxProcessingFunctional2D {
 public:
  GrowthCoupling2D(Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics_,
                   Dynamics<T, AD_Descriptor> *emptyDynamics_, T Kh_,
                   std::map<plint, BubbleInfo2D> const &bubbles_,
                   bool surfaceDiffusion_);
  ~GrowthCoupling2D();
  GrowthCoupling2D(
      GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> const &rhs);
  GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> &operator=(
      GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> const &rhs);
  void swap(GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> &rhs);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> blocks);
  virtual GrowthCoupling2D<T, AD_Descriptor, FS_Descriptor> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics;
  Dynamics<T, AD_Descriptor> *emptyDynamics;
  T Kh;
  std::map<plint, BubbleInfo2D> bubbles;
  bool surfaceDiffusion;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // ADVECTION_DIFFUSION_2D_H
