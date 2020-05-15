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

#ifndef BODY_FORCE_2D_HH
#define BODY_FORCE_2D_HH

#include "lbfoam/models/bodyForce2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"

namespace plb {

namespace lbfoam {

/* ****************** AddConstForceToMomentum2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor>
AddConstForceToMomentum2D<T, Descriptor>::AddConstForceToMomentum2D(
    Array<T, 2> const &force_)
    : force(force_) {}

template <typename T, template <typename U> class Descriptor>
void AddConstForceToMomentum2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 3);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);

  Dot2D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      // Momentum correction.
      Cell<T, Descriptor> &cell = lattice->get(iX, iY);
      if (cell.getDynamics().hasMoments()) {
        T dynamicOmega = cell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, cell);
        T tau = 0.0;
        if (!util::isZero(dynamicOmega)) {
          tau = (T)1 / dynamicOmega;
        } else {
          tau = (T)1 / cell.getDynamics().getOmega();
        }
        T rho = Descriptor<T>::fullRho(rhoBar->get(iX + ofsR.x, iY + ofsR.y));
        j->get(iX + ofsJ.x, iY + ofsJ.y) += rho * tau * force;
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
AddConstForceToMomentum2D<T, Descriptor>
    *AddConstForceToMomentum2D<T, Descriptor>::clone() const {
  return new AddConstForceToMomentum2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AddConstForceToMomentum2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::staticVariables;  // j
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AddConstForceToMomentum2D<T, Descriptor>::appliesTo()
    const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void addConstForceToMomentum(MultiBlockLattice2D<T, Descriptor> &lattice,
                             MultiScalarField2D<T> &rhoBar,
                             MultiTensorField2D<T, 2> &j,
                             Array<T, 2> const &force, Box2D const &domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  applyProcessingFunctional(new AddConstForceToMomentum2D<T, Descriptor>(force),
                            domain, args);
}

/* ****************** AddCustomForceToMomentum2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
AddCustomForceToMomentum2D<
    T, Descriptor, ForceFunction>::AddCustomForceToMomentum2D(ForceFunction f_)
    : f(f_) {}

template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
void AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 3);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);

  Dot2D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);

  Dot2D location = lattice->getLocation();

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    plint x = iX + location.x;
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      plint y = iY + location.y;
      // Momentum correction.
      Cell<T, Descriptor> &cell = lattice->get(iX, iY);
      if (cell.getDynamics().hasMoments()) {
        T dynamicOmega = cell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, cell);
        T tau = 0.0;
        if (!util::isZero(dynamicOmega)) {
          tau = (T)1 / dynamicOmega;
        } else {
          tau = (T)1 / cell.getDynamics().getOmega();
        }
        T rho = Descriptor<T>::fullRho(rhoBar->get(iX + ofsR.x, iY + ofsR.y));
        Array<T, 2> force;
        f(x, y, force);
        j->get(iX + ofsJ.x, iY + ofsJ.y) += rho * tau * force;
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>
    *AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>::clone() const {
  return new AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>(*this);
}

template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
void AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::staticVariables;  // j
}

template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
BlockDomain::DomainT
AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T, template <class U> class Descriptor, class ForceFunction>
void addCustomForceToMomentum(MultiBlockLattice2D<T, Descriptor> &lattice,
                              MultiScalarField2D<T> &rhoBar,
                              MultiTensorField2D<T, 2> &j, ForceFunction f,
                              Box2D const &domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  applyProcessingFunctional(
      new AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>(f), domain,
      args);
}

/* ****************** AddForceToMomentum2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor>
void AddForceToMomentum2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 4);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);
  TensorField2D<T, 2> *force = dynamic_cast<TensorField2D<T, 2> *>(blocks[3]);
  PLB_ASSERT(force);

  Dot2D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);
  Dot2D ofsF = computeRelativeDisplacement(*lattice, *force);

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      // Momentum correction.
      Cell<T, Descriptor> &cell = lattice->get(iX, iY);
      if (cell.getDynamics().hasMoments()) {
        T dynamicOmega = cell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, cell);
        T tau = 0.0;
        if (!util::isZero(dynamicOmega)) {
          tau = (T)1 / dynamicOmega;
        } else {
          tau = (T)1 / cell.getDynamics().getOmega();
        }
        T rho = Descriptor<T>::fullRho(rhoBar->get(iX + ofsR.x, iY + ofsR.y));
        Array<T, 2> const &f = force->get(iX + ofsF.x, iY + ofsF.y);
        j->get(iX + ofsJ.x, iY + ofsJ.y) += rho * tau * f;
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
AddForceToMomentum2D<T, Descriptor>
    *AddForceToMomentum2D<T, Descriptor>::clone() const {
  return new AddForceToMomentum2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AddForceToMomentum2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::staticVariables;  // j
  modified[3] = modif::nothing;          // Force
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AddForceToMomentum2D<T, Descriptor>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void addForceToMomentum(MultiBlockLattice2D<T, Descriptor> &lattice,
                        MultiScalarField2D<T> &rhoBar,
                        MultiTensorField2D<T, 2> &j,
                        MultiTensorField2D<T, 2> &force, Box2D const &domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  args.push_back(&force);
  applyProcessingFunctional(new AddForceToMomentum2D<T, Descriptor>, domain,
                            args);
}

/* ****************** FreeSurfaceAddForceToMomentum2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddForceToMomentum2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 5);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);
  ScalarField2D<int> *flag = dynamic_cast<ScalarField2D<int> *>(blocks[3]);
  PLB_ASSERT(flag);
  TensorField2D<T, 2> *force = dynamic_cast<TensorField2D<T, 2> *>(blocks[4]);
  PLB_ASSERT(force);

  Dot2D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);
  Dot2D ofsFlag = computeRelativeDisplacement(*lattice, *flag);
  Dot2D ofsForce = computeRelativeDisplacement(*lattice, *force);

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      if (lbfoam::freeSurfaceFlag2D::isWet(
              flag->get(iX + ofsFlag.x, iY + ofsFlag.y))) {
        // Momentum correction.
        Cell<T, Descriptor> &cell = lattice->get(iX, iY);
        T dynamicOmega = cell.getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega, cell);
        T tau = 0.0;
        if (!util::isZero(dynamicOmega)) {
          tau = (T)1 / dynamicOmega;
        } else {
          tau = (T)1 / cell.getDynamics().getOmega();
        }
        T rho = Descriptor<T>::fullRho(
            rhoBar->get(iX + ofsRhoBar.x, iY + ofsRhoBar.y));
        Array<T, 2> const &f = force->get(iX + ofsForce.x, iY + ofsForce.y);
        j->get(iX + ofsJ.x, iY + ofsJ.y) += rho * tau * f;
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceAddForceToMomentum2D<T, Descriptor>
    *FreeSurfaceAddForceToMomentum2D<T, Descriptor>::clone() const {
  return new FreeSurfaceAddForceToMomentum2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddForceToMomentum2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::staticVariables;  // j
  modified[3] = modif::nothing;          // Flag
  modified[4] = modif::nothing;          // Force
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT FreeSurfaceAddForceToMomentum2D<T, Descriptor>::appliesTo()
    const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceAddForceToMomentum(MultiBlockLattice2D<T, Descriptor> &lattice,
                                   MultiScalarField2D<T> &rhoBar,
                                   MultiTensorField2D<T, 2> &j,
                                   MultiScalarField2D<int> &flag,
                                   MultiTensorField2D<T, 2> &force,
                                   Box2D const &domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  args.push_back(&flag);
  args.push_back(&force);
  applyProcessingFunctional(new FreeSurfaceAddForceToMomentum2D<T, Descriptor>,
                            domain, args);
}

/* ****************** ComputeRotatingFrameForce2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor>
ComputeRotatingFrameForce2D<T, Descriptor>::ComputeRotatingFrameForce2D(
    Array<T, 2> const &constantForce_, Array<T, 2> const &angularVelocity_,
    Array<T, 2> const &origin_, bool incompressibleModel_)
    : constantForce(constantForce_),
      angularVelocity(angularVelocity_),
      origin(origin_),
      incompressibleModel(incompressibleModel_) {}

template <typename T, template <typename U> class Descriptor>
void ComputeRotatingFrameForce2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 4);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);
  TensorField2D<T, 2> *force = dynamic_cast<TensorField2D<T, 2> *>(blocks[3]);
  PLB_ASSERT(force);

  Dot2D loc = lattice->getLocation();

  Dot2D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);
  Dot2D ofsForce = computeRelativeDisplacement(*lattice, *force);

  if (!util::isZero(norm(constantForce))) {
    T angularVelocityNorm = norm(angularVelocity);
    if (!util::isZero(angularVelocityNorm)) {
      Array<T, 2> normedAxis = angularVelocity / angularVelocityNorm;
      plint t = lattice->getTimeCounter().getTime();
      T theta = -angularVelocityNorm * t;
      constantForce = rotateAtOrigin(constantForce, normedAxis, theta);
    }
  }

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    T x = iX + loc.x - origin[0];
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      T y = iY + loc.y - origin[1];

      // Constant force.
      Array<T, 2> newForce(constantForce);

      // Coriolis and centripetal forces.
      Array<T, 2> r(x, y);
      Array<T, 2> velocity(j->get(iX + ofsJ.x, iY + ofsJ.y));
      if (!incompressibleModel) {
        T rho = Descriptor<T>::fullRho(
            rhoBar->get(iX + ofsRhoBar.x, iY + ofsRhoBar.y));
        velocity /= rho;
      }
      Array<T, 2> const &oldForce =
          force->get(iX + ofsForce.x, iY + ofsForce.y);
      velocity += (T)0.5 * oldForce;
      newForce +=
          -((T)2 * crossProduct(angularVelocity, velocity) +
            crossProduct(angularVelocity, crossProduct(angularVelocity, r)));

      force->get(iX + ofsForce.x, iY + ofsForce.y) = newForce;
    }
  }
}

template <typename T, template <typename U> class Descriptor>
ComputeRotatingFrameForce2D<T, Descriptor>
    *ComputeRotatingFrameForce2D<T, Descriptor>::clone() const {
  return new ComputeRotatingFrameForce2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ComputeRotatingFrameForce2D<T, Descriptor>::appliesTo()
    const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void ComputeRotatingFrameForce2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::nothing;          // j
  modified[3] = modif::staticVariables;  // Force
}

template <typename T, template <typename U> class Descriptor>
void computeRotatingFrameForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar,
    MultiTensorField2D<T, 2> &j, MultiTensorField2D<T, 2> &force,
    Array<T, 2> const &constantForce, Array<T, 2> const &angularVelocity,
    Array<T, 2> const &origin, bool incompressibleModel, Box2D domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  args.push_back(&force);
  applyProcessingFunctional(
      new ComputeRotatingFrameForce2D<T, Descriptor>(
          constantForce, angularVelocity, origin, incompressibleModel),
      domain, args);
}

/* ****************** FreeSurfaceComputeRotatingFrameForce2D
 * *************************************************** */

template <typename T, template <typename U> class Descriptor>
FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>::
    FreeSurfaceComputeRotatingFrameForce2D(Array<T, 2> const &constantForce_,
                                           Array<T, 2> const &angularVelocity_,
                                           Array<T, 2> const &origin_,
                                           bool incompressibleModel_)
    : constantForce(constantForce_),
      angularVelocity(angularVelocity_),
      origin(origin_),
      incompressibleModel(incompressibleModel_) {}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeRotatingFrameForce2D<
    T, Descriptor>::processGenericBlocks(Box2D domain,
                                         std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 5);
  BlockLattice2D<T, Descriptor> *lattice =
      dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
  PLB_ASSERT(lattice);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(rhoBar);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[2]);
  PLB_ASSERT(j);
  ScalarField2D<int> *flag = dynamic_cast<ScalarField2D<int> *>(blocks[3]);
  PLB_ASSERT(flag);
  TensorField2D<T, 2> *force = dynamic_cast<TensorField2D<T, 2> *>(blocks[4]);
  PLB_ASSERT(force);

  Dot2D loc = lattice->getLocation();

  Dot2D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
  Dot2D ofsJ = computeRelativeDisplacement(*lattice, *j);
  Dot2D ofsFlag = computeRelativeDisplacement(*lattice, *flag);
  Dot2D ofsForce = computeRelativeDisplacement(*lattice, *force);

  if (!util::isZero(norm(constantForce))) {
    T angularVelocityNorm = norm(angularVelocity);
    if (!util::isZero(angularVelocityNorm)) {
      Array<T, 2> normedAxis = angularVelocity / angularVelocityNorm;
      plint t = lattice->getTimeCounter().getTime();
      T theta = -angularVelocityNorm * t;
      constantForce = rotateAtOrigin(constantForce, normedAxis, theta);
    }
  }

  for (plint iX = domain.x0; iX <= domain.x1; iX++) {
    T x = iX + loc.x - origin[0];
    for (plint iY = domain.y0; iY <= domain.y1; iY++) {
      T y = iY + loc.y - origin[1];

      Array<T, 2> newForce(Array<T, 2>::zero());

      if (lbfoam::freeSurfaceFlag2D::isWet(
              flag->get(iX + ofsFlag.x, iY + ofsFlag.y))) {
        // Constant force.
        newForce = constantForce;

        // Coriolis and centripetal forces.
        Array<T, 2> r(x, y);
        Array<T, 2> velocity(j->get(iX + ofsJ.x, iY + ofsJ.y));
        if (!incompressibleModel) {
          T rho = Descriptor<T>::fullRho(
              rhoBar->get(iX + ofsRhoBar.x, iY + ofsRhoBar.y));
          velocity /= rho;
        }
        Array<T, 2> const &oldForce =
            force->get(iX + ofsForce.x, iY + ofsForce.y);
        velocity += (T)0.5 * oldForce;
        newForce +=
            -((T)2 * crossProduct(angularVelocity, velocity) +
              crossProduct(angularVelocity, crossProduct(angularVelocity, r)));
      }

      force->get(iX + ofsForce.x, iY + ofsForce.y) = newForce;
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>
    *FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>::clone() const {
  return new FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT
FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::nothing;          // Lattice
  modified[1] = modif::nothing;          // rhoBar
  modified[2] = modif::nothing;          // j
  modified[3] = modif::nothing;          // Flag
  modified[4] = modif::staticVariables;  // Force
}

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeRotatingFrameForce(
    MultiBlockLattice2D<T, Descriptor> &lattice, MultiScalarField2D<T> &rhoBar,
    MultiTensorField2D<T, 2> &j, MultiScalarField2D<int> &flag,
    MultiTensorField2D<T, 2> &force, Array<T, 2> const &constantForce,
    Array<T, 2> const &angularVelocity, Array<T, 2> const &origin,
    bool incompressibleModel, Box2D domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&lattice);
  args.push_back(&rhoBar);
  args.push_back(&j);
  args.push_back(&flag);
  args.push_back(&force);
  applyProcessingFunctional(
      new FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>(
          constantForce, angularVelocity, origin, incompressibleModel),
      domain, args);
}

}  // namespace lbfoam
}  // namespace plb

#endif  // BODY_FORCE_2D_HH
