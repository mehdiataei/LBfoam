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

#ifndef IMMERSED_WALLS_2D_HH
#define IMMERSED_WALLS_2D_HH

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "immersedWalls2D.h"

namespace plb {

namespace lbfoam {

/* ******** ReduceAxialTorqueImmersed2D ************************************ */

template <typename T>
ReduceAxialTorqueImmersed2D<T>::ReduceAxialTorqueImmersed2D(
    Array<T, 2> const& center_, Array<T, 2> const& unitaryAxis_,
    int reductionFlag_)
    : center(center_),
      unitaryAxis(unitaryAxis_),
      sum_torque_ids(Array<plint, 2>(this->getStatistics().subscribeSum(),
                                     this->getStatistics().subscribeSum())),
      reductionFlag(reductionFlag_) {}

template <typename T>
void ReduceAxialTorqueImmersed2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<Array<T, 2> > const& g = wallData->g;
  std::vector<int> const& flags = wallData->flags;
  Array<T, 2> offset = wallData->offset;
  PLB_ASSERT(vertices.size() == g.size());
  PLB_ASSERT(vertices.size() == flags.size());

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i];
    if (flags[i] == reductionFlag && closedOpenContained(vertex, domain)) {
      Array<T, 2> physVertex = vertex + offset;
      Array<T, 2> r(physVertex - center);
      r -= dot(r, unitaryAxis) * unitaryAxis;
      Array<T, 2> torque(crossProduct(r, g[i]));
      this->getStatistics().gatherSum(sum_torque_ids[0], torque[0]);
      this->getStatistics().gatherSum(sum_torque_ids[1], torque[1]);
    }
  }
}

template <typename T>
ReduceAxialTorqueImmersed2D<T>* ReduceAxialTorqueImmersed2D<T>::clone() const {
  return new ReduceAxialTorqueImmersed2D<T>(*this);
}

template <typename T>
void ReduceAxialTorqueImmersed2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceAxialTorqueImmersed2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T>
Array<T, 2> ReduceAxialTorqueImmersed2D<T>::getSumTorque() const {
  return Array<T, 2>(this->getStatistics().getSum(sum_torque_ids[0]),
                     this->getStatistics().getSum(sum_torque_ids[1]));
}

/* ******** ReduceImmersedForce2D ************************************ */

template <typename T>
ReduceImmersedForce2D<T>::ReduceImmersedForce2D(int reductionFlag_)
    : sum_g_ids(Array<plint, 2>(this->getStatistics().subscribeSum(),
                                this->getStatistics().subscribeSum())),
      reductionFlag(reductionFlag_) {}

template <typename T>
void ReduceImmersedForce2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<Array<T, 2> > const& g = wallData->g;
  std::vector<int> const& flags = wallData->flags;
  PLB_ASSERT(vertices.size() == g.size());
  PLB_ASSERT(vertices.size() == flags.size());

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i];
    if (flags[i] == reductionFlag && closedOpenContained(vertex, domain)) {
      this->getStatistics().gatherSum(sum_g_ids[0], g[i][0]);
      this->getStatistics().gatherSum(sum_g_ids[1], g[i][1]);
      ;
    }
  }
}

template <typename T>
ReduceImmersedForce2D<T>* ReduceImmersedForce2D<T>::clone() const {
  return new ReduceImmersedForce2D<T>(*this);
}

template <typename T>
void ReduceImmersedForce2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceImmersedForce2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T>
Array<T, 2> ReduceImmersedForce2D<T>::getSumG() const {
  return Array<T, 2>(this->getStatistics().getSum(sum_g_ids[0]),
                     this->getStatistics().getSum(sum_g_ids[1]));
}

/* ******** ReduceImmersedArea2D ************************************ */

template <typename T>
ReduceImmersedArea2D<T>::ReduceImmersedArea2D(int reductionFlag_)
    : sum_area_id(this->getStatistics().subscribeSum()),
      reductionFlag(reductionFlag_) {}

template <typename T>
void ReduceImmersedArea2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  std::vector<int> const& flags = wallData->flags;
  PLB_ASSERT(vertices.size() == areas.size());
  PLB_ASSERT(vertices.size() == flags.size());

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i];
    if (flags[i] == reductionFlag && closedOpenContained(vertex, domain)) {
      this->getStatistics().gatherSum(sum_area_id, areas[i]);
    }
  }
}

template <typename T>
ReduceImmersedArea2D<T>* ReduceImmersedArea2D<T>::clone() const {
  return new ReduceImmersedArea2D<T>(*this);
}

template <typename T>
void ReduceImmersedArea2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceImmersedArea2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

template <typename T>
T ReduceImmersedArea2D<T>::getSumArea() const {
  return this->getStatistics().getSum(sum_area_id);
}

/* ******** InamuroIteration2D ************************************ */

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction>::InamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_)
    : velFunction(velFunction_),
      tau(tau_),
      incompressibleModel(incompressibleModel_) {}

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 3);
  ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  TensorField2D<T, 2>* j = dynamic_cast<TensorField2D<T, 2>*>(blocks[1]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(j);
  PLB_ASSERT(container);

  Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  Array<T, 2> absOffset = wallData->offset;

  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  PLB_ASSERT(vertices.size() == areas.size());
  std::vector<Array<T, 2> > deltaG(vertices.size());
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());

  // In this iteration, the force is computed for every vertex.
  if (incompressibleModel) {
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      // Use the weighting function to compute the average momentum
      // and the average density on the surface vertex.
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunction2D<T>().W(r);
          averageJ += W * nextJ;
        }
      }
      // averageJ += (T)0.5*g[i];
      Array<T, 2> wallVelocity = velFunction(vertex + absOffset);
      deltaG[i] = areas[i] * (wallVelocity - averageJ);
      g[i] += deltaG[i];
    }
  } else {  // Compressible model.
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      T averageRhoBar = T();
      // Use the weighting function to compute the average momentum
      // and the average density on the surface vertex.
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          T nextRhoBar = rhoBar->get(pos[0], pos[1]);
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunction2D<T>().W(r);
          averageJ += W * nextJ;
          averageRhoBar += W * nextRhoBar;
        }
      }
      // averageJ += (T)0.5*g[i];
      Array<T, 2> wallVelocity = velFunction(vertex + absOffset);
      deltaG[i] =
          areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
      // g[i] += deltaG[i];
      g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
    }
  }

  // In this iteration, the force is applied from every vertex to the grid
  // nodes.
  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        nextJ += tau * W * deltaG[i];
        j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
      }
    }
  }
}

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction>* InamuroIteration2D<T, VelFunction>::clone()
    const {
  return new InamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::staticVariables;  // J
  modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT InamuroIteration2D<T, VelFunction>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** IndexedInamuroIteration2D ************************************ */

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction>::IndexedInamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_)
    : velFunction(velFunction_),
      tau(tau_),
      incompressibleModel(incompressibleModel_) {}

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 3);
  ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  TensorField2D<T, 2>* j = dynamic_cast<TensorField2D<T, 2>*>(blocks[1]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(j);
  PLB_ASSERT(container);

  Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  PLB_ASSERT(vertices.size() == areas.size());
  std::vector<Array<T, 2> > deltaG(vertices.size());
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());
  std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
  PLB_ASSERT(vertices.size() == globalVertexIds.size());

  if (incompressibleModel) {
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunctionInline<T>().W(r);
          averageJ += W * nextJ;
        }
      }
      // averageJ += (T)0.5*g[i];
      Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
      deltaG[i] = areas[i] * (wallVelocity - averageJ);
      g[i] += deltaG[i];
    }
  } else {  // Compressible model.
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      T averageRhoBar = T();
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          T nextRhoBar = rhoBar->get(pos[0], pos[1]);
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunctionInline<T>().W(r);
          averageJ += W * nextJ;
          averageRhoBar += W * nextRhoBar;
        }
      }
      // averageJ += (T)0.5*g[i];
      Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
      deltaG[i] =
          areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
      // g[i] += deltaG[i];
      g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
    }
  }

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunctionInline<T>().W(r);
        nextJ += tau * W * deltaG[i];
        j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
      }
    }
  }
}

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction>*
IndexedInamuroIteration2D<T, VelFunction>::clone() const {
  return new IndexedInamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::staticVariables;  // J
  modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT IndexedInamuroIteration2D<T, VelFunction>::appliesTo()
    const {
  return BlockDomain::bulk;
}

/* ******** ConstVelInamuroIteration2D ************************************ */

template <typename T>
ConstVelInamuroIteration2D<T>::ConstVelInamuroIteration2D(
    Array<T, 2> const& wallVelocity_, T tau_, bool incompressibleModel_)
    : wallVelocity(wallVelocity_),
      tau(tau_),
      incompressibleModel(incompressibleModel_) {}

template <typename T>
void ConstVelInamuroIteration2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 3);
  ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  TensorField2D<T, 2>* j = dynamic_cast<TensorField2D<T, 2>*>(blocks[1]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(j);
  PLB_ASSERT(container);

  Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  PLB_ASSERT(vertices.size() == areas.size());
  std::vector<Array<T, 2> > deltaG(vertices.size());
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());

  if (incompressibleModel) {
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunction2D<T>().W(r);
          averageJ += W * nextJ;
        }
      }
      // averageJ += (T)0.5*g[i];
      deltaG[i] = areas[i] * (wallVelocity - averageJ);
      g[i] += deltaG[i];
    }
  } else {  // Compressible model.
    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const& vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      Array<T, 2> averageJ;
      averageJ.resetToZero();
      T averageRhoBar = T();
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          T nextRhoBar = rhoBar->get(pos[0], pos[1]);
          Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
          Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
          T W = inamuroDeltaFunction2D<T>().W(r);
          averageJ += W * nextJ;
          averageRhoBar += W * nextRhoBar;
        }
      }
      // averageJ += (T)0.5*g[i];
      deltaG[i] =
          areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
      // g[i] += deltaG[i];
      g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
    }
  }

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        nextJ += tau * W * deltaG[i];
        j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
      }
    }
  }
}

template <typename T>
ConstVelInamuroIteration2D<T>* ConstVelInamuroIteration2D<T>::clone() const {
  return new ConstVelInamuroIteration2D<T>(*this);
}

template <typename T>
void ConstVelInamuroIteration2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::staticVariables;  // J
  modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ConstVelInamuroIteration2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** ComputeImmersedBoundaryForce2D ************************************
 */

template <typename T>
void ComputeImmersedBoundaryForce2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 2);
  TensorField2D<T, 2>* force = dynamic_cast<TensorField2D<T, 2>*>(blocks[0]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[1]);
  PLB_ASSERT(force);
  PLB_ASSERT(container);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());

  for (plint iX = 0; iX < force->getNx(); iX++) {
    for (plint iY = 0; iY < force->getNy(); iY++) {
      force->get(iX, iY).resetToZero();
    }
  }

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        force->get(pos[0], pos[1]) += W * g[i];
      }
    }
  }
}

template <typename T>
ComputeImmersedBoundaryForce2D<T>* ComputeImmersedBoundaryForce2D<T>::clone()
    const {
  return new ComputeImmersedBoundaryForce2D<T>(*this);
}

template <typename T>
void ComputeImmersedBoundaryForce2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Force
  modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ComputeImmersedBoundaryForce2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallData2D ************************************
 */

template <typename T>
InstantiateImmersedWallData2D<T>::InstantiateImmersedWallData2D(
    std::vector<Array<T, 2> > const& vertices_, std::vector<T> const& areas_,
    std::vector<Array<T, 2> > const& normals_)
    : vertices(vertices_), areas(areas_), normals(normals_) {
  PLB_ASSERT(vertices.size() == areas.size());
  PLB_ASSERT(normals.size() == 0 || normals.size() == areas.size());
}

template <typename T>
void InstantiateImmersedWallData2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);
  bool useNormals = normals.size() > 0;
  Dot2D location = container->getLocation();
  Array<T, 2> offset(location.x, location.y);
  ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
  Box2D extendedEnvelope(domain.enlarge(2));

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i] - offset;
    // Vertices which are close to the boundaries of the extendedEnvelope
    // are irrelevant, because they will act upon the bulk of the computational
    // domain through an Inamuro kernel, which at this distance is close to
    // zero. It is therefore OK, numerically speaking to exclude an
    // epsilon-margin close to these boundaries. Plus, it is required for
    // technical reasons, because if later on we pass across the boundaries of
    // the extendedEnvelope because of roundoff errors, the code will crash.
    static const T epsilon = 1.e-4;
    if (contained(vertex, extendedEnvelope, epsilon)) {
      wallData->vertices.push_back(vertex);
      wallData->areas.push_back(areas[i]);
      if (useNormals) {
        wallData->normals.push_back(normals[i]);
      }
      wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
      wallData->globalVertexIds.push_back(i);
    }
  }
  wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
  wallData->offset = offset;
  container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallData2D<T>* InstantiateImmersedWallData2D<T>::clone()
    const {
  return new InstantiateImmersedWallData2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallData2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallData2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallDataWithTagging2D
 * ************************************ */

template <typename T>
InstantiateImmersedWallDataWithTagging2D<T>::
    InstantiateImmersedWallDataWithTagging2D(
        std::vector<Array<T, 2> > const& vertices_,
        std::vector<T> const& areas_, int fluidFlag_)
    : vertices(vertices_), areas(areas_), fluidFlag(fluidFlag_) {
  PLB_ASSERT(vertices.size() == areas.size());
}

template <typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 2);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);
  Dot2D location = container->getLocation();
  Array<T, 2> offset(location.x, location.y);

  ScalarField2D<int>* flagMatrix = dynamic_cast<ScalarField2D<int>*>(blocks[1]);
  PLB_ASSERT(flagMatrix);
  Dot2D ofsFlag = computeRelativeDisplacement(*container, *flagMatrix);
  Array<plint, 2> flagDispl(ofsFlag.x, ofsFlag.y);

  ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
  Box2D extendedEnvelope(domain.enlarge(2));

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i] - offset;
    // Vertices which are close to the boundaries of the extendedEnvelope
    // are irrelevant, because they will act upon the bulk of the computational
    // domain through an Inamuro kernel, which at this distance is close to
    // zero. It is therefore OK, numerically speaking to exclude an
    // epsilon-margin close to these boundaries. Plus, it is required for
    // technical reasons, because if later on we pass across the boundaries of
    // the extendedEnvelope because of roundoff errors, the code will crash.
    static const T epsilon = 1.e-4;
    if (contained(vertex, extendedEnvelope, epsilon)) {
      wallData->vertices.push_back(vertex);
      wallData->areas.push_back(areas[i]);
      wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
      wallData->globalVertexIds.push_back(i);
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      bool hasFluidNeighbor = false;
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy) + flagDispl);
          if (flagMatrix->get(pos[0], pos[1]) == fluidFlag) {
            hasFluidNeighbor = true;
          }
        }
      }
      if (hasFluidNeighbor) {
        wallData->flags.push_back(0);
      } else {
        wallData->flags.push_back(1);
      }
    }
  }
  wallData->offset = offset;
  container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallDataWithTagging2D<T>*
InstantiateImmersedWallDataWithTagging2D<T>::clone() const {
  return new InstantiateImmersedWallDataWithTagging2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Container Block with triangle data.
  modified[1] = modif::nothing;          // Flag matrix.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithTagging2D<T>::appliesTo()
    const {
  return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallDataWithIndexedTagging2D
 * ************************************ */

template <typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>::
    InstantiateImmersedWallDataWithIndexedTagging2D(
        std::vector<Array<T, 2> > const& vertices_,
        std::vector<T> const& areas_, std::vector<int> const& flags_)
    : vertices(vertices_), areas(areas_), flags(flags_) {
  PLB_ASSERT(vertices.size() == areas.size());
  PLB_ASSERT(vertices.size() == flags.size());
}

template <typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);
  Dot2D location = container->getLocation();
  Array<T, 2> offset(location.x, location.y);

  ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
  Box2D extendedEnvelope(domain.enlarge(2));

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i] - offset;
    // Vertices which are close to the boundaries of the extendedEnvelope
    // are irrelevant, because they will act upon the bulk of the computational
    // domain through an Inamuro kernel, which at this distance is close to
    // zero. It is therefore OK, numerically speaking to exclude an
    // epsilon-margin close to these boundaries. Plus, it is required for
    // technical reasons, because if later on we pass across the boundaries of
    // the extendedEnvelope because of roundoff errors, the code will crash.
    static const T epsilon = 1.e-4;
    if (contained(vertex, extendedEnvelope, epsilon)) {
      wallData->vertices.push_back(vertex);
      wallData->areas.push_back(areas[i]);
      wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
      wallData->flags.push_back(flags[i]);
      wallData->globalVertexIds.push_back(i);
    }
  }
  wallData->offset = offset;
  container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>*
InstantiateImmersedWallDataWithIndexedTagging2D<T>::clone() const {
  return new InstantiateImmersedWallDataWithIndexedTagging2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT
InstantiateImmersedWallDataWithIndexedTagging2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** InstantiateSurfaceBlockData2D ************************************
 */

template <typename T>
InstantiateSurfaceBlockData2D<T>::InstantiateSurfaceBlockData2D(
    plint envelopeWidth_, std::vector<Array<T, 2> > const& vertices_)
    : envelopeWidth(envelopeWidth_), vertices(vertices_) {
  PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T>
void InstantiateSurfaceBlockData2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);

  Dot2D location = container->getLocation();
  Array<T, 2> offset(location.x, location.y);

  SurfaceBlockData2D<T>* surfaceData = new SurfaceBlockData2D<T>;
  Box2D extendedEnvelope(domain.enlarge(envelopeWidth));

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> vertex = vertices[i] - offset;
    if (containedInclusive(vertex, extendedEnvelope)) {
      surfaceData->vertices.push_back(vertex);
    }
  }
  surfaceData->offset = offset;
  container->setData(surfaceData);
}

template <typename T>
InstantiateSurfaceBlockData2D<T>* InstantiateSurfaceBlockData2D<T>::clone()
    const {
  return new InstantiateSurfaceBlockData2D<T>(*this);
}

template <typename T>
void InstantiateSurfaceBlockData2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateSurfaceBlockData2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** SurfaceOnLattice2D *************************************** */

template <typename T, typename U>
SurfaceOnLattice2D<T, U>::SurfaceOnLattice2D(U value_, plint envelopeWidth_)
    : value(value_), envelopeWidth(envelopeWidth_) {
  PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T, typename U>
void SurfaceOnLattice2D<T, U>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 2);
  ScalarField2D<U>* surfaceOnLattice =
      dynamic_cast<ScalarField2D<U>*>(blocks[0]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[1]);
  PLB_ASSERT(surfaceOnLattice);
  PLB_ASSERT(container);

  SurfaceBlockData2D<T>* surfaceData =
      dynamic_cast<SurfaceBlockData2D<T>*>(container->getData());
  PLB_ASSERT(surfaceData);

  std::vector<Array<T, 2> > const& vertices = surfaceData->vertices;

  plint w = envelopeWidth;
  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    // x  x  x . x  x  x
    for (plint dx = -(w - 1); dx <= +w; ++dx) {
      for (plint dy = -(w - 1); dy <= +w; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        if (contained(pos[0], pos[1], domain)) {
          surfaceOnLattice->get(pos[0], pos[1]) = value;
        }
      }
    }
  }
}

template <typename T, typename U>
SurfaceOnLattice2D<T, U>* SurfaceOnLattice2D<T, U>::clone() const {
  return new SurfaceOnLattice2D<T, U>(*this);
}

template <typename T, typename U>
void SurfaceOnLattice2D<T, U>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Surface values (ScalarField2D<U>).
  modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice2D<T, U>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** SurfaceOnLattice2D_N *************************************** */

template <typename T, typename U>
SurfaceOnLattice2D_N<T, U>::SurfaceOnLattice2D_N(U value_, plint envelopeWidth_)
    : value(value_), envelopeWidth(envelopeWidth_) {
  PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T, typename U>
void SurfaceOnLattice2D_N<T, U>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 2);
  NTensorField2D<U>* surfaceOnLattice =
      dynamic_cast<NTensorField2D<U>*>(blocks[0]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[1]);
  PLB_ASSERT(surfaceOnLattice);
  PLB_ASSERT(container);

  SurfaceBlockData2D<T>* surfaceData =
      dynamic_cast<SurfaceBlockData2D<T>*>(container->getData());
  PLB_ASSERT(surfaceData);

  std::vector<Array<T, 2> > const& vertices = surfaceData->vertices;

  plint w = envelopeWidth;
  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    // x  x  x . x  x  x
    for (plint dx = -(w - 1); dx <= +w; ++dx) {
      for (plint dy = -(w - 1); dy <= +w; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        if (contained(pos[0], pos[1], domain)) {
          *surfaceOnLattice->get(pos[0], pos[1]) = value;
        }
      }
    }
  }
}

template <typename T, typename U>
SurfaceOnLattice2D_N<T, U>* SurfaceOnLattice2D_N<T, U>::clone() const {
  return new SurfaceOnLattice2D_N<T, U>(*this);
}

template <typename T, typename U>
void SurfaceOnLattice2D_N<T, U>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::staticVariables;  // Surface values (NTensorField2D<U>).
  modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice2D_N<T, U>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** ResetForceStatistics2D ************************************ */

template <typename T>
void ResetForceStatistics2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 1);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
  PLB_ASSERT(container);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2> >& g = wallData->g;

  for (pluint i = 0; i < g.size(); i++) {
    g[i].resetToZero();
  }
}

template <typename T>
ResetForceStatistics2D<T>* ResetForceStatistics2D<T>::clone() const {
  return new ResetForceStatistics2D<T>(*this);
}

template <typename T>
void ResetForceStatistics2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ResetForceStatistics2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** RecomputeImmersedForce2D ************************************ */

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
RecomputeImmersedForce2D<T, Descriptor, NormalFunction>::
    RecomputeImmersedForce2D(NormalFunction normalFunction_, T omega_,
                             T densityOffset_, bool incompressibleModel_)
    : normalFunction(normalFunction_),
      omega(omega_),
      rho0(densityOffset_),
      incompressibleModel(incompressibleModel_) {
  PLB_ASSERT(densityOffset_ > (T)0);
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
void RecomputeImmersedForce2D<T, Descriptor, NormalFunction>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 3);

  ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  TensorField2D<T, SymmetricTensorImpl<T, 2>::n>* PiNeq =
      dynamic_cast<TensorField2D<T, SymmetricTensorImpl<T, 2>::n>*>(blocks[1]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(PiNeq);
  PLB_ASSERT(container);

  Dot2D ofsPN = computeRelativeDisplacement(*rhoBar, *PiNeq);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  PLB_ASSERT(vertices.size() == areas.size());
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());
  std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
  PLB_ASSERT(vertices.size() == globalVertexIds.size());

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> normal = normalFunction(globalVertexIds[i]);

    // Interpolate rhoBar and PiNeq on the vertex position.

    Array<T, 2> const& vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);

    T averageRhoBar = 0.0;
    Array<T, SymmetricTensorImpl<T, 2>::n> averagePiNeq;
    averagePiNeq.resetToZero();

    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        T nextRhoBar = rhoBar->get(pos[0], pos[1]);
        Array<T, SymmetricTensorImpl<T, 2>::n>& nextPiNeq =
            PiNeq->get(pos[0] + ofsPN.x, pos[1] + ofsPN.y);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        averageRhoBar += W * nextRhoBar;
        averagePiNeq += W * nextPiNeq;
      }
    }

    // Compute the force on the fluid at the vertex position.

    T averageRho = Descriptor<T>::fullRho(averageRhoBar);

    Array<T, 2> averagePi_n;
    SymmetricTensorImpl<T, 2>::matVectMult(averagePiNeq, normal, averagePi_n);

    // Here we want the force acting on the fluid from the solid. "normal"
    // points towards the fluid, this is why the minus sign in front of the area
    // is needed.
    if (incompressibleModel) {
      g[i] = -areas[i] *
             (-(averageRho - rho0) * Descriptor<T>::cs2 * normal +
              (omega / (T)2. - (T)1.) * averagePi_n);  // Incompressible vision
    } else {
      g[i] = -areas[i] *
             (-(averageRho - rho0) * Descriptor<T>::cs2 * normal +
              Descriptor<T>::invRho(averageRhoBar) * (omega / (T)2. - (T)1.) *
                  averagePi_n);  // Compressible vision
    }
  }
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
RecomputeImmersedForce2D<T, Descriptor, NormalFunction>*
RecomputeImmersedForce2D<T, Descriptor, NormalFunction>::clone() const {
  return new RecomputeImmersedForce2D<T, Descriptor, NormalFunction>(*this);
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
void RecomputeImmersedForce2D<T, Descriptor, NormalFunction>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // RhoBar
  modified[1] = modif::nothing;  // PiNeq
  modified[2] = modif::nothing;  // Container with triangle data
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
BlockDomain::DomainT
RecomputeImmersedForce2D<T, Descriptor, NormalFunction>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** OpenSurfaceImmersedForce2D ************************************ */

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>::
    OpenSurfaceImmersedForce2D(NormalFunction normalFunction_, T omega_,
                               T densityOffset_, bool incompressibleModel_)
    : normalFunction(normalFunction_),
      omega(omega_),
      rho0(densityOffset_),
      incompressibleModel(incompressibleModel_) {
  PLB_ASSERT(densityOffset_ > (T)0);
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
void OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  PLB_PRECONDITION(blocks.size() == 3);

  ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  TensorField2D<T, SymmetricTensorImpl<T, 2>::n>* PiNeq =
      dynamic_cast<TensorField2D<T, SymmetricTensorImpl<T, 2>::n>*>(blocks[1]);
  AtomicContainerBlock2D* container =
      dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(PiNeq);
  PLB_ASSERT(container);

  Dot2D ofsPN = computeRelativeDisplacement(*rhoBar, *PiNeq);

  ImmersedWallData2D<T>* wallData =
      dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2> > const& vertices = wallData->vertices;
  std::vector<T> const& areas = wallData->areas;
  PLB_ASSERT(vertices.size() == areas.size());
  std::vector<Array<T, 2> >& g = wallData->g;
  PLB_ASSERT(vertices.size() == g.size());
  std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
  PLB_ASSERT(vertices.size() == globalVertexIds.size());

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> normal = normalFunction(globalVertexIds[i]);

    Array<T, 2> const& vertex = vertices[i];

    // Interpolate rhoBar and PiNeq on the vertex position.
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);

    T averageRhoBar1 = 0.0, averageRhoBar2 = 0.0;
    Array<T, SymmetricTensorImpl<T, 2>::n> averagePiNeq1, averagePiNeq2;
    averagePiNeq1.resetToZero();
    averagePiNeq2.resetToZero();
    plint n1 = 0, n2 = 0;
    T w1 = T(), w2 = T();

    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        T nextRhoBar = rhoBar->get(pos[0], pos[1]);
        Array<T, SymmetricTensorImpl<T, 2>::n>& nextPiNeq =
            PiNeq->get(pos[0] + ofsPN.x, pos[1] + ofsPN.y);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        if (dot(r, normal) > 0) {
          averageRhoBar1 += W * nextRhoBar;
          averagePiNeq1 += W * nextPiNeq;
          w1 += W;
          ++n1;
        } else {
          averageRhoBar2 += W * nextRhoBar;
          averagePiNeq2 += W * nextPiNeq;
          w2 += W;
          ++n2;
        }
      }
    }

    T averageRho1 = rho0;
    T averageRho2 = rho0;
    if (n1 > 0) {
      averageRhoBar1 /= w1;
      averageRho1 = Descriptor<T>::fullRho(averageRhoBar1);
      averagePiNeq1 /= w1;
    }
    if (n2 > 0) {
      averageRhoBar2 /= w2;
      averageRho2 = Descriptor<T>::fullRho(averageRhoBar2);
      averagePiNeq2 /= w2;
    }

    // Compute the force on the fluid at the vertex position.
    Array<T, 2> averagePi_n1;
    SymmetricTensorImpl<T, 2>::matVectMult(averagePiNeq1, normal, averagePi_n1);
    // Important: on the reverse side, take the negative normal.
    Array<T, 2> averagePi_n2;
    SymmetricTensorImpl<T, 2>::matVectMult(averagePiNeq2, -normal,
                                           averagePi_n2);

    // Here we want the force acting on the fluid from the solid. "normal"
    // points towards the fluid, this is why the minus sign in front of the area
    // is needed.
    if (incompressibleModel) {
      g[i] = -areas[i] *
             ((averageRho2 - averageRho1) * Descriptor<T>::cs2 * normal +
              (omega / (T)2. - (T)1.) *
                  (averagePi_n1 + averagePi_n2));  // Incompressible vision
    } else {
      g[i] = -areas[i] *
             ((averageRho2 - averageRho1) * Descriptor<T>::cs2 * normal +
              Descriptor<T>::invRho(averageRhoBar1) * (omega / (T)2. - (T)1.) *
                  averagePi_n1 +  // Compressible vision
              Descriptor<T>::invRho(averageRhoBar2) * (omega / (T)2. - (T)1.) *
                  averagePi_n2);
    }
  }
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>*
OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>::clone() const {
  return new OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>(*this);
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
void OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;  // RhoBar
  modified[1] = modif::nothing;  // PiNeq
  modified[2] = modif::nothing;  // Container with triangle data
}

template <typename T, template <typename U> class Descriptor,
          class NormalFunction>
BlockDomain::DomainT
OpenSurfaceImmersedForce2D<T, Descriptor, NormalFunction>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** TwoPhaseInamuroParam2D ************************************ */

template <typename T>
T TwoPhaseInamuroParam2D<T>::area(plint i) const {
  PLB_ASSERT((pluint)i < numVertices);
  return wallData->areas[i];
}

template <typename T>
Array<T, 2>& TwoPhaseInamuroParam2D<T>::g(plint i) {
  PLB_ASSERT((pluint)i < numVertices);
  return wallData->g[i];
}

template <typename T>
Array<T, 2> TwoPhaseInamuroParam2D<T>::vertex(plint i) const {
  PLB_ASSERT((pluint)i < numVertices);
  return wallData->vertices[i];
}

template <typename T>
Array<T, 2> TwoPhaseInamuroParam2D<T>::absoluteVertex(plint i) const {
  PLB_ASSERT((pluint)i < numVertices);
  return wallData->vertices[i] + absOffset;
}

template <typename T>
Array<plint, 2> TwoPhaseInamuroParam2D<T>::intVertex(plint i) const {
  PLB_ASSERT((pluint)i < numVertices);
  Array<T, 2> vertex = wallData->vertices[i];
  return Array<plint, 2>((plint)vertex[0], (plint)vertex[1]);
}

template <typename T>
T TwoPhaseInamuroParam2D<T>::rhoBar(plint iX, plint iY) const {
  int flag = getFlag(iX, iY);
  if (flag == freeSurfaceFlag2D::empty) {
    return rhoBar2_->get(iX + ofsRhoBar2.x, iY + ofsRhoBar2.y);
  } else {
    return rhoBar_->get(iX, iY);
  }
}

template <typename T>
Array<T, 2> TwoPhaseInamuroParam2D<T>::j(plint iX, plint iY) const {
  int flag = getFlag(iX, iY);
  if (flag == freeSurfaceFlag2D::empty) {
    return j2_->get(iX + ofsJ2.x, iY + ofsJ2.y);
  } else {
    return j_->get(iX + ofsJ.x, iY + ofsJ.y);
  }
}

template <typename T>
void TwoPhaseInamuroParam2D<T>::addToJ(plint iX, plint iY, Array<T, 2> deltaJ) {
  int flag = getFlag(iX, iY);
  if (flag == freeSurfaceFlag2D::interface) {
    j_->get(iX + ofsJ.x, iY + ofsJ.y) += deltaJ;
    j2_->get(iX + ofsJ2.x, iY + ofsJ2.y) += deltaJ;
  } else if (flag == freeSurfaceFlag2D::empty) {
    j2_->get(iX + ofsJ2.x, iY + ofsJ2.y) += deltaJ;
  } else {
    j_->get(iX + ofsJ.x, iY + ofsJ.y) += deltaJ;
  }
}

template <typename T>
int TwoPhaseInamuroParam2D<T>::getFlag(plint iX, plint iY) const {
  return flag_->get(iX + ofsFlag.x, iY + ofsFlag.y);
}

template <typename T>
pluint TwoPhaseInamuroParam2D<T>::getGlobalVertexId(plint i) const {
  PLB_ASSERT((pluint)i < numVertices);
  return wallData->globalVertexIds[i];
}

template <typename T>
T TwoPhaseInamuroParam2D<T>::getTau(plint iX, plint iY) const {
  T vf = volumeFraction_->get(iX + ofsVF.x, iY + ofsVF.y);
  return tau * vf + tau2 * (1. - vf);
}

template <typename T>
TwoPhaseInamuroParam2D<T>::TwoPhaseInamuroParam2D(
    std::vector<AtomicBlock2D*>& blocks, T tau_, T tau2_)
    : tau(tau_), tau2(tau2_) {
  PLB_PRECONDITION(blocks.size() == 7);
  rhoBar_ = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
  rhoBar2_ = dynamic_cast<ScalarField2D<T>*>(blocks[1]);
  j_ = dynamic_cast<TensorField2D<T, 2>*>(blocks[2]);
  j2_ = dynamic_cast<TensorField2D<T, 2>*>(blocks[3]);
  flag_ = dynamic_cast<ScalarField2D<int>*>(blocks[4]);
  volumeFraction_ = dynamic_cast<ScalarField2D<T>*>(blocks[5]);
  container = dynamic_cast<AtomicContainerBlock2D*>(blocks[6]);

  PLB_ASSERT(rhoBar_);
  PLB_ASSERT(rhoBar2_);
  PLB_ASSERT(j_);
  PLB_ASSERT(j2_);
  PLB_ASSERT(flag_);
  PLB_ASSERT(volumeFraction_);
  PLB_ASSERT(container);

  ofsRhoBar2 = computeRelativeDisplacement(*rhoBar_, *rhoBar2_);
  ofsJ = computeRelativeDisplacement(*rhoBar_, *j_);
  ofsJ2 = computeRelativeDisplacement(*rhoBar_, *j2_);
  ofsFlag = computeRelativeDisplacement(*rhoBar_, *flag_);
  ofsVF = computeRelativeDisplacement(*rhoBar_, *volumeFraction_);

  wallData = dynamic_cast<ImmersedWallData2D<T>*>(container->getData());
  PLB_ASSERT(wallData);
  absOffset = wallData->offset;

  numVertices = wallData->vertices.size();
  PLB_ASSERT(numVertices == wallData->areas.size());
  PLB_ASSERT(numVertices == wallData->g.size());
}

/* ******** TwoPhaseInamuroIteration2D ************************************ */

template <typename T, class VelFunction>
TwoPhaseInamuroIteration2D<T, VelFunction>::TwoPhaseInamuroIteration2D(
    VelFunction velFunction_, T tau_, T tau2_)
    : velFunction(velFunction_), tau(tau_), tau2(tau2_) {}

template <typename T, class VelFunction>
void TwoPhaseInamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
  std::vector<Array<T, 2> > deltaG(param.getNumVertices());

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    Array<T, 2> averageJ;
    averageJ.resetToZero();
    T averageRhoBar = T();
    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        T nextRhoBar = param.rhoBar(pos[0], pos[1]);
        Array<T, 2> nextJ = param.j(pos[0], pos[1]);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        averageJ += W * nextJ;
        averageRhoBar += W * nextRhoBar;
      }
    }
    // averageJ += (T)0.5*param.g(i);
    Array<T, 2> wallVelocity = velFunction(param.absoluteVertex(i));
    deltaG[i] =
        param.area(i) * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
    // param.g(i) += deltaG[i];
    param.g(i) += deltaG[i] / ((T)1.0 + averageRhoBar);
  }

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        param.addToJ(pos[0], pos[1],
                     param.getTau(pos[0], pos[1]) * W * deltaG[i]);
      }
    }
  }
}

template <typename T, class VelFunction>
TwoPhaseInamuroIteration2D<T, VelFunction>*
TwoPhaseInamuroIteration2D<T, VelFunction>::clone() const {
  return new TwoPhaseInamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void TwoPhaseInamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::nothing;          // RhoBar2
  modified[2] = modif::staticVariables;  // j
  modified[3] = modif::staticVariables;  // j2
  modified[4] = modif::nothing;          // flag
  modified[5] = modif::nothing;          // volume fraction
  modified[6] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT TwoPhaseInamuroIteration2D<T, VelFunction>::appliesTo()
    const {
  return BlockDomain::bulk;
}

/* ******** TwoPhaseIndexedInamuroIteration2D
 * ************************************ */

template <typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration2D<
    T, VelFunction>::TwoPhaseIndexedInamuroIteration2D(VelFunction velFunction_,
                                                       T tau_, T tau2_)
    : velFunction(velFunction_), tau(tau_), tau2(tau2_) {}

template <typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
  std::vector<Array<T, 2> > deltaG(param.getNumVertices());

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    Array<T, 2> averageJ;
    averageJ.resetToZero();
    T averageRhoBar = T();
    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        T nextRhoBar = param.rhoBar(pos[0], pos[1]);
        Array<T, 2> nextJ = param.j(pos[0], pos[1]);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        averageJ += W * nextJ;
        averageRhoBar += W * nextRhoBar;
      }
    }
    // averageJ += (T)0.5*param.g(i);
    Array<T, 2> wallVelocity = velFunction(param.getGlobalVertexId(i));
    deltaG[i] =
        param.area(i) * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
    // param.g(i) += deltaG[i];
    param.g(i) += deltaG[i] / ((T)1.0 + averageRhoBar);
  }

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        param.addToJ(pos[0], pos[1],
                     param.getTau(pos[0], pos[1]) * W * deltaG[i]);
      }
    }
  }
}

template <typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration2D<T, VelFunction>*
TwoPhaseIndexedInamuroIteration2D<T, VelFunction>::clone() const {
  return new TwoPhaseIndexedInamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::nothing;          // RhoBar2
  modified[2] = modif::staticVariables;  // j
  modified[3] = modif::staticVariables;  // j2
  modified[4] = modif::nothing;          // flag
  modified[5] = modif::nothing;          // volume fraction
  modified[6] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT
TwoPhaseIndexedInamuroIteration2D<T, VelFunction>::appliesTo() const {
  return BlockDomain::bulk;
}

/* ******** TwoPhaseConstVelInamuroIteration2D
 * ************************************ */

template <typename T>
TwoPhaseConstVelInamuroIteration2D<T>::TwoPhaseConstVelInamuroIteration2D(
    Array<T, 2> const& wallVelocity_, T tau_, T tau2_)
    : wallVelocity(wallVelocity_), tau(tau_), tau2(tau2_) {}

template <typename T>
void TwoPhaseConstVelInamuroIteration2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D*> blocks) {
  TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
  std::vector<Array<T, 2> > deltaG(param.getNumVertices());

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    Array<T, 2> averageJ;
    averageJ.resetToZero();
    T averageRhoBar = T();
    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        T nextRhoBar = param.rhoBar(pos[0], pos[1]);
        Array<T, 2> nextJ = param.j(pos[0], pos[1]);
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        averageJ += W * nextJ;
        averageRhoBar += W * nextRhoBar;
      }
    }
    // averageJ += (T)0.5*param.g(i);
    deltaG[i] =
        param.area(i) * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
    // param.g(i) += deltaG[i];
    param.g(i) += deltaG[i] / ((T)1.0 + averageRhoBar);
  }

  for (pluint i = 0; i < param.getNumVertices(); ++i) {
    Array<T, 2> vertex(param.vertex(i));
    Array<plint, 2> intPos(param.intVertex(i));

    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        Array<T, 2> r(pos[0] - vertex[0], pos[1] - vertex[1]);
        T W = inamuroDeltaFunction2D<T>().W(r);
        param.addToJ(pos[0], pos[1],
                     param.getTau(pos[0], pos[1]) * W * deltaG[i]);
      }
    }
  }
}

template <typename T>
TwoPhaseConstVelInamuroIteration2D<T>*
TwoPhaseConstVelInamuroIteration2D<T>::clone() const {
  return new TwoPhaseConstVelInamuroIteration2D<T>(*this);
}

template <typename T>
void TwoPhaseConstVelInamuroIteration2D<T>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;          // RhoBar
  modified[1] = modif::nothing;          // RhoBar2
  modified[2] = modif::staticVariables;  // j
  modified[3] = modif::staticVariables;  // j2
  modified[4] = modif::nothing;          // flag
  modified[5] = modif::nothing;          // Volume fraction
  modified[6] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT TwoPhaseConstVelInamuroIteration2D<T>::appliesTo() const {
  return BlockDomain::bulk;
}

}  // namespace lbfoam
}  // namespace plb

#endif  // IMMERSED_WALLS_2D_HH
