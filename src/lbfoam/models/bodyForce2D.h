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

#ifndef BODY_FORCE_2D_H
#define BODY_FORCE_2D_H

#include <vector>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

// Implementation of the momentum correction algorithm for applying a constant
// body force.
template <typename T, template <typename U> class Descriptor>
class AddConstForceToMomentum2D : public BoxProcessingFunctional2D {
 public:
  AddConstForceToMomentum2D(Array<T, 2> const& force_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual AddConstForceToMomentum2D<T, Descriptor>* clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;
  virtual BlockDomain::DomainT appliesTo() const;

 private:
  Array<T, 2> force;
};

template <typename T, template <typename U> class Descriptor>
void addConstForceToMomentum(MultiBlockLattice2D<T, Descriptor>& lattice,
                             MultiScalarField2D<T>& rhoBar,
                             MultiTensorField2D<T, 2>& j,
                             Array<T, 2> const& force, Box2D const& domain);

// Implementation of the momentum correction algorithm for applying a body force
// given by a user-provided function.
template <typename T, template <typename U> class Descriptor,
          class ForceFunction>
class AddCustomForceToMomentum2D : public BoxProcessingFunctional2D {
 public:
  AddCustomForceToMomentum2D(ForceFunction f_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual AddCustomForceToMomentum2D<T, Descriptor, ForceFunction>* clone()
      const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;
  virtual BlockDomain::DomainT appliesTo() const;

 private:
  ForceFunction f;
};

template <typename T, template <class U> class Descriptor, class ForceFunction>
void addCustomForceToMomentum(MultiBlockLattice2D<T, Descriptor>& lattice,
                              MultiScalarField2D<T>& rhoBar,
                              MultiTensorField2D<T, 2>& j, ForceFunction f,
                              Box2D const& domain);

// Implementation of the momentum correction algorithm for applying a body
// force.
template <typename T, template <typename U> class Descriptor>
class AddForceToMomentum2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual AddForceToMomentum2D<T, Descriptor>* clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;
  virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
void addForceToMomentum(MultiBlockLattice2D<T, Descriptor>& lattice,
                        MultiScalarField2D<T>& rhoBar,
                        MultiTensorField2D<T, 2>& j,
                        MultiTensorField2D<T, 2>& force, Box2D const& domain);

// Implementation of the momentum correction algorithm for applying a body force
// only to the wet nodes.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddForceToMomentum2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual FreeSurfaceAddForceToMomentum2D<T, Descriptor>* clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;
  virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
void freeSurfaceAddForceToMomentum(MultiBlockLattice2D<T, Descriptor>& lattice,
                                   MultiScalarField2D<T>& rhoBar,
                                   MultiTensorField2D<T, 2>& j,
                                   MultiScalarField2D<int>& flag,
                                   MultiTensorField2D<T, 2>& force,
                                   Box2D const& domain);

// Transform a constant force in a rotating reference frame, and add the
// Coriolis and centripetal forces. The Coriolis force depends on the local
// velocity. To compute this velocity we need to know the body force. For this,
// we use the body force at the previous iteration, which is expected to be
// contained in the "force" block passed to processGenericBlocks (fourth block
// in the vector). This data processor overwrites the "force" block with the new
// body force computed.
template <typename T, template <typename U> class Descriptor>
class ComputeRotatingFrameForce2D : public BoxProcessingFunctional2D {
 public:
  ComputeRotatingFrameForce2D(Array<T, 2> const& constantForce_,
                              Array<T, 2> const& angularVelocity_,
                              Array<T, 2> const& origin_,
                              bool incompressibleModel_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual ComputeRotatingFrameForce2D<T, Descriptor>* clone() const;
  virtual BlockDomain::DomainT appliesTo() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;

 private:
  Array<T, 2> constantForce;
  Array<T, 2> angularVelocity;
  Array<T, 2> origin;
  bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor>
void computeRotatingFrameForce(
    MultiBlockLattice2D<T, Descriptor>& lattice, MultiScalarField2D<T>& rhoBar,
    MultiTensorField2D<T, 2>& j, MultiTensorField2D<T, 2>& force,
    Array<T, 2> const& constantForce, Array<T, 2> const& angularVelocity,
    Array<T, 2> const& origin, bool incompressibleModel, Box2D domain);

// Transform a constant force in a rotating reference frame, and add the
// Coriolis and centripetal forces only on wet nodes (on the rest of the cells,
// the force is set to zero). The Coriolis force depends on the local velocity.
// To compute this velocity we need to know the body force. For this, we use the
// body force at the previous iteration, which is expected to be contained in
// the "force" block passed to processGenericBlocks (fifth block in the vector).
// This data processor overwrites the "force" block with the new body force
// computed.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeRotatingFrameForce2D
    : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceComputeRotatingFrameForce2D(Array<T, 2> const& constantForce_,
                                         Array<T, 2> const& angularVelocity_,
                                         Array<T, 2> const& origin_,
                                         bool incompressibleModel_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D*> blocks);
  virtual FreeSurfaceComputeRotatingFrameForce2D<T, Descriptor>* clone() const;
  virtual BlockDomain::DomainT appliesTo() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT>& modified) const;

 private:
  Array<T, 2> constantForce;
  Array<T, 2> angularVelocity;
  Array<T, 2> origin;
  bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor>
void freeSurfaceComputeRotatingFrameForce(
    MultiBlockLattice2D<T, Descriptor>& lattice, MultiScalarField2D<T>& rhoBar,
    MultiTensorField2D<T, 2>& j, MultiScalarField2D<int>& flag,
    MultiTensorField2D<T, 2>& force, Array<T, 2> const& constantForce,
    Array<T, 2> const& angularVelocity, Array<T, 2> const& origin,
    bool incompressibleModel, Box2D domain);

}  // namespace plb

#endif  // BODY_FORCE_2D_H
