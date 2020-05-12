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

#ifndef FREE_SURFACE_MODEL_2D_H
#define FREE_SURFACE_MODEL_2D_H

#include <algorithm>
#include <map>

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "basicDynamics/dynamicsProcessor2D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/models/freeSurfaceInitializer2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "lbfoam/models/immersedWalls2D.h"
#include "multiBlock/coupling2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/group2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeNormals2D : public BoxProcessingFunctional2D {
 public:
  virtual FreeSurfaceComputeNormals2D<T, Descriptor> *clone() const {
    return new FreeSurfaceComputeNormals2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::staticVariables;  // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceGeometry2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceGeometry2D(T contactAngle_) : contactAngle(contactAngle_) {
    // The contact angle must take values between 0 and 180 degrees. If it is
    // negative, this means that contact angle effects will not be modeled.
    PLB_ASSERT(util::lessEqual(contactAngle, (T)180));

    if (util::lessThan(contactAngle, (T)0)) {
      useContactAngle = 0;
    } else {
      useContactAngle = 1;
    }
  }
  virtual FreeSurfaceGeometry2D<T, Descriptor> *clone() const {
    return new FreeSurfaceGeometry2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::staticVariables;  // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::staticVariables;  // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  ScalarField2D<int> *getInterfaceFlags(
      Box2D domain, FreeSurfaceProcessorParam2D<T, Descriptor> &param);
  void computeHeights2D(FreeSurfaceProcessorParam2D<T, Descriptor> &param,
                        int integrationDirection, plint iX, plint iY, T h[3]);
  //    void computeHeights2DC(FreeSurfaceProcessorParam2D<T,Descriptor>& param,
  //    Array<int,2>& wallTangent0, //disabled for CA mehdi
  //             int integrationDirection, plint iX, plint iY, T h[3]);
 private:
  enum {
    unTagged = 0,
    notInterface = 1,
    regular = 2,
    contactLine = 4,
    adjacent = 10
  };

 private:
  T contactAngle;
  int useContactAngle;
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeCurvature2D : public BoxProcessingFunctional2D {
 public:
  typedef T (*ContactAngleFunction)(
      T x, T y);  // Returns the contact angle in degrees.
 public:
  FreeSurfaceComputeCurvature2D(T contactAngle_)
      : contactAngle(contactAngle_), contactAngleFunction(0) {
    // The contact angle must take values between 0 and 180 degrees. If it is
    // negative, this means that contact angle effects will not be modeled.
    PLB_ASSERT(util::lessEqual(contactAngle, (T)180));

    if (util::lessThan(contactAngle, (T)0)) {
      useContactAngle = 0;
    } else {
      useContactAngle = 1;
    }

    if (useContactAngle) {
      T pi = std::acos((T)-1);
      contactAngle *= pi / (T)180;
    }
  }
  FreeSurfaceComputeCurvature2D(ContactAngleFunction contactAngleFunction_)
      : contactAngle(-1.0), contactAngleFunction(contactAngleFunction_) {
    // The contact angle must take values between 0 and 180 degrees.
    // If the function pointer is 0, this means that contact angle effects will
    // not be modeled.
    if (contactAngleFunction == 0) {
      useContactAngle = 0;
    } else {
      useContactAngle = 1;
    }
  }
  virtual FreeSurfaceComputeCurvature2D<T, Descriptor> *clone() const {
    return new FreeSurfaceComputeCurvature2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::staticVariables;  // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T contactAngle;
  ContactAngleFunction contactAngleFunction;
  int useContactAngle;
};

/// Compute the mass balance on every node in the domain, and store in mass
/// matrix.
/** Input:
 *   - Flag-status:   needed in bulk+1
 *   - Mass:          needed in bulk
 *   - Volume fraction: needed in bulk
 *   - Populations:   needed in bulk+1
 * Output:
 *   - mass.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceMassChange2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceMassChange2D<T, Descriptor> *clone() const {
    return new FreeSurfaceMassChange2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;          // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }
};

/// Completion scheme on the post-collide populations on interface cells.
/** Input:
 *   - Flag-status:   needed in bulk+1
 *   - Volume fraction: needed in bulk+1
 *   - Populations:   needed in bulk+1
 *   - Momentum:      needed in bulk+1
 *   - Density:       needed in bulk+1
 * Output:
 *   - Populations.
 **/
// ASK: This data processor loops over the whole volume. Is this really
//      necessary, or could one of the lists be used instead?
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceCompletion2D : public BoxProcessingFunctional2D {
 public:
  virtual FreeSurfaceCompletion2D<T, Descriptor> *clone() const {
    return new FreeSurfaceCompletion2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid. Should be: staticVariables.
    modified[1] = modif::nothing;  // rhoBar.
    modified[2] = modif::nothing;  // j.
    modified[3] = modif::nothing;  // Mass.
    modified[4] = modif::nothing;  // Volume fraction.
    modified[5] = modif::nothing;  // Flag-status.
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface-lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }
};

/// Compute and store mass, volume-fraction and macroscopic variables.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - mass, volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceMacroscopic2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceMacroscopic2D(bool incompressibleModel_)
      : incompressibleModel(incompressibleModel_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceMacroscopic2D<T, Descriptor> *clone() const {
    return new FreeSurfaceMacroscopic2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] =
        modif::nothing;  // Fluid. Should be: staticVariables (maybe not...).
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass. Should be: staticVariables.
    modified[4] = modif::staticVariables;  // Volume fraction.
    modified[5] =
        modif::nothing;  // Flag-status. TODO: used to be staticVariables...
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface-lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }

 private:
  bool incompressibleModel;
};

/// Add the surface tension contribution.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddSurfaceTension2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceAddSurfaceTension2D(T surfaceTension_, bool incompressibleModel_)
      : surfaceTension(surfaceTension_),
        incompressibleModel(incompressibleModel_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceAddSurfaceTension2D<T, Descriptor> *clone() const {
    return new FreeSurfaceAddSurfaceTension2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] =
        modif::nothing;  // Fluid. Should be: staticVariables (maybe not...).
    modified[1] = modif::staticVariables;  // rhoBar.
    if (incompressibleModel) {
      modified[2] = modif::nothing;  // j.
    } else {
      modified[2] = modif::staticVariables;  // j.
    }
    modified[3] = modif::nothing;  // Mass. TODO: used to be staticVariables...
    modified[4] = modif::staticVariables;  // Volume fraction.
    modified[5] =
        modif::nothing;  // Flag-status. TODO: used to be staticVariables...
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface-lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }

 private:
  T surfaceTension;
  bool incompressibleModel;
};

/// Add the surface tension contribution. The surface tension coefficient is
/// read from a scalar field.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddSurfaceTensionFromScalarField2D
    : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceAddSurfaceTensionFromScalarField2D(bool incompressibleModel_)
      : incompressibleModel(incompressibleModel_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceAddSurfaceTensionFromScalarField2D<T, Descriptor> *clone()
      const {
    return new FreeSurfaceAddSurfaceTensionFromScalarField2D<T, Descriptor>(
        *this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] =
        modif::nothing;  // Fluid. Should be: staticVariables (maybe not...).
    modified[1] = modif::staticVariables;  // rhoBar.
    if (incompressibleModel) {
      modified[2] = modif::nothing;  // j.
    } else {
      modified[2] = modif::staticVariables;  // j.
    }
    modified[3] = modif::nothing;  // Mass. TODO: used to be staticVariables...
    modified[4] = modif::staticVariables;  // Volume fraction.
    modified[5] =
        modif::nothing;  // Flag-status. TODO: used to be staticVariables...
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface-lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.

    modified[10] = modif::nothing;  // Surface tension coefficient.
  }

 private:
  bool incompressibleModel;
};

/// Stabilization scheme on the post-collide populations.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Populations:   needed in bulk
 *   - Momentum:      needed in bulk
 *   - Density:       needed in bulk
 * Output:
 *   - Populations.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceStabilize2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceStabilize2D() {}
  virtual FreeSurfaceStabilize2D<T, Descriptor> *clone() const {
    return new FreeSurfaceStabilize2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid (should be staticVariables).
    modified[1] = modif::nothing;  // rhoBar.
    modified[2] = modif::nothing;  // j.
    modified[3] = modif::nothing;  // Mass.
    modified[4] = modif::nothing;  // Volume fraction.
    modified[5] = modif::nothing;  // Flag-status.
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }
};

/// Based on the current flag status, decide, upon the value of mass fraction,
/// which nodes shall
///   switch state.
/** Input:
 *   - Volume fraction: needed in bulk+2
 *   - Flag-status:   needed in bulk+2
 * Output:
 *   - interface-to-fluid list: defined in bulk+2
 *   - interface-to-empty list: defined in bulk+1
 *   - empty-to-interface list: defined in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeInterfaceLists2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceComputeInterfaceLists2D<T, Descriptor> *clone() const {
    return new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid (not used in this processor).
    modified[1] = modif::nothing;  // rhoBar.
    modified[2] = modif::nothing;  // j.
    modified[3] = modif::nothing;  // Mass (not used in this processor).
    modified[4] = modif::nothing;  // Volume fraction, read-only.
    modified[5] = modif::nothing;  // Flag-status, read-only.
    modified[6] = modif::nothing;  // Normal.
    modified[7] =
        modif::staticVariables;  // Interface-lists; all lists are created here.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }

 private:
  static T kappa;  // Safety threshold for state-change, to prevent
                   // back-and-forth oscillations.
};

/** Input:
 *   - interface-to-fluid list: needed in bulk+1
 *   - interface-to-empty list: needed in bulk+1
 *   - density: needed in bulk+1
 *   - mass:    needed in bulk+1
 *   - flag:    needed in bulk+1
 * Output:
 *   - flag, dynamics, mass, volumeFraction, density, force, momentum
 *   - mass-excess-list: defined in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceIniInterfaceToAnyNodes2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceIniInterfaceToAnyNodes2D(T rhoDefault_);
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);

  virtual FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor> *clone() const {
    return new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid. Gets assigned new dynamics. Should
                                   // be: dataStructure
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::nothing;          // j. Should be: staticVariables.
    modified[3] =
        modif::staticVariables;  // Mass. Is redistributed and initialized from
                                 // neighborying density.
    modified[4] = modif::nothing;  // Volume fraction. Is default-initialized.
                                   // Should be: staticVariables.
    modified[5] = modif::staticVariables;  // Flag-status. Is adapted according
                                           // to cell-change lists.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists. Read-only.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T rhoDefault;
};

/// Based on the previously computed empty->interface list, initialize flow
/// variables for
///   new interface cells.
/** Input:
 *   - Populations: needed in bulk+0
 *   - Momentum:    needed in bulk+1
 *   - Density:     needed in bulk+1
 *   - Flag-status: needed in bulk+0
 * Output:
 *   - flag-status:   initialized to "interface" on corresponding cells.
 *   - lattice:       initialized from neighbor averages on new interface cells.
 *   - mass:          initialized to zero on new interface cells.
 *   - mass-fraction: initialized to zero on new interface cells.
 *   - momentum
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceIniEmptyToInterfaceNodes2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceIniEmptyToInterfaceNodes2D(
      Dynamics<T, Descriptor> *dynamicsTemplate_,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_)
      : dynamicsTemplate(dynamicsTemplate_), force(force_) {}
  FreeSurfaceIniEmptyToInterfaceNodes2D(
      FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor> const &rhs)
      : dynamicsTemplate(rhs.dynamicsTemplate->clone()), force(rhs.force) {}
  FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor> *operator=(
      FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor> const &rhs) {
    FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(rhs).swap(*this);
    return *this;
  }
  void swap(FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor> &rhs) {
    std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
    std::swap(force, rhs.force);
  }
  ~FreeSurfaceIniEmptyToInterfaceNodes2D() { delete dynamicsTemplate; }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor> *clone() const {
    return new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid. Should be: dataStructure
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::nothing;          // j. Should be: staticVariables.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;  // Volume fraction, read-only. Should be:
                                   // staticVariables
    modified[5] = modif::staticVariables;  // Flag-status, read-only.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists. Read access to
                                           // gasCellToInitializeData.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  Dynamics<T, Descriptor> *dynamicsTemplate;
  Array<T, Descriptor<T>::ExternalField::sizeOfForce>
      force;  // Body force, for initialization of the new interface cell.
};

/// Isolated cells cannot be part of the interface. This data processor spots
/// and removes them.
/** Input:
 *   - Flag-status: needed in bulk+2
 *   - mass:        needed in bulk+1
 *   - density:     needed in bulk+1
 * Output:
 *   - interfaceToFluidNodes:   initialized in bulk+1
 *   - interfaceToEmptyNodes:   initialized in bulk+1
 *   - massExcess list:         initialized in bulk+1
 *   - mass, density, mass-fraction, dynamics, force, momentum, flag: in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceRemoveFalseInterfaceCells2D
    : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceRemoveFalseInterfaceCells2D(T rhoDefault_)
      : rhoDefault(rhoDefault_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor> *clone() const {
    return new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid: Gets NoDynamics when node changes
                                   // to empty. Should be: dataStructure.
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::nothing;          // j. Should be: staticVariables.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] =
        modif::nothing;  // Volume fraction. Should be: staticVariables.
    modified[5] = modif::staticVariables;  // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T rhoDefault;
};

/// Enforce exact mass balance when interface cells become fluid or empty.
/** Input:
 *   - mass-excess list: needed in bulk+1
 *   - Flag-status: needed in bulk+2
 *   - mass:        needed in bulk+2
 *   - density:     needed in bulk+2
 * Output:
 *   - mass, mass-fraction
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceEqualMassExcessReDistribution2D
    : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor> *clone()
      const {
    return new FreeSurfaceEqualMassExcessReDistribution2D(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::dataStructure;    // Fluid.
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::staticVariables;  // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeStatistics2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual FreeSurfaceComputeStatistics2D<T, Descriptor> *clone() const {
    return new FreeSurfaceComputeStatistics2D(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;  // Fluid.
    modified[1] = modif::nothing;  // rhoBar.
    modified[2] = modif::nothing;  // j.
    modified[3] = modif::nothing;  // Mass.
    modified[4] = modif::nothing;  // Volume fraction.
    modified[5] = modif::nothing;  // Flag-status.
    modified[6] = modif::nothing;  // Normal.
    modified[7] = modif::nothing;  // Interface lists.
    modified[8] = modif::nothing;  // Curvature.
    modified[9] = modif::nothing;  // Outside density.
  }
};

template <typename T, template <typename U> class Descriptor>
class InitializeInterfaceLists2D : public BoxProcessingFunctional2D {
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks) {
    PLB_ASSERT(atomicBlocks.size() == 1);

    AtomicContainerBlock2D *containerInterfaceLists =
        dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[0]);
    PLB_ASSERT(containerInterfaceLists);
    InterfaceLists2D<T, Descriptor> *InterfaceLists =
        new InterfaceLists2D<T, Descriptor>;
    containerInterfaceLists->setData(InterfaceLists);
  }
  virtual InitializeInterfaceLists2D<T, Descriptor> *clone() const {
    return new InitializeInterfaceLists2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    // Default-assign potential other parameters present in a multi-fluid
    // system.
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::staticVariables;
  }
};

/// Wrapper for execution of InitializeInterfaceLists2D.
template <typename T, template <typename U> class Descriptor>
void initializeInterfaceLists2D(MultiContainerBlock2D &interfaceListBlock) {
  std::vector<MultiBlock2D *> arg;
  arg.push_back(&interfaceListBlock);
  applyProcessingFunctional(new InitializeInterfaceLists2D<T, Descriptor>,
                            interfaceListBlock.getBoundingBox(), arg);
}

/// Addition of the external forces.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Momentum:      needed in bulk
 * Output:
 *   - Momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddExternalForce2D : public BoxProcessingFunctional2D {
 public:
  FreeSurfaceAddExternalForce2D(T rhoDefault_) : rhoDefault(rhoDefault_) {}
  virtual FreeSurfaceAddExternalForce2D<T, Descriptor> *clone() const {
    return new FreeSurfaceAddExternalForce2D<T, Descriptor>(*this);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T rhoDefault;
};

/// Repelling interface nodes kinematically from the immersed walls.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Density:       needed in bulk
 *   - Momentum:      needed in bulk
 *   - Container:     needed in bulk and envelope
 * Output:
 *   - Density.
 *   - Momentum.
 **/
// TODO: This data processor changes the momentum in the vicinity of the
//       immersed boundary. Maybe it affects strongly the measurement of
//       force and torque on the immersed surface. Needs to be checked.
template <typename T, class VelFunction>
class RepelInterfaceFromImmersedWalls2D : public BoxProcessingFunctional2D {
 public:
  RepelInterfaceFromImmersedWalls2D(VelFunction velFunction_, T rhoDefault_,
                                    bool strongRepelling_)
      : velFunction(velFunction_),
        rhoDefault(rhoDefault_),
        strongRepelling(strongRepelling_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual RepelInterfaceFromImmersedWalls2D<T, VelFunction> *clone() const {
    return new RepelInterfaceFromImmersedWalls2D<T, VelFunction>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Flag
    modified[3] = modif::nothing;  // Container Block with triangle data.
  }
  virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

 private:
  VelFunction velFunction;
  T rhoDefault;
  bool strongRepelling;
};

/// Protect temporarily fluid nodes close to immersed walls from turing to
/// interface.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Container:     needed in bulk and envelope
 * Output:
 *   - Flag-status.
 **/
template <typename T>
class TemporarilyProtectImmersedWalls2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual TemporarilyProtectImmersedWalls2D<T> *clone() const {
    return new TemporarilyProtectImmersedWalls2D<T>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // Flag
    modified[1] = modif::nothing;  // Container Block with triangle data.
  }
  virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
};

/// Remove the above protection.
/** Input:
 *   - Flag-status:   needed in bulk
 * Output:
 *   - Flag-status.
 **/
template <typename T>
class RemoveProtectionFromImmersedWalls2D : public BoxProcessingFunctional2D {
 public:
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual RemoveProtectionFromImmersedWalls2D<T> *clone() const {
    return new RemoveProtectionFromImmersedWalls2D<T>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // Flag
  }
  virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
};

template <typename T, template <typename U> class Descriptor>
struct FreeSurfaceFields2D {
  static const int envelopeWidth;
  static const int smallEnvelopeWidth;
  static const int envelopeWidthForImmersedWalls;

  FreeSurfaceFields2D(
      SparseBlockStructure2D const &blockStructure,
      Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_,
      T contactAngle_,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_)
      : dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        lattice(MultiBlockManagement2D(
                    blockStructure,
                    defaultMultiBlockPolicy2D().getThreadAttribution(),
                    smallEnvelopeWidth),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
                dynamics->clone()),
        helperLists(lattice),
        mass(lattice),
        flag(MultiBlockManagement2D(
                 blockStructure,
                 defaultMultiBlockPolicy2D().getThreadAttribution(),
                 envelopeWidth),
             defaultMultiBlockPolicy2D().getBlockCommunicator(),
             defaultMultiBlockPolicy2D().getCombinedStatistics(),
             defaultMultiBlockPolicy2D().getMultiScalarAccess<int>()),
        volumeFraction((MultiBlock2D &)flag),
        curvature(MultiBlockManagement2D(
                      blockStructure,
                      defaultMultiBlockPolicy2D().getThreadAttribution(),
                      envelopeWidth),
                  defaultMultiBlockPolicy2D().getBlockCommunicator(),
                  defaultMultiBlockPolicy2D().getCombinedStatistics(),
                  defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()),
        outsideDensity((MultiBlock2D &)curvature),
        rhoBar(MultiBlockManagement2D(
                   blockStructure,
                   defaultMultiBlockPolicy2D().getThreadAttribution(),
                   smallEnvelopeWidth),
               defaultMultiBlockPolicy2D().getBlockCommunicator(),
               defaultMultiBlockPolicy2D().getCombinedStatistics(),
               defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()),
        j(MultiBlockManagement2D(
              blockStructure,
              defaultMultiBlockPolicy2D().getThreadAttribution(),
              smallEnvelopeWidth),
          defaultMultiBlockPolicy2D().getBlockCommunicator(),
          defaultMultiBlockPolicy2D().getCombinedStatistics(),
          defaultMultiBlockPolicy2D().getMultiTensorAccess<T, 2>()),
        normal((MultiBlock2D &)curvature),
        container(0) {
    incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
    if (incompressibleModel) {
      // Incompressible: rho0=1
      PLB_ASSERT(util::isOne(rhoDefault));
    }
#endif

    useSurfaceTension = !util::isZero(surfaceTension);

    freeSurfaceArgs = aggregateFreeSurfaceParams(
        lattice, rhoBar, j, mass, volumeFraction, flag, normal, helperLists,
        curvature, outsideDensity);

    initializeInterfaceLists2D<T, Descriptor>(helperLists);
    lattice.periodicity().toggleAll(true);
    mass.periodicity().toggleAll(true);
    flag.periodicity().toggleAll(true);
    volumeFraction.periodicity().toggleAll(true);
    curvature.periodicity().toggleAll(true);
    outsideDensity.periodicity().toggleAll(true);
    rhoBar.periodicity().toggleAll(true);
    j.periodicity().toggleAll(true);
    normal.periodicity().toggleAll(true);
    // setToConstant(flag, flag.getBoundingBox(),
    // (int)freeSurfaceFlag2D::empty); setToConstant(outsideDensity,
    // outsideDensity.getBoundingBox(), rhoDefault);
    rhoBarJparam.push_back(&lattice);
    rhoBarJparam.push_back(&rhoBar);
    rhoBarJparam.push_back(&j);

    lattice.internalStatSubscription().subscribeSum();  // Total mass.
    lattice.internalStatSubscription().subscribeSum();  // Lost mass.
    lattice.internalStatSubscription()
        .subscribeIntSum();  // Num interface cells.

    freeSurfaceDataProcessors();
  }

  // TODO: The default argument repelInterface is to test different
  //       methods of repelling bubbles (and droplets) from immersed boundaries.
  //       0: no repelling
  //       1: use rhoBar-j kinematic repelling
  //       2: use repelling by changing the "temporarilyProtect" flag
  //       3: use rhoBar-j kinematic repelling, but a very strong one
  //       When tests are over, this argument is meant to be disposed of.
  template <class VelFunction>
  FreeSurfaceFields2D(
      SparseBlockStructure2D const &blockStructure,
      Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_,
      T contactAngle_,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
      plint numIBIterations, std::vector<Array<T, 2>> const &vertices,
      std::vector<T> const &areas, std::vector<int> const &flags,
      VelFunction velFunction, int repelInterface = 0)
      : dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        lattice(MultiBlockManagement2D(
                    blockStructure,
                    defaultMultiBlockPolicy2D().getThreadAttribution(),
                    smallEnvelopeWidth),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
                dynamics->clone()),
        helperLists(lattice),
        mass(lattice),
        flag(MultiBlockManagement2D(
                 blockStructure,
                 defaultMultiBlockPolicy2D().getThreadAttribution(),
                 envelopeWidthForImmersedWalls),
             defaultMultiBlockPolicy2D().getBlockCommunicator(),
             defaultMultiBlockPolicy2D().getCombinedStatistics(),
             defaultMultiBlockPolicy2D().getMultiScalarAccess<int>()),
        volumeFraction((MultiBlock2D &)flag),
        curvature(MultiBlockManagement2D(
                      blockStructure,
                      defaultMultiBlockPolicy2D().getThreadAttribution(),
                      envelopeWidth),
                  defaultMultiBlockPolicy2D().getBlockCommunicator(),
                  defaultMultiBlockPolicy2D().getCombinedStatistics(),
                  defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()),
        outsideDensity((MultiBlock2D &)curvature),
        rhoBar(MultiBlockManagement2D(
                   blockStructure,
                   defaultMultiBlockPolicy2D().getThreadAttribution(),
                   envelopeWidthForImmersedWalls),
               defaultMultiBlockPolicy2D().getBlockCommunicator(),
               defaultMultiBlockPolicy2D().getCombinedStatistics(),
               defaultMultiBlockPolicy2D().getMultiScalarAccess<T>()),
        j(MultiBlockManagement2D(
              blockStructure,
              defaultMultiBlockPolicy2D().getThreadAttribution(),
              envelopeWidthForImmersedWalls),
          defaultMultiBlockPolicy2D().getBlockCommunicator(),
          defaultMultiBlockPolicy2D().getCombinedStatistics(),
          defaultMultiBlockPolicy2D().getMultiTensorAccess<T, 2>()),
        normal((MultiBlock2D &)curvature) {
    container = new MultiContainerBlock2D((MultiBlock2D &)rhoBar);
    PLB_ASSERT(container);

    incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
    if (incompressibleModel) {
      // Incompressible: rho0=1
      PLB_ASSERT(util::isOne(rhoDefault));
    }
#endif

    useSurfaceTension = !util::isZero(surfaceTension);

    freeSurfaceArgs = aggregateFreeSurfaceParams(
        lattice, rhoBar, j, mass, volumeFraction, flag, normal, helperLists,
        curvature, outsideDensity);

    initializeInterfaceLists2D<T, Descriptor>(helperLists);
    lattice.periodicity().toggleAll(true);
    mass.periodicity().toggleAll(true);
    flag.periodicity().toggleAll(true);
    volumeFraction.periodicity().toggleAll(true);
    curvature.periodicity().toggleAll(true);
    outsideDensity.periodicity().toggleAll(true);
    rhoBar.periodicity().toggleAll(true);
    j.periodicity().toggleAll(true);
    normal.periodicity().toggleAll(true);
    // setToConstant(flag, flag.getBoundingBox(),
    // (int)freeSurfaceFlag2D::empty); setToConstant(outsideDensity,
    // outsideDensity.getBoundingBox(), rhoDefault);
    rhoBarJparam.push_back(&lattice);
    rhoBarJparam.push_back(&rhoBar);
    rhoBarJparam.push_back(&j);

    lattice.internalStatSubscription().subscribeSum();  // Total mass.
    lattice.internalStatSubscription().subscribeSum();  // Lost mass.
    lattice.internalStatSubscription()
        .subscribeIntSum();  // Num interface cells.

    freeSurfaceDataProcessorsForImmersedWalls(
        numIBIterations, vertices, areas, flags, velFunction, repelInterface);
  }

  FreeSurfaceFields2D(FreeSurfaceFields2D<T, Descriptor> const &rhs)
      : dynamics(rhs.dynamics->clone()),
        incompressibleModel(rhs.incompressibleModel),
        rhoDefault(rhs.rhoDefault),
        surfaceTension(rhs.surfaceTension),
        contactAngle(rhs.contactAngle),
        useSurfaceTension(rhs.useSurfaceTension),
        force(rhs.force),
        lattice(rhs.lattice),
        helperLists(rhs.helperLists),
        mass(rhs.mass),
        flag(rhs.flag),
        volumeFraction(rhs.volumeFraction),
        curvature(rhs.curvature),
        outsideDensity(rhs.outsideDensity),
        rhoBar(rhs.rhoBar),
        j(rhs.j),
        normal(rhs.normal),
        rhoBarJparam(rhs.rhoBarJparam),
        freeSurfaceArgs(rhs.freeSurfaceArgs) {
    container = 0;
    if (rhs.container) {
      container = rhs.container->clone();
      PLB_ASSERT(container);
    }
  }

  void swap(FreeSurfaceFields2D<T, Descriptor> &rhs) {
    std::swap(dynamics, rhs.dynamics);
    std::swap(incompressibleModel, rhs.incompressibleModel);
    std::swap(rhoDefault, rhs.rhoDefault);
    std::swap(surfaceTension, rhs.surfaceTension);
    std::swap(contactAngle, rhs.contactAngle);
    std::swap(useSurfaceTension, rhs.useSurfaceTension);
    std::swap(force, rhs.force);
    std::swap(lattice, rhs.lattice);
    std::swap(helperLists, rhs.helperLists);
    std::swap(mass, rhs.mass);
    std::swap(flag, rhs.flag);
    std::swap(volumeFraction, rhs.volumeFraction);
    std::swap(curvature, rhs.curvature);
    std::swap(outsideDensity, rhs.outsideDensity);
    std::swap(rhoBar, rhs.rhoBar);
    std::swap(j, rhs.j);
    std::swap(normal, rhs.normal);
    std::swap(container, rhs.container);
    std::swap(rhoBarJparam, rhs.rhoBarJparam);
    std::swap(freeSurfaceArgs, rhs.freeSurfaceArgs);
  }

  FreeSurfaceFields2D<T, Descriptor> &operator=(
      FreeSurfaceFields2D<T, Descriptor> const &rhs) {
    FreeSurfaceFields2D<T, Descriptor>(rhs).swap(*this);
    return *this;
  }

  FreeSurfaceFields2D<T, Descriptor> *clone() const {
    return new FreeSurfaceFields2D<T, Descriptor>(*this);
  }

  ~FreeSurfaceFields2D() {
    delete dynamics;
    delete container;
  }

  void periodicityToggle(plint direction, bool periodic) {
    PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

    lattice.periodicity().toggle(direction, periodic);
    mass.periodicity().toggle(direction, periodic);
    flag.periodicity().toggle(direction, periodic);
    volumeFraction.periodicity().toggle(direction, periodic);
    curvature.periodicity().toggle(direction, periodic);
    outsideDensity.periodicity().toggle(direction, periodic);
    rhoBar.periodicity().toggle(direction, periodic);
    j.periodicity().toggle(direction, periodic);
    normal.periodicity().toggle(direction, periodic);
    if (container) {
      container->periodicity().toggle(direction, periodic);
    }
  }

  void periodicityToggleAll(bool periodic) {
    lattice.periodicity().toggleAll(periodic);
    mass.periodicity().toggleAll(periodic);
    flag.periodicity().toggleAll(periodic);
    volumeFraction.periodicity().toggleAll(periodic);
    curvature.periodicity().toggleAll(periodic);
    outsideDensity.periodicity().toggleAll(periodic);
    rhoBar.periodicity().toggleAll(periodic);
    j.periodicity().toggleAll(periodic);
    normal.periodicity().toggleAll(periodic);
    if (container) {
      container->periodicity().toggleAll(periodic);
    }
  }

  void defaultInitialize(bool useConstRho = true, bool useZeroMomentum = true,
                         bool initializeCell = true) {
    applyProcessingFunctional(new DefaultInitializeFreeSurface2D<T, Descriptor>(
                                  dynamics->clone(), force, rhoDefault,
                                  useConstRho, useZeroMomentum, initializeCell),
                              lattice.getBoundingBox(), freeSurfaceArgs);
  }

  void partiallyDefaultInitialize(bool useConstRho = true,
                                  bool useZeroMomentum = true,
                                  bool initializeCell = true) {
    applyProcessingFunctional(
        new PartiallyDefaultInitializeFreeSurface2D<T, Descriptor>(
            dynamics->clone(), force, rhoDefault, useConstRho, useZeroMomentum,
            initializeCell),
        lattice.getBoundingBox(), freeSurfaceArgs);
  }

  void freeSurfaceDataProcessors() {
    plint pl;  // Processor level.

    /***** Initial level ******/
    pl = 0;

    integrateProcessingFunctional(
        new ExternalRhoJcollideAndStream2D<T, Descriptor>,
        lattice.getBoundingBox(), rhoBarJparam, pl);

    integrateProcessingFunctional(
        new FreeSurfaceComputeNormals2D<T, Descriptor>,
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    if (useSurfaceTension) {
      integrateProcessingFunctional(
          new FreeSurfaceComputeCurvature2D<T, Descriptor>(contactAngle),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);

      // To change to the curvature calculation with height functions, uncomment
      // the next data processor and comment out the two previous ones. If only
      // the next data processor is used and there is no surface tension, the
      // normals are not computed at all. Be careful if you intent to use the
      // normals and do not have the surface tension algorithm enabled.
      // integrateProcessingFunctional (
      //        new FreeSurfaceGeometry2D<T,Descriptor>(contactAngle),
      //        lattice.getBoundingBox(), freeSurfaceArgs, pl );
    }

    integrateProcessingFunctional(new FreeSurfaceMassChange2D<T, Descriptor>,
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    integrateProcessingFunctional(new FreeSurfaceCompletion2D<T, Descriptor>,
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    integrateProcessingFunctional(
        new FreeSurfaceMacroscopic2D<T, Descriptor>(incompressibleModel),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    if (useSurfaceTension) {
      integrateProcessingFunctional(
          new FreeSurfaceAddSurfaceTension2D<T, Descriptor>(
              surfaceTension, incompressibleModel),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);
    }

    integrateProcessingFunctional(new FreeSurfaceStabilize2D<T, Descriptor>(),
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoDefault),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
            dynamics->clone(), force),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoDefault),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceComputeStatistics2D<T, Descriptor>,
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    bool useForce = !util::isZero(norm(force));
    if (useForce) {
      integrateProcessingFunctional(
          new FreeSurfaceAddExternalForce2D<T, Descriptor>(rhoDefault),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);
    }
  }

  template <class VelFunction>
  void freeSurfaceDataProcessorsForImmersedWalls(
      plint numIBIterations, std::vector<Array<T, 2>> const &vertices,
      std::vector<T> const &areas, std::vector<int> const &flags,
      VelFunction velFunction, int repelInterface) {
    plint pl;  // Processor level.

    /***** Initial level ******/
    pl = 0;

    integrateProcessingFunctional(
        new ExternalRhoJcollideAndStream2D<T, Descriptor>,
        lattice.getBoundingBox(), rhoBarJparam, pl);

    integrateProcessingFunctional(
        new FreeSurfaceComputeNormals2D<T, Descriptor>,
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    if (useSurfaceTension) {
      integrateProcessingFunctional(
          new FreeSurfaceComputeCurvature2D<T, Descriptor>(contactAngle),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);

      // To change to the curvature calculation with height functions, uncomment
      // the next data processor and comment out the two previous ones. If only
      // the next data processor is used and there is no surface tension, the
      // normals are not computed at all. Be careful if you intent to use the
      // normals and do not have the surface tension algorithm enabled.
      // integrateProcessingFunctional (
      //        new FreeSurfaceGeometry2D<T,Descriptor>(contactAngle),
      //        lattice.getBoundingBox(), freeSurfaceArgs, pl );
    }

    integrateProcessingFunctional(new FreeSurfaceMassChange2D<T, Descriptor>,
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    integrateProcessingFunctional(new FreeSurfaceCompletion2D<T, Descriptor>,
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    integrateProcessingFunctional(
        new FreeSurfaceMacroscopic2D<T, Descriptor>(incompressibleModel),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    if (useSurfaceTension) {
      integrateProcessingFunctional(
          new FreeSurfaceAddSurfaceTension2D<T, Descriptor>(
              surfaceTension, incompressibleModel),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);
    }

    integrateProcessingFunctional(new FreeSurfaceStabilize2D<T, Descriptor>(),
                                  lattice.getBoundingBox(), freeSurfaceArgs,
                                  pl);

    std::vector<MultiBlock2D *> immersedWallDataArgs;
    immersedWallDataArgs.push_back(container);
    integrateProcessingFunctional(
        new InstantiateImmersedWallDataWithIndexedTagging2D<T>(vertices, areas,
                                                               flags),
        container->getBoundingBox(), lattice, immersedWallDataArgs, pl);

    /***** New level ******/
    if (repelInterface == 2) {
      pl++;

      std::vector<MultiBlock2D *> tmpProtectionArgs;
      tmpProtectionArgs.push_back(&flag);
      tmpProtectionArgs.push_back(container);
      integrateProcessingFunctional(new TemporarilyProtectImmersedWalls2D<T>(),
                                    flag.getBoundingBox(), lattice,
                                    tmpProtectionArgs, pl);
    }

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoDefault),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
            dynamics->clone(), force),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoDefault),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceComputeStatistics2D<T, Descriptor>,
        lattice.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/

    for (int i = 0; i < numIBIterations; i++) {
      pl++;
      T tau = (T)1 / dynamics->getOmega();
      std::vector<MultiBlock2D *> args;
      args.push_back(&rhoBar);
      args.push_back(&j);
      args.push_back(container);
      integrateProcessingFunctional(
          new IndexedInamuroIteration2D<T, VelFunction>(velFunction, tau,
                                                        incompressibleModel),
          rhoBar.getBoundingBox(), lattice, args, pl);
    }

    if (repelInterface == 1 || repelInterface == 3) {
      bool strongRepelling = repelInterface == 1 ? false : true;

      // TODO: This data processor changes the momentum in the vicinity of the
      //       immersed boundary. Maybe it affects strongly the measurement of
      //       force and torque on the immersed surface. Needs to be checked.
      pl++;
      std::vector<MultiBlock2D *> args;
      args.push_back(&rhoBar);
      args.push_back(&j);
      args.push_back(&flag);
      args.push_back(container);
      integrateProcessingFunctional(
          new RepelInterfaceFromImmersedWalls2D<T, VelFunction>(
              velFunction, rhoDefault, strongRepelling),
          rhoBar.getBoundingBox(), lattice, args, pl);
    }

    /***** New level ******/
    if (repelInterface == 2) {
      pl++;

      std::vector<MultiBlock2D *> rmProtectionArgs;
      rmProtectionArgs.push_back(&flag);
      integrateProcessingFunctional(
          new RemoveProtectionFromImmersedWalls2D<T>(), flag.getBoundingBox(),
          lattice, rmProtectionArgs, pl);
    }

    /***** New level ******/
    bool useForce = !util::isZero(norm(force));
    if (useForce) {
      pl++;

      integrateProcessingFunctional(
          new FreeSurfaceAddExternalForce2D<T, Descriptor>(rhoDefault),
          lattice.getBoundingBox(), freeSurfaceArgs, pl);
    }
  }

  void appendBlocksToCheckpointVector(
      std::vector<MultiBlock2D *> &checkpointBlocks) {
    checkpointBlocks.push_back(&lattice);
    checkpointBlocks.push_back(&mass);
    checkpointBlocks.push_back(&flag);
    checkpointBlocks.push_back(&volumeFraction);
    checkpointBlocks.push_back(&outsideDensity);
    checkpointBlocks.push_back(&rhoBar);
    checkpointBlocks.push_back(&j);
  }

  Dynamics<T, Descriptor> *dynamics;
  bool incompressibleModel;
  T rhoDefault;
  T surfaceTension;
  T contactAngle;
  bool useSurfaceTension;
  Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
  MultiBlockLattice2D<T, Descriptor> lattice;
  MultiContainerBlock2D helperLists;
  MultiScalarField2D<T> mass;
  MultiScalarField2D<int> flag;
  MultiScalarField2D<T> volumeFraction;
  MultiScalarField2D<T> curvature;
  MultiScalarField2D<T> outsideDensity;
  MultiScalarField2D<T> rhoBar;
  MultiTensorField2D<T, 2> j;
  MultiTensorField2D<T, 2> normal;
  MultiContainerBlock2D *container;
  std::vector<MultiBlock2D *> rhoBarJparam;
  std::vector<MultiBlock2D *> freeSurfaceArgs;
};

template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields2D<T, Descriptor>::envelopeWidth =
    4;  // Mehdi: 2? Necessary when we use height functions to compute the
        // curvature,
// or when double smoothing is used at the data processor that
// computes the normals from the volume fraction.
// template<typename T, template<typename U> class Descriptor>
// const int FreeSurfaceFields2D<T,Descriptor>::envelopeWidth = 4; // Necessary
// when we use height functions to compute the curvature and use the old contact
// angle algorithm.
template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields2D<T, Descriptor>::smallEnvelopeWidth =
    4;  // Mehdi: Changed to 3 for disjoining pressure

template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields2D<T, Descriptor>::envelopeWidthForImmersedWalls = 4;

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceSetup2D {
 public:
  typedef T (*ContactAngleFunction)(
      T x, T y);  // Returns the contact angle in degrees.
 public:
  FreeSurfaceSetup2D(Group2D &group_, Dynamics<T, Descriptor> const &dynamics_,
                     T rhoDefault_, T surfaceTension_, T contactAngle_,
                     Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
                     std::string fluidname_ = "fluid1",
                     std::string fsPrefix_ = "",
                     std::string rhoBarJprefix_ = "",
                     bool hasImmersedWalls_ = false)
      : group(group_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_) {
    freeSurfaceArgs.clear();
    rhoBarJparam.clear();
    incompressibleModel = dynamics.velIsJ();
#ifdef PLB_DEBUG
    if (incompressibleModel) {
      // Incompressible: rho0=1
      PLB_ASSERT(util::isOne(rhoDefault));
    }
#endif
  }

  FreeSurfaceSetup2D(Group2D &group_, Dynamics<T, Descriptor> const &dynamics_,
                     T rhoDefault_, T surfaceTension_,
                     ContactAngleFunction contactAngleFunction_,
                     Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
                     std::string fluidname_ = "fluid1",
                     std::string fsPrefix_ = "",
                     std::string rhoBarJprefix_ = "",
                     bool hasImmersedWalls_ = false)
      : group(group_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_) {
    freeSurfaceArgs.clear();
    rhoBarJparam.clear();
    incompressibleModel = dynamics.velIsJ();
#ifdef PLB_DEBUG
    if (incompressibleModel) {
      // Incompressible: rho0=1
      PLB_ASSERT(util::isOne(rhoDefault));
    }
#endif
  }

  void periodicityToggle(plint direction, bool periodic) {
    PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

    group.get(fluidname).periodicity().toggle(direction, periodic);
    group.get(fsPrefix + "mass").periodicity().toggle(direction, periodic);
    group.get(fsPrefix + "flag").periodicity().toggle(direction, periodic);
    group.get(fsPrefix + "volumeFraction")
        .periodicity()
        .toggle(direction, periodic);
    group.get(fsPrefix + "curvature").periodicity().toggle(direction, periodic);
    group.get(fsPrefix + "outsideDensity")
        .periodicity()
        .toggle(direction, periodic);
    group.get(rhoBarJprefix + "rhoBar")
        .periodicity()
        .toggle(direction, periodic);
    group.get(rhoBarJprefix + "j").periodicity().toggle(direction, periodic);
    group.get(fsPrefix + "normal").periodicity().toggle(direction, periodic);
  }

  void periodicityToggleAll(bool periodic) {
    group.get(fluidname).periodicity().toggleAll(periodic);
    group.get(fsPrefix + "mass").periodicity().toggleAll(periodic);
    group.get(fsPrefix + "flag").periodicity().toggleAll(periodic);
    group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(periodic);
    group.get(fsPrefix + "curvature").periodicity().toggleAll(periodic);
    group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(periodic);
    group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(periodic);
    group.get(rhoBarJprefix + "j").periodicity().toggleAll(periodic);
    group.get(fsPrefix + "normal").periodicity().toggleAll(periodic);
  }

  void createFreeSurfaceFields() {
    const int envelopeWidth = FreeSurfaceFields2D<T, Descriptor>::envelopeWidth;
    const int smallEnvelopeWidth =
        FreeSurfaceFields2D<T, Descriptor>::smallEnvelopeWidth;
    const int smallOrLargeEnvelopeWidth =
        hasImmersedWalls
            ? FreeSurfaceFields2D<T, Descriptor>::envelopeWidthForImmersedWalls
            : smallEnvelopeWidth;
    const int mediumOrLargeEnvelopeWidth =
        hasImmersedWalls
            ? FreeSurfaceFields2D<T, Descriptor>::envelopeWidthForImmersedWalls
            : envelopeWidth;
#ifdef PLB_DEBUG
    bool hasFluid = fieldExists(fluidname, smallEnvelopeWidth);
#endif
    PLB_ASSERT(hasFluid);
    group.get(fluidname).periodicity().toggleAll(true);
    if (!fieldExists(fsPrefix + "helperLists", 0)) {
      group.generateContainer(fsPrefix + "helperLists");
    }
    if (!fieldExists(fsPrefix + "mass", smallEnvelopeWidth)) {
      group.generateScalar<T>(fsPrefix + "mass", smallEnvelopeWidth);
    }
    group.get(fsPrefix + "mass").periodicity().toggleAll(true);
    if (!fieldExists("flag", mediumOrLargeEnvelopeWidth)) {
      group.generateScalar<int>(fsPrefix + "flag", mediumOrLargeEnvelopeWidth);
    }
    group.get(fsPrefix + "flag").periodicity().toggleAll(true);
    if (!fieldExists(fsPrefix + "volumeFraction", envelopeWidth)) {
      group.generateScalar<T>(fsPrefix + "volumeFraction", envelopeWidth);
    }
    group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(true);
    if (!fieldExists(fsPrefix + "curvature", envelopeWidth)) {
      group.generateScalar<T>(fsPrefix + "curvature", envelopeWidth);
    }
    group.get(fsPrefix + "curvature").periodicity().toggleAll(true);
    if (!fieldExists(fsPrefix + "outsideDensity", envelopeWidth)) {
      group.generateScalar<T>(fsPrefix + "outsideDensity", envelopeWidth);
    }
    group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(true);
    if (!fieldExists(rhoBarJprefix + "rhoBar", smallOrLargeEnvelopeWidth)) {
      group.generateScalar<T>(rhoBarJprefix + "rhoBar",
                              smallOrLargeEnvelopeWidth);
    }
    group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(true);
    if (!fieldExists(rhoBarJprefix + "j", smallOrLargeEnvelopeWidth)) {
      group.generateTensor<T, 2>(rhoBarJprefix + "j",
                                 smallOrLargeEnvelopeWidth);
    }
    group.get(rhoBarJprefix + "j").periodicity().toggleAll(true);
    if (!fieldExists(fsPrefix + "normal", envelopeWidth)) {
      group.generateTensor<T, 2>(fsPrefix + "normal", envelopeWidth);
    }
    group.get(fsPrefix + "normal").periodicity().toggleAll(true);
    if (hasImmersedWalls) {
      if (!fieldExists(fsPrefix + "ibm_container", smallOrLargeEnvelopeWidth)) {
        group.generateContainer(fsPrefix + "ibm_container",
                                smallOrLargeEnvelopeWidth);
      }
    }

    freeSurfaceArgs.push_back(&group.get(fluidname));
    freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "rhoBar"));
    freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "j"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "mass"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "volumeFraction"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "flag"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "normal"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "helperLists"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "curvature"));
    freeSurfaceArgs.push_back(&group.get(fsPrefix + "outsideDensity"));

    rhoBarJparam.push_back(&group.get(fluidname));
    rhoBarJparam.push_back(&group.get(rhoBarJprefix + "rhoBar"));
    rhoBarJparam.push_back(&group.get(rhoBarJprefix + "j"));
  }

  Actions2D freeSurfaceActions(std::map<std::string, plint> *ids = 0) {
    Actions2D actions;

    bool useSurfaceTension = !util::isZero(surfaceTension);

    plint id = -1;
    id = actions.addBlock(group.get(fluidname));  // 0
    if (ids != 0) (*ids)[fluidname] = id;
    id = actions.addBlock(group.get(rhoBarJprefix + "rhoBar"));  // 1
    if (ids != 0) (*ids)[rhoBarJprefix + "rhoBar"] = id;
    id = actions.addBlock(group.get(rhoBarJprefix + "j"));  // 2
    if (ids != 0) (*ids)[rhoBarJprefix + "j"] = id;
    id = actions.addBlock(group.get(fsPrefix + "mass"));  // 3
    if (ids != 0) (*ids)[fsPrefix + "mass"] = id;
    id = actions.addBlock(group.get(fsPrefix + "volumeFraction"));  // 4
    if (ids != 0) (*ids)[fsPrefix + "volumeFraction"] = id;
    id = actions.addBlock(group.get(fsPrefix + "flag"));  // 5
    if (ids != 0) (*ids)[fsPrefix + "flag"] = id;
    id = actions.addBlock(group.get(fsPrefix + "normal"));  // 6
    if (ids != 0) (*ids)[fsPrefix + "normal"] = id;
    id = actions.addBlock(group.get(fsPrefix + "helperLists"));  // 7
    if (ids != 0) (*ids)[fsPrefix + "helperLists"] = id;
    id = actions.addBlock(group.get(fsPrefix + "curvature"));  // 8
    if (ids != 0) (*ids)[fsPrefix + "curvature"] = id;
    id = actions.addBlock(group.get(fsPrefix + "outsideDensity"));  // 9
    if (ids != 0) (*ids)[fsPrefix + "outsideDensity"] = id;

    std::vector<plint> freeSurfaceBlocks;
    for (plint i = 0; i < 10; ++i) freeSurfaceBlocks.push_back(i);

    initializeInterfaceLists2D<T, Descriptor>(
        group.getContainer(fsPrefix + "helperLists"));
    // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
    // (int)freeSurfaceFlag2D::empty);
    // setToConstant(group.getScalar<T>("outsideDensity"),
    // group.getBoundingBox(), rhoDefault);

    std::vector<plint> rhoBarJblocks;
    for (plint i = 0; i < 3; ++i) rhoBarJblocks.push_back(i);  // mehdi 2?

    group.get(fluidname)
        .internalStatSubscription()
        .subscribeSum();  // Total mass.
    group.get(fluidname)
        .internalStatSubscription()
        .subscribeSum();  // Lost mass.
    group.get(fluidname)
        .internalStatSubscription()
        .subscribeIntSum();  // Num interface cells.

    actions.addProcessor(new ExternalRhoJcollideAndStream2D<T, Descriptor>,
                         rhoBarJblocks, group.getBoundingBox());
    actions.addProcessor(new FreeSurfaceComputeNormals2D<T, Descriptor>,
                         freeSurfaceBlocks, group.getBoundingBox());
    actions.addCommunication(0, modif::staticVariables);  // fluid1
    actions.addCommunication(6, modif::staticVariables);  // normal

    if (useSurfaceTension) {
      if (contactAngleFunction == 0) {
        actions.addProcessor(
            new FreeSurfaceComputeCurvature2D<T, Descriptor>(contactAngle),
            freeSurfaceBlocks, group.getBoundingBox());
      } else {
        actions.addProcessor(new FreeSurfaceComputeCurvature2D<T, Descriptor>(
                                 contactAngleFunction),
                             freeSurfaceBlocks, group.getBoundingBox());
      }
    }
    actions.addProcessor(new FreeSurfaceMassChange2D<T, Descriptor>,
                         freeSurfaceBlocks, group.getBoundingBox());
    actions.addProcessor(new FreeSurfaceCompletion2D<T, Descriptor>,
                         freeSurfaceBlocks, group.getBoundingBox());
    actions.addProcessor(
        new FreeSurfaceMacroscopic2D<T, Descriptor>(incompressibleModel),
        freeSurfaceBlocks, group.getBoundingBox());
    if (useSurfaceTension) {
      actions.addProcessor(new FreeSurfaceAddSurfaceTension2D<T, Descriptor>(
                               surfaceTension, incompressibleModel),
                           freeSurfaceBlocks, group.getBoundingBox());
    }
    actions.addProcessor(new FreeSurfaceStabilize2D<T, Descriptor>(),
                         freeSurfaceBlocks, group.getBoundingBox());
    actions.addCommunication(0, modif::staticVariables);  // fluid1
    actions.addCommunication(1, modif::staticVariables);  // rhoBar.
    actions.addCommunication(2, modif::staticVariables);  // j.
    actions.addCommunication(3, modif::staticVariables);  // Mass.
    actions.addCommunication(4, modif::staticVariables);  // Volume fraction.
    if (useSurfaceTension) {
      actions.addCommunication(8, modif::staticVariables);  // Curvature.
    }

    actions.addProcessor(
        new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
        freeSurfaceBlocks, group.getBoundingBox());
    actions.addProcessor(
        new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoDefault),
        freeSurfaceBlocks, group.getBoundingBox());
    actions.addProcessor(
        new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
            dynamics.clone(), force),
        freeSurfaceBlocks, group.getBoundingBox());

    actions.addCommunication(1, modif::staticVariables);  // rhoBar.
    actions.addCommunication(3, modif::staticVariables);  // Mass.
    actions.addCommunication(5, modif::staticVariables);  // Flag-status.
    actions.addCommunication(7, modif::staticVariables);  // Interface-lists.

    actions.addProcessor(
        new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoDefault),
        freeSurfaceBlocks, group.getBoundingBox());
    actions.addCommunication(1, modif::staticVariables);  // rhoBar.
    actions.addCommunication(3, modif::staticVariables);  // Mass.
    actions.addCommunication(5, modif::staticVariables);  // Flag-status.

    actions.addProcessor(
        new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>,
        freeSurfaceBlocks, group.getBoundingBox());
    actions.addProcessor(new FreeSurfaceComputeStatistics2D<T, Descriptor>,
                         freeSurfaceBlocks, group.getBoundingBox());
    bool useForce = !util::isZero(norm(force));
    if (useForce) {
      actions.addProcessor(
          new FreeSurfaceAddExternalForce2D<T, Descriptor>(rhoDefault),
          freeSurfaceBlocks, group.getBoundingBox());
    }
    actions.addCommunication(0, modif::dataStructure);    // Fluid.
    actions.addCommunication(1, modif::staticVariables);  // rhoBar.
    actions.addCommunication(2, modif::staticVariables);  // j.
    actions.addCommunication(3, modif::staticVariables);  // Mass.
    actions.addCommunication(4, modif::staticVariables);  // Volume fraction.

    actions.addEvaluateStats(0);                 // Fluid.
    actions.addIncrementTime<T, Descriptor>(0);  // Fluid.

    return actions;
  }

  template <class VelFunction>
  Actions2D immersedWallActions(plint numIBIterations,
                                MultiContainerBlock2D &container,
                                std::vector<Array<T, 2>> const &vertices,
                                std::vector<T> const &areas,
                                std::vector<int> const &flags,
                                VelFunction velFunction, bool strongRepelling) {
    Actions2D actions;

    plint containerID = actions.addBlock(group.get(fsPrefix + "ibm_container"));
    plint rhoBarID = actions.addBlock(group.get(rhoBarJprefix + "rhoBar"));
    plint jID = actions.addBlock(group.get(rhoBarJprefix + "j"));
    plint flagID = actions.addBlock(group.get(fsPrefix + "flag"));

    actions.addProcessor(new InstantiateImmersedWallDataWithIndexedTagging2D<T>(
                             vertices, areas, flags),
                         containerID, group.getBoundingBox());
    T tau = (T)1 / dynamics.getOmega();
    for (int i = 0; i < numIBIterations; i++) {
      actions.addProcessor(new IndexedInamuroIteration2D<T, VelFunction>(
                               velFunction, tau, incompressibleModel),
                           rhoBarID, jID, containerID, group.getBoundingBox());
      actions.addCommunication(jID, modif::staticVariables);
    }
    actions.addProcessor(new RepelInterfaceFromImmersedWalls2D<T, VelFunction>(
                             velFunction, rhoDefault, strongRepelling),
                         rhoBarID, jID, flagID, containerID,
                         group.getBoundingBox());
    actions.addCommunication(rhoBarID, modif::staticVariables);
    actions.addCommunication(jID, modif::staticVariables);

    return actions;
  }

  // If surfaceTensionField != 0, then the
  // FreeSurfaceAddSurfaceTensionFromScalarField2D processor is integrated
  // instead of the FreeSurfaceAddSurfaceTension2D one. Contrary to the Palabos
  // convention, this function does NOT take ownership of the
  // surfaceTensionField despite the fact that a pointer is passed. The caller
  // is responsible for its memory management.
  void createFreeSurfaceProcessors(
      plint initialProcessorLevel = 0,
      MultiScalarField2D<T> *surfaceTensionField = 0) {
    PLB_ASSERT(initialProcessorLevel >= 0);

    bool useSurfaceTension =
        (surfaceTensionField != 0 || !util::isZero(surfaceTension));

    initializeInterfaceLists2D<T, Descriptor>(
        group.getContainer(fsPrefix + "helperLists"));
    // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
    // (int)freeSurfaceFlag2D::empty);
    // setToConstant(group.getScalar<T>("outsideDensity"),
    // group.getBoundingBox(), rhoDefault);

    group.get(fluidname)
        .internalStatSubscription()
        .subscribeSum();  // Total mass.
    group.get(fluidname)
        .internalStatSubscription()
        .subscribeSum();  // Lost mass.
    group.get(fluidname)
        .internalStatSubscription()
        .subscribeIntSum();  // Num interface cells.

    plint pl = initialProcessorLevel;  // Processor level.

    /***** Initial level ******/
    if (pl == 0) {
      integrateProcessingFunctional(
          new ExternalRhoJcollideAndStream2D<T, Descriptor>,
          group.getBoundingBox(), rhoBarJparam, pl);
    }

    integrateProcessingFunctional(
        new FreeSurfaceComputeNormals2D<T, Descriptor>, group.getBoundingBox(),
        freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    if (useSurfaceTension) {
      if (contactAngleFunction == 0) {
        integrateProcessingFunctional(
            new FreeSurfaceComputeCurvature2D<T, Descriptor>(contactAngle),
            group.getBoundingBox(), freeSurfaceArgs, pl);
      } else {
        integrateProcessingFunctional(
            new FreeSurfaceComputeCurvature2D<T, Descriptor>(
                contactAngleFunction),
            group.getBoundingBox(), freeSurfaceArgs, pl);
      }

      // To change to the curvature calculation with height functions, uncomment
      // the next data processor and comment out the two previous ones. If only
      // the next data processor is used and there is no surface tension, the
      // normals are not computed at all. Be careful if you intent to use the
      // normals and do not have the surface tension algorithm enabled.
      // integrateProcessingFunctional (
      //        new FreeSurfaceGeometry2D<T,Descriptor>(contactAngle),
      //        group.getBoundingBox(), freeSurfaceArgs, pl );
    }

    integrateProcessingFunctional(new FreeSurfaceMassChange2D<T, Descriptor>,
                                  group.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(new FreeSurfaceCompletion2D<T, Descriptor>,
                                  group.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceMacroscopic2D<T, Descriptor>(incompressibleModel),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    if (useSurfaceTension) {
      if (surfaceTensionField == 0) {
        integrateProcessingFunctional(
            new FreeSurfaceAddSurfaceTension2D<T, Descriptor>(
                surfaceTension, incompressibleModel),
            group.getBoundingBox(), freeSurfaceArgs, pl);
      } else {
        std::vector<MultiBlock2D *> args(freeSurfaceArgs);
        args.push_back(surfaceTensionField);

        integrateProcessingFunctional(
            new FreeSurfaceAddSurfaceTensionFromScalarField2D<T, Descriptor>(
                incompressibleModel),
            group.getBoundingBox(), args, pl);
      }
    }

    integrateProcessingFunctional(new FreeSurfaceStabilize2D<T, Descriptor>(),
                                  group.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoDefault),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
            dynamics.clone(), force),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoDefault),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    /***** New level ******/
    pl++;

    integrateProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
        group.getBoundingBox(), freeSurfaceArgs, pl);

    integrateProcessingFunctional(
        new FreeSurfaceComputeStatistics2D<T, Descriptor>,
        group.getBoundingBox(), freeSurfaceArgs, pl);

    bool useForce = !util::isZero(norm(force));
    if (useForce) {
      integrateProcessingFunctional(
          new FreeSurfaceAddExternalForce2D<T, Descriptor>(rhoDefault),
          group.getBoundingBox(), freeSurfaceArgs, pl);
    }
  }

  void defaultInitialize(bool useConstRho = true, bool useZeroMomentum = true,
                         bool initializeCell = true) {
    applyProcessingFunctional(new DefaultInitializeFreeSurface2D<T, Descriptor>(
                                  dynamics.clone(), force, rhoDefault,
                                  useConstRho, useZeroMomentum, initializeCell),
                              group.get(fluidname).getBoundingBox(),
                              freeSurfaceArgs);
  }

  void partiallyDefaultInitialize(bool useConstRho = true,
                                  bool useZeroMomentum = true,
                                  bool initializeCell = true) {
    applyProcessingFunctional(
        new PartiallyDefaultInitializeFreeSurface2D<T, Descriptor>(
            dynamics.clone(), force, rhoDefault, useConstRho, useZeroMomentum,
            initializeCell),
        group.get(fluidname).getBoundingBox(), freeSurfaceArgs);
  }

  Group2D &getGroup() { return group; }

  Group2D const &getGroup() const { return group; }

  Dynamics<T, Descriptor> const &getDynamics() const { return dynamics; }

  T getRhoDefault() const { return rhoDefault; }

  T getSurfaceTension() const { return surfaceTension; }

  T getContactAngle() const { return contactAngle; }

  ContactAngleFunction getContactAngleFunction() const {
    return contactAngleFunction;
  }

  Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &getForce() const {
    return force;
  }

  bool getIncompressibleModel() const { return incompressibleModel; }

 private:
  bool fieldExists(std::string name, plint envelopeWidth) {
    if (group.hasBlock(name)) {
      PLB_ASSERT(group.get(name).getMultiBlockManagement().getEnvelopeWidth() >=
                 envelopeWidth);
      return true;
    }
    return false;
  }

 private:
  Group2D &group;
  Dynamics<T, Descriptor> const &dynamics;
  bool incompressibleModel;
  T rhoDefault, surfaceTension, contactAngle;
  ContactAngleFunction contactAngleFunction;
  Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
  std::string fluidname, fsPrefix, rhoBarJprefix;
  bool hasImmersedWalls;

 public:
  std::vector<MultiBlock2D *> freeSurfaceArgs;
  std::vector<MultiBlock2D *> rhoBarJparam;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_MODEL_2D_H
