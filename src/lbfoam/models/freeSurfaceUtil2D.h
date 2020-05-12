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

#ifndef FREE_SURFACE_UTIL_2D_H
#define FREE_SURFACE_UTIL_2D_H

#include <set>
#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "finiteDifference/fdWrapper2D.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

namespace lbfoam {

/// Constants used in a free surface flag matrix for cell tagging.
namespace freeSurfaceFlag2D {
enum Flag {
  empty = 0,
  interface = 1,
  fluid = 2,
  wall = 4,
  protect = 5,
  temporarilyProtect = 6,
  protectEmpty = 7,
  slipWall = 8
};
inline std::string flagToString(int flag) {
  switch (flag) {
    case empty:
      return "empty";
    case interface:
      return "interface";
    case fluid:
      return "fluid";
    case wall:
      return "wall";
    case protect:
      return "protect";
    case temporarilyProtect:
      return "temporarilyProtect";
    case protectEmpty:
      return "protectEmpty";
    case slipWall:
      return "slipWall";
    default:
      PLB_ASSERT(false);
  }
  return std::string();
}
inline Flag invert(int flag) {
  switch (flag) {
    case empty:
      return fluid;
    case interface:
      return interface;
    case fluid:
      return empty;
    case wall:
      return wall;
    case protect:
      return protect;
    case temporarilyProtect:
      return temporarilyProtect;
    case protectEmpty:
      return protectEmpty;
    case slipWall:
      return slipWall;
    default:
      PLB_ASSERT(false);
  }
  return (Flag)(-1);
}
inline bool isWet(int flag) {
  return flag == interface || flag == fluid || flag == protect ||
         flag == temporarilyProtect;
}
inline bool isFullWet(int flag) {
  return flag == fluid || flag == protect || flag == temporarilyProtect;
}
inline bool isProtected(int flag) {
  return flag == protect || flag == temporarilyProtect;
}
inline bool isEmpty(int flag) { return flag == empty || flag == protectEmpty; }
inline bool isAnyWall(int flag) { return flag == wall || flag == slipWall; }
inline bool isWall(int flag) { return flag == wall; }
inline bool isSlipWall(int flag) { return flag == slipWall; }
}  // namespace freeSurfaceFlag2D

/// Create a parameter-list for most free-surface data processors.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlock2D *> aggregateFreeSurfaceParams(
    MultiBlockLattice2D<T, Descriptor> &fluid, MultiScalarField2D<T> &rhoBar,
    MultiTensorField2D<T, 2> &j, MultiScalarField2D<T> &mass,
    MultiScalarField2D<T> &volumeFraction, MultiScalarField2D<int> &flag,
    MultiTensorField2D<T, 2> &normal, MultiContainerBlock2D &InterfaceLists2D,
    MultiScalarField2D<T> &curvature, MultiScalarField2D<T> &outsideDensity) {
  std::vector<MultiBlock2D *> aggregation;

  aggregation.push_back(&fluid);
  aggregation.push_back(&rhoBar);
  aggregation.push_back(&j);
  aggregation.push_back(&mass);
  aggregation.push_back(&volumeFraction);
  aggregation.push_back(&flag);
  aggregation.push_back(&normal);
  aggregation.push_back(&InterfaceLists2D);
  aggregation.push_back(&curvature);
  aggregation.push_back(&outsideDensity);

  return aggregation;
}

/// Data structure for holding lists of cells along the free surface in an
/// AtomicContainerBlock.
template <typename T, template <typename U> class Descriptor>
struct InterfaceLists2D : public ContainerBlockData {
  typedef Array<plint, Descriptor<T>::d> Node;
  /// Holds all nodes which have excess mass.
  std::map<Node, T> massExcess;
  /// Holds all nodes that need to change status from interface to fluid.
  std::set<Node> interfaceToFluid;
  /// Holds all nodes that need to change status from interface to empty.
  std::set<Node> interfaceToEmpty;
  /// Holds all nodes that need to change status from empty to interface.
  std::set<Node> emptyToInterface;

  virtual InterfaceLists2D<T, Descriptor> *clone() const {
    return new InterfaceLists2D<T, Descriptor>(*this);
  }
};

/// A wrapper offering convenient access to the free-surface data provided to
/// data processors. Avoids verbous casting, asserting, etc.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceProcessorParam2D {
 public:
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  FreeSurfaceProcessorParam2D(std::vector<AtomicBlock2D *> &atomicBlocks) {
    PLB_ASSERT(atomicBlocks.size() >= 10);

    fluid_ = dynamic_cast<BlockLattice2D<T, Descriptor> *>(atomicBlocks[0]);
    PLB_ASSERT(fluid_);

    rhoBar_ = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
    PLB_ASSERT(rhoBar_);

    j_ = dynamic_cast<TensorField2D<T, 2> *>(atomicBlocks[2]);
    PLB_ASSERT(j_);

    mass_ = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
    PLB_ASSERT(mass_);

    volumeFraction_ = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[4]);
    PLB_ASSERT(volumeFraction_);

    flag_ = dynamic_cast<ScalarField2D<int> *>(atomicBlocks[5]);
    PLB_ASSERT(flag_);

    normal_ = dynamic_cast<TensorField2D<T, 2> *>(atomicBlocks[6]);
    PLB_ASSERT(normal_);

    containerInterfaceLists_ =
        dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[7]);
    PLB_ASSERT(containerInterfaceLists_);

    interfaceLists_ = dynamic_cast<InterfaceLists2D<T, Descriptor> *>(
        containerInterfaceLists_->getData());
    // PLB_ASSERT(interfaceLists_); // This assertion must be commented out to
    // work with TwoPhase.

    curvature_ = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[8]);
    PLB_ASSERT(curvature_);

    outsideDensity_ = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[9]);
    PLB_ASSERT(outsideDensity_);

    absoluteOffset = fluid_->getLocation();
    relativeOffsetRhoBar = computeRelativeDisplacement(*fluid_, *rhoBar_);
    relativeOffsetJ = computeRelativeDisplacement(*fluid_, *j_);
    relativeOffsetMass = computeRelativeDisplacement(*fluid_, *mass_);
    relativeOffsetVF = computeRelativeDisplacement(*fluid_, *volumeFraction_);
    relativeOffsetFS = computeRelativeDisplacement(*fluid_, *flag_);
    relativeOffsetNormal = computeRelativeDisplacement(*fluid_, *normal_);
    relativeOffsetC = computeRelativeDisplacement(*fluid_, *curvature_);
    relativeOffsetOD = computeRelativeDisplacement(*fluid_, *outsideDensity_);
  }

  Cell<T, Descriptor> &cell(plint iX, plint iY) { return fluid_->get(iX, iY); }
  T &mass(plint iX, plint iY) {
    return mass_->get(iX + relativeOffsetMass.x, iY + relativeOffsetMass.y);
  }
  T &volumeFraction(plint iX, plint iY) {
    return volumeFraction_->get(iX + relativeOffsetVF.x,
                                iY + relativeOffsetVF.y);
  }
  T &curvature(plint iX, plint iY) {
    return curvature_->get(iX + relativeOffsetC.x, iY + relativeOffsetC.y);
  }
  T &outsideDensity(plint iX, plint iY) {
    return outsideDensity_->get(iX + relativeOffsetOD.x,
                                iY + relativeOffsetOD.y);
  }
  int &flag(plint iX, plint iY) {
    return flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y);
  }
  void setForce(
      plint iX, plint iY,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force) {
    setExternalForce(cell(iX, iY), force);
  }
  Array<T, Descriptor<T>::ExternalField::sizeOfForce> getForce(plint iX,
                                                               plint iY) {
    return getExternalForce(cell(iX, iY));
  }
  void setMomentum(plint iX, plint iY, Array<T, 2> const &momentum) {
    j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y) = momentum;
  }
  Array<T, 2> getMomentum(plint iX, plint iY) {
    return j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y);
  }
  T getDensity(plint iX, plint iY) {
    return Descriptor<T>::fullRho(
        rhoBar_->get(iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y));
  }
  void setDensity(plint iX, plint iY, T rho) {
    rhoBar_->get(iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y) =
        Descriptor<T>::rhoBar(rho);
  }
  void setNormal(plint iX, plint iY, Array<T, 2> const &normal) {
    normal_->get(iX + relativeOffsetNormal.x, iY + relativeOffsetNormal.y) =
        normal;
  }
  Array<T, 2> getNormal(plint iX, plint iY) {
    return normal_->get(iX + relativeOffsetNormal.x,
                        iY + relativeOffsetNormal.y);
  }

  void attributeDynamics(plint iX, plint iY,
                         Dynamics<T, Descriptor> *dynamics) {
    fluid_->attributeDynamics(iX, iY, dynamics);
  }

  bool isBoundary(plint iX, plint iY) {
    return cell(iX, iY).getDynamics().isBoundary();
  }

  void addToTotalMass(T addedTotalMass) {
    fluid_->getInternalStatistics().gatherSum(0, addedTotalMass);
  }
  void addToLostMass(T addedLostMass) {
    fluid_->getInternalStatistics().gatherSum(1, addedLostMass);
  }
  void addToInterfaceCells(plint addedInterfaceCells) {
    fluid_->getInternalStatistics().gatherIntSum(0, addedInterfaceCells);
  }
  T getSumMassMatrix() const {
    return fluid_->getInternalStatistics().getSum(0);
  }
  T getSumLostMass() const { return fluid_->getInternalStatistics().getSum(1); }
  T getTotalMass() const { return getSumMassMatrix() + getSumLostMass(); }
  plint getNumInterfaceCells() const {
    return fluid_->getInternalStatistics().getIntSum(0);
  }

  T smooth(ScalarField2D<T> const &scalar, Dot2D const &ofs, plint iX,
           plint iY) {
    using namespace freeSurfaceFlag2D;

    if (isAnyWall(
            flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y))) {
      return (scalar.get(iX + ofs.x, iY + ofs.y));
    }

    T val = 0.0;
    int n = 0;
    for (int i = -1; i < 2; i++) {
      plint nextX = iX + i;
      for (int j = -1; j < 2; j++) {
        plint nextY = iY + j;

        if (!(i == 0 && j == 0) &&
            !isAnyWall(flag_->get(nextX + relativeOffsetFS.x,
                                  nextY + relativeOffsetFS.y))) {
          n++;
          val += scalar.get(nextX + ofs.x, nextY + ofs.y);
        }
      }
    }

    if (n != 0) {
      val /= (T)n;
    } else {
      val = scalar.get(iX + ofs.x, iY + ofs.y);
    }

    return (val);
  }

  // The following function can be called with a lattice descriptor
  // "Descriptor2", which is different from the lattice descriptor "Descriptor"
  // used in the simulation.
  template <typename T2, template <typename U2> class Descriptor2>
  T lbmSmooth(ScalarField2D<T> const &scalar, Dot2D const &ofs, plint iX,
              plint iY) {
    typedef Descriptor2<T2> D;
    using namespace freeSurfaceFlag2D;

    if (isAnyWall(
            flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y))) {
      return (scalar.get(iX + ofs.x, iY + ofs.y));
    }

    T val = 0.0, sum = 0.0;
    for (plint iPop = 1; iPop < D::q; iPop++) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];
      sum += D::t[iPop];

      T nextVal = 0.0;

      // First, extrapolate the scalar field on the wall (if necessary).
      if (isAnyWall(flag_->get(nextX + relativeOffsetFS.x,
                               nextY + relativeOffsetFS.y))) {
        T locVal = scalar.get(iX + ofs.x, iY + ofs.y);
        plint opp = indexTemplates::opposite<D>(iPop);
        plint prevX = iX + D::c[opp][0];
        plint prevY = iY + D::c[opp][1];
        if (isAnyWall(flag_->get(prevX + relativeOffsetFS.x,
                                 prevY + relativeOffsetFS.y))) {
          nextVal = locVal;
        } else {
          T prevVal = scalar.get(prevX + ofs.x, prevY + ofs.y);
          nextVal = (T)2 * locVal - prevVal;
        }
      } else {
        nextVal = scalar.get(nextX + ofs.x, nextY + ofs.y);
      }

      val += D::t[iPop] * nextVal;
    }
    val /= sum;

    return (val);
  }

  // The following function can be called with a lattice descriptor
  // "Descriptor2", which is different from the lattice descriptor "Descriptor"
  // used in the simulation.
  template <typename T2, template <typename U2> class Descriptor2>
  Array<T, 2> lbmSmooth(TensorField2D<T, 2> const &tensor, Dot2D const &ofs,
                        plint iX, plint iY) {
    typedef Descriptor2<T2> D;
    using namespace freeSurfaceFlag2D;

    if (isAnyWall(
            flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y))) {
      return (tensor.get(iX + ofs.x, iY + ofs.y));
    }

    T sum = 0.0;
    Array<T, 2> val((T)0, (T)0);
    for (plint iPop = 1; iPop < D::q; iPop++) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];
      sum += D::t[iPop];

      Array<T, 2> nextVal((T)0, (T)0);

      // First, extrapolate the tensor field on the wall (if necessary).
      if (isAnyWall(flag_->get(nextX + relativeOffsetFS.x,
                               nextY + relativeOffsetFS.y))) {
        Array<T, 2> locVal = tensor.get(iX + ofs.x, iY + ofs.y);
        plint opp = indexTemplates::opposite<D>(iPop);
        plint prevX = iX + D::c[opp][0];
        plint prevY = iY + D::c[opp][1];
        if (isAnyWall(flag_->get(prevX + relativeOffsetFS.x,
                                 prevY + relativeOffsetFS.y))) {
          nextVal = locVal;
        } else {
          Array<T, 2> prevVal = tensor.get(prevX + ofs.x, prevY + ofs.y);
          nextVal = (T)2 * locVal - prevVal;
        }
      } else {
        nextVal = tensor.get(nextX + ofs.x, nextY + ofs.y);
      }

      val += D::t[iPop] * nextVal;
    }
    val /= sum;

    return (val);
  }

  // In the free-surface algorithm there is the need to compute derivatives with
  // finite differences. The following function, computes the x/y/z widths and
  // positions for finite difference stencils depending on the local point in
  // the simulation domain. The stencils are tuned so that lattice nodes inside
  // walls are excluded. The maximum stencil width is always an odd number of
  // the form: 2 * h + 1, with h being less or equal to the envelope width of
  // the flag matrix.
  void getFdStencilWidthsAndPositions(plint iX, plint iY, plint h,
                                      Array<int, 2> &widths,
                                      Array<int, 2> &positions) {
    using namespace freeSurfaceFlag2D;

    plint x = iX + relativeOffsetFS.x;
    plint y = iY + relativeOffsetFS.y;

    // PLB_ASSERT(!isAnyWall(flag_->get(x, y, z)));

    int left, right;

    // x-direction.

    left = 0;
    for (plint d = 1; d <= h; d++, left++) {
      if (isAnyWall(flag_->get(x - d, y))) {
        break;
      }
    }
    right = 0;
    for (plint d = 1; d <= h; d++, right++) {
      if (isAnyWall(flag_->get(x + d, y))) {
        break;
      }
    }

    widths[0] = left + right + 1;
    positions[0] = left;

    // y-direction.

    left = 0;
    for (plint d = 1; d <= h; d++, left++) {
      if (isAnyWall(flag_->get(x, y - d))) {
        break;
      }
    }
    right = 0;
    for (plint d = 1; d <= h; d++, right++) {
      if (isAnyWall(flag_->get(x, y + d))) {
        break;
      }
    }

    widths[1] = left + right + 1;
    positions[1] = left;
  }

  // Compute the gradient of a scalar field with finite differences excluding
  // the wall cells.
  Array<T, 2> computeGradient(ScalarField2D<T> const &scalar, Dot2D const &ofs,
                              plint h, plint iX, plint iY) {
    using namespace freeSurfaceFlag2D;

    // PLB_ASSERT(h >= 0 && h <= 3);
    // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY +
    // relativeOffsetFS.y)));

    Array<int, 2> widths, positions;
    getFdStencilWidthsAndPositions(iX, iY, h, widths, positions);

    plint i = iX + ofs.x;
    plint j = iY + ofs.y;
    Array<T, 2> gradient;

    // Template parameters:
    // 1: because we compute first order derivatives
    // 7: this is at least (2 * h + 1), with "h" being at most equal to the
    // envelope width of the flag matrix
    //    We set this to its current possible maximum value 7, to be more
    //    efficient in case we want to use high-order finite differences. This
    //    must be changed if "h" is ever expected to be larger than 3.
    gradient[0] = computeScalarXderivative<T, 1, 7>(scalar, widths[0],
                                                    positions[0], i, j);
    gradient[1] = computeScalarYderivative<T, 1, 7>(scalar, widths[1],
                                                    positions[1], i, j);

    return (gradient);
  }

  // Compute the divergence of a tensor field with finite differences excluding
  // the wall cells.
  T computeDivergence(TensorField2D<T, 2> const &tensor, Dot2D const &ofs,
                      plint h, plint iX, plint iY) {
    using namespace freeSurfaceFlag2D;

    // PLB_ASSERT(h >= 0 && h <= 3);
    // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY +
    // relativeOffsetFS.y)));

    Array<int, 2> widths, positions;
    getFdStencilWidthsAndPositions(iX, iY, h, widths, positions);

    plint i = iX + ofs.x;
    plint j = iY + ofs.y;

    // Template parameters:
    // 3: because "tensor" contains 2D vectors
    // 1: because we compute first order derivatives
    // 7: this is at least (2 * h + 1), with "h" being at most equal to the
    // envelope width of the flag matrix
    //    We set this to its current possible maximum value 7, to be more
    //    efficient in case we want to use high-order finite differences. This
    //    must be changed if "h" is ever expected to be larger than 3.
    Array<T, 2> d_dx = computeTensorXderivative<T, 2, 1, 7>(tensor, widths[0],
                                                            positions[0], i, j);
    Array<T, 2> d_dy = computeTensorYderivative<T, 2, 1, 7>(tensor, widths[1],
                                                            positions[1], i, j);

    T divergence = d_dx[0] + d_dy[1];

    return (divergence);
  }

  // Compute the gradient of a scalar field "the lattice Boltzmann way".
  // With this method the scalar field values are first extrapolated on the
  // neighboring wall cells (if necessary). So, in contrast to the above finite
  // difference implementation, the wall cells are not excluded here. The
  // following function can be called with a lattice descriptor "Descriptor2",
  // which is different from the lattice descriptor "Descriptor" used in the
  // simulation.
  template <typename T2, template <typename U2> class Descriptor2>
  Array<T, 2> lbmComputeGradient(ScalarField2D<T> const &scalar,
                                 Dot2D const &ofs, plint iX, plint iY) {
    typedef Descriptor2<T2> D;
    using namespace freeSurfaceFlag2D;

    // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY +
    // relativeOffsetFS.y)));

    Array<T, 2> gradient((T)0, (T)0);
    for (plint iPop = 1; iPop < D::q; ++iPop) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];

      T nextVal = 0.0;

      // First, extrapolate the scalar field on the wall (if necessary).
      if (isAnyWall(flag_->get(nextX + relativeOffsetFS.x,
                               nextY + relativeOffsetFS.y))) {
        T locVal = scalar.get(iX + ofs.x, iY + ofs.y);
        plint opp = indexTemplates::opposite<D>(iPop);
        plint prevX = iX + D::c[opp][0];
        plint prevY = iY + D::c[opp][1];
        if (isAnyWall(flag_->get(prevX + relativeOffsetFS.x,
                                 prevY + relativeOffsetFS.y))) {
          nextVal = locVal;
        } else {
          T prevVal = scalar.get(prevX + ofs.x, prevY + ofs.y);
          nextVal = (T)2 * locVal - prevVal;
        }
      } else {
        nextVal = scalar.get(nextX + ofs.x, nextY + ofs.y);
      }

      gradient[0] += D::t[iPop] * D::c[iPop][0] * nextVal;
      gradient[1] += D::t[iPop] * D::c[iPop][1] * nextVal;
    }
    gradient *= D::invCs2;

    return (gradient);
  }

  // Compute the divergence of a tensor field "the lattice Boltzmann way".
  // With this method the tensor field values are first extrapolated on the
  // neighboring wall cells (if necessary). So, in contrast to the above finite
  // difference implementation, the wall cells are not excluded here. The
  // following function can be called with a lattice descriptor "Descriptor2",
  // which is different from the lattice descriptor "Descriptor" used in the
  // simulation.
  template <typename T2, template <typename U2> class Descriptor2>
  T lbmComputeDivergence(TensorField2D<T, 2> const &tensor, Dot2D const &ofs,
                         plint iX, plint iY) {
    typedef Descriptor2<T2> D;
    using namespace freeSurfaceFlag2D;

    // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY +
    // relativeOffsetFS.y)));

    T divergence = 0.0;
    for (plint iPop = 1; iPop < D::q; ++iPop) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];

      Array<T, 2> nextVal((T)0, (T)0);

      // First, extrapolate the scalar field on the wall (if necessary).
      if (isAnyWall(flag_->get(nextX + relativeOffsetFS.x,
                               nextY + relativeOffsetFS.y))) {
        Array<T, 2> locVal = tensor.get(iX + ofs.x, iY + ofs.y);
        plint opp = indexTemplates::opposite<D>(iPop);
        plint prevX = iX + D::c[opp][0];
        plint prevY = iY + D::c[opp][1];
        if (isAnyWall(flag_->get(prevX + relativeOffsetFS.x,
                                 prevY + relativeOffsetFS.y))) {
          nextVal = locVal;
        } else {
          Array<T, 2> prevVal = tensor.get(prevX + ofs.x, prevY + ofs.y);
          nextVal = (T)2 * locVal - prevVal;
        }
      } else {
        nextVal = tensor.get(nextX + ofs.x, nextY + ofs.y);
      }

      divergence += D::t[iPop] *
                    (D::c[iPop][0] * nextVal[0] + D::c[iPop][1] * nextVal[1]);
    }
    divergence *= D::invCs2;

    return (divergence);
  }

  std::map<Node, T> &massExcess() { return interfaceLists_->massExcess; }
  std::set<Node> &interfaceToFluid() {
    return interfaceLists_->interfaceToFluid;
  }
  std::set<Node> &interfaceToEmpty() {
    return interfaceLists_->interfaceToEmpty;
  }
  std::set<Node> &emptyToInterface() {
    return interfaceLists_->emptyToInterface;
  }

  BlockLattice2D<T, Descriptor> *latticeP() { return fluid_; }
  ScalarField2D<T> *rhoBarP() { return rhoBar_; }
  TensorField2D<T, 2> *jP() { return j_; }
  ScalarField2D<T> *massP() { return mass_; }
  ScalarField2D<T> *volumeFractionP() { return volumeFraction_; }
  ScalarField2D<int> *flagP() { return flag_; }
  TensorField2D<T, 2> *normalP() { return normal_; }
  ScalarField2D<T> *curvatureP() { return curvature_; }
  ScalarField2D<T> *outsideDensityP() { return outsideDensity_; }

  Dot2D const &absOffset() const { return absoluteOffset; }
  Dot2D const &rhoBarOffset() const { return relativeOffsetRhoBar; }
  Dot2D const &jOffset() const { return relativeOffsetJ; }
  Dot2D const &massOffset() const { return relativeOffsetMass; }
  Dot2D const &volumeFractionOffset() const { return relativeOffsetVF; }
  Dot2D const &flagOffset() const { return relativeOffsetFS; }
  Dot2D const &normalOffset() const { return relativeOffsetNormal; }
  Dot2D const &curvatureOffset() const { return relativeOffsetC; }
  Dot2D const &outsideDensityOffset() const { return relativeOffsetOD; }

  Box2D getBoundingBox() const { return volumeFraction_->getBoundingBox(); }

 private:
  BlockLattice2D<T, Descriptor> *fluid_;
  ScalarField2D<T> *rhoBar_;
  TensorField2D<T, 2> *j_;
  ScalarField2D<T> *mass_;
  ScalarField2D<T> *volumeFraction_;
  ScalarField2D<int> *flag_;
  TensorField2D<T, 2> *normal_;
  AtomicContainerBlock2D *containerInterfaceLists_;
  InterfaceLists2D<T, Descriptor> *interfaceLists_;
  ScalarField2D<T> *curvature_;
  ScalarField2D<T> *outsideDensity_;

  Dot2D absoluteOffset, relativeOffsetRhoBar, relativeOffsetJ,
      relativeOffsetMass, relativeOffsetVF, relativeOffsetFS,
      relativeOffsetNormal, relativeOffsetC, relativeOffsetOD;
};

template <typename T>
std::auto_ptr<MultiScalarField2D<T>> computeNormalizedVolumeFraction(
    MultiScalarField2D<T> &volumeFraction, Box2D domain) {
  std::auto_ptr<MultiScalarField2D<T>> result =
      extractSubDomain<T>(volumeFraction, domain);
  boundScalarField(*result, domain, (T)0, (T)1);
  return result;
}

template <typename T>
std::auto_ptr<MultiScalarField2D<T>> computeNormalizedVolumeFraction(
    MultiScalarField2D<T> &volumeFraction) {
  return computeNormalizedVolumeFraction(volumeFraction,
                                         volumeFraction.getBoundingBox());
}

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_UTIL_2D_H
