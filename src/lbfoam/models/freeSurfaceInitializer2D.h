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

#ifndef FREE_SURFACE_INITIALIZER_2D_H
#define FREE_SURFACE_INITIALIZER_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
class DefaultInitializeFreeSurface2D : public BoxProcessingFunctional2D {
 public:
  DefaultInitializeFreeSurface2D(
      Dynamics<T, Descriptor> *dynamicsTemplate_,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
      T rhoIni_ = (T)1., bool useRhoIni_ = true, bool useZeroMomentum_ = true,
      bool initializeCell_ = true)
      : dynamicsTemplate(dynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_) {}
  DefaultInitializeFreeSurface2D(
      DefaultInitializeFreeSurface2D<T, Descriptor> const &rhs)
      : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        force(rhs.force),
        rhoIni(rhs.rhoIni),
        useRhoIni(rhs.useRhoIni),
        useZeroMomentum(rhs.useZeroMomentum),
        initializeCell(rhs.initializeCell) {}
  DefaultInitializeFreeSurface2D<T, Descriptor> *operator=(
      DefaultInitializeFreeSurface2D<T, Descriptor> const &rhs) {
    DefaultInitializeFreeSurface2D<T, Descriptor>(rhs).swap(*this);
    return *this;
  }
  void swap(DefaultInitializeFreeSurface2D<T, Descriptor> &rhs) {
    std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
    std::swap(force, rhs.force);
    std::swap(rhoIni, rhs.rhoIni);
    std::swap(useRhoIni, rhs.useRhoIni);
    std::swap(useZeroMomentum, rhs.useZeroMomentum);
    std::swap(initializeCell, rhs.initializeCell);
  }
  virtual ~DefaultInitializeFreeSurface2D() { delete dynamicsTemplate; }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual DefaultInitializeFreeSurface2D<T, Descriptor> *clone() const {
    return new DefaultInitializeFreeSurface2D(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::dataStructure;    // Fluid
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass
    modified[4] = modif::staticVariables;  // Volume-fraction
    modified[5] = modif::staticVariables;  // Flag-status
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface lists
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::staticVariables;  // Outside density.
  }

 private:
  Dynamics<T, Descriptor> *dynamicsTemplate;
  Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
  T rhoIni;
  bool useRhoIni, useZeroMomentum, initializeCell;
};

// Same as DefaultInitializeFreeSurface, but without initializing the
// Volume-fraction.
template <typename T, template <typename U> class Descriptor>
class PartiallyDefaultInitializeFreeSurface2D
    : public BoxProcessingFunctional2D {
 public:
  PartiallyDefaultInitializeFreeSurface2D(
      Dynamics<T, Descriptor> *dynamicsTemplate_,
      Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
      T rhoIni_ = (T)1., bool useRhoIni_ = true, bool useZeroMomentum_ = true,
      bool initializeCell_ = true)
      : dynamicsTemplate(dynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_) {}
  PartiallyDefaultInitializeFreeSurface2D(
      PartiallyDefaultInitializeFreeSurface2D<T, Descriptor> const &rhs)
      : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        force(rhs.force),
        rhoIni(rhs.rhoIni),
        useRhoIni(rhs.useRhoIni),
        useZeroMomentum(rhs.useZeroMomentum),
        initializeCell(rhs.initializeCell) {}
  PartiallyDefaultInitializeFreeSurface2D<T, Descriptor> *operator=(
      PartiallyDefaultInitializeFreeSurface2D<T, Descriptor> const &rhs) {
    PartiallyDefaultInitializeFreeSurface2D<T, Descriptor>(rhs).swap(*this);
    return *this;
  }
  void swap(PartiallyDefaultInitializeFreeSurface2D<T, Descriptor> &rhs) {
    std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
    std::swap(force, rhs.force);
    std::swap(rhoIni, rhs.rhoIni);
    std::swap(useRhoIni, rhs.useRhoIni);
    std::swap(useZeroMomentum, rhs.useZeroMomentum);
    std::swap(initializeCell, rhs.initializeCell);
  }
  virtual ~PartiallyDefaultInitializeFreeSurface2D() {
    delete dynamicsTemplate;
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual PartiallyDefaultInitializeFreeSurface2D<T, Descriptor> *clone()
      const {
    return new PartiallyDefaultInitializeFreeSurface2D(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::dataStructure;    // Fluid
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass
    modified[4] = modif::staticVariables;  // Volume-fraction
    modified[5] = modif::staticVariables;  // Flag-status
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface lists
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::staticVariables;  // Outside density.
  }

 private:
  Dynamics<T, Descriptor> *dynamicsTemplate;
  Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
  T rhoIni;
  bool useRhoIni, useZeroMomentum, initializeCell;
};

template <typename T, class InsideFunction>
class AnalyticalIniVolumeFraction2D : public BoxProcessingFunctional2D {
 public:
  AnalyticalIniVolumeFraction2D(InsideFunction const &insideFunction_,
                                plint subDivision_ = 5)
      : insideFunction(insideFunction_), subDivision(subDivision_) {
    PLB_ASSERT(subDivision > 1);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual AnalyticalIniVolumeFraction2D<T, InsideFunction> *clone() const {
    return new AnalyticalIniVolumeFraction2D<T, InsideFunction>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;  // Volume-fraction
    modified[1] = modif::staticVariables;  // Flag-status
  }

 private:
  void subDomainVolumeFraction(plint iX, plint iY, int &flag,
                               T &volumeFraction);

 private:
  InsideFunction const &insideFunction;
  plint subDivision;
};

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(MultiScalarField2D<T> &volumeFraction,
                                 MultiScalarField2D<int> &flagStatus,
                                 InsideFunction const &insideFunction,
                                 Box2D domain, plint subDivision = 5);

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(MultiScalarField2D<T> &volumeFraction,
                                 MultiScalarField2D<int> &flagStatus,
                                 InsideFunction const &insideFunction,
                                 plint subDivision = 5);

template <typename T, template <typename U> class Descriptor>
class ConstantIniVelocityFreeSurface2D : public BoxProcessingFunctional2D {
 public:
  ConstantIniVelocityFreeSurface2D(Array<T, 2> velocity_, T rhoIni_)
      : velocity(velocity_), rhoIni(rhoIni_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual ConstantIniVelocityFreeSurface2D<T, Descriptor> *clone() const {
    return new ConstantIniVelocityFreeSurface2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::staticVariables;  // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  Array<T, 2> velocity;
  T rhoIni;
};

template <typename T, template <typename U> class Descriptor>
class InletConstVolumeFraction2D : public BoxProcessingFunctional2D {
 public:
  InletConstVolumeFraction2D(T volumeFraction_)
      : volumeFraction(volumeFraction_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual InletConstVolumeFraction2D<T, Descriptor> *clone() const {
    return new InletConstVolumeFraction2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class MaskedInletConstVolumeFraction2D : public BoxProcessingFunctional2D {
 public:
  MaskedInletConstVolumeFraction2D(T volumeFraction_, int whichFlag_)
      : volumeFraction(volumeFraction_), whichFlag(whichFlag_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual MaskedInletConstVolumeFraction2D<T, Descriptor> *clone() const {
    return new MaskedInletConstVolumeFraction2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::nothing;          // rhoBar.
    modified[1] = modif::staticVariables;  // mass.
    modified[2] = modif::nothing;          // mask.
  }

 private:
  T volumeFraction;
  int whichFlag;
};

template <typename T, template <typename U> class Descriptor>
void maskedInletConstVolumeFraction2D(MultiScalarField2D<T> &rhoBar,
                                      MultiScalarField2D<T> &mass,
                                      MultiScalarField2D<int> &mask,
                                      T volumeFraction, int whichFlag,
                                      Box2D domain);

template <typename T, template <typename U> class Descriptor>
class OutletMaximumVolumeFraction2D : public BoxProcessingFunctional2D {
 public:
  OutletMaximumVolumeFraction2D(T volumeFraction_)
      : volumeFraction(volumeFraction_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual OutletMaximumVolumeFraction2D<T, Descriptor> *clone() const {
    return new OutletMaximumVolumeFraction2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class OutletVolumeFractionInRange2D : public BoxProcessingFunctional2D {
 public:
  OutletVolumeFractionInRange2D(T minFraction_, T maxFraction_)
      : minFraction(minFraction_), maxFraction(maxFraction_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual OutletVolumeFractionInRange2D<T, Descriptor> *clone() const {
    return new OutletVolumeFractionInRange2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T minFraction, maxFraction;
};

template <typename T, template <typename U> class Descriptor>
class OutletMaximumVolumeFraction2_2D : public BoxProcessingFunctional2D {
 public:
  OutletMaximumVolumeFraction2_2D(T volumeFraction_)
      : volumeFraction(volumeFraction_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual OutletMaximumVolumeFraction2_2D<T, Descriptor> *clone() const {
    return new OutletMaximumVolumeFraction2_2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::staticVariables;  // Fluid.
    modified[1] = modif::nothing;          // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class NoSlipMaximumVolumeFraction2D : public BoxProcessingFunctional2D {
 public:
  NoSlipMaximumVolumeFraction2D(T volumeFraction_)
      : volumeFraction(volumeFraction_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual NoSlipMaximumVolumeFraction2D<T, Descriptor> *clone() const {
    return new NoSlipMaximumVolumeFraction2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::nothing;          // j.
    modified[3] = modif::nothing;          // Mass.
    modified[4] = modif::nothing;          // Volume-fraction.
    modified[5] = modif::nothing;          // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::nothing;          // Outside density.
  }

 private:
  T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class PunchSphere2D : public BoxProcessingFunctional2D {
 public:
  PunchSphere2D(Array<T, 2> const &center_, T radius_, T rho0_)
      : center(center_), radius(radius_), rho0(rho0_) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual PunchSphere2D<T, Descriptor> *clone() const {
    return new PunchSphere2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::dataStructure;    // Fluid.
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::staticVariables;  // Volume-fraction.
    modified[5] = modif::staticVariables;  // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::staticVariables;  // Outside density.
  }

 private:
  Array<T, 2> center;
  T radius;
  T rho0;
};

template <typename T, template <typename U> class Descriptor>
class AnalyticalPunchSphere2D : public BoxProcessingFunctional2D {
 public:
  AnalyticalPunchSphere2D(Array<T, 2> const &center_, T radius_, T rho0_,
                          plint subDivision_ = 5)
      : center(center_),
        radius(radius_),
        rho0(rho0_),
        subDivision(subDivision_) {
    PLB_ASSERT(subDivision > 1);
  }
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual AnalyticalPunchSphere2D<T, Descriptor> *clone() const {
    return new AnalyticalPunchSphere2D<T, Descriptor>(*this);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
    modified[0] = modif::dataStructure;    // Fluid.
    modified[1] = modif::staticVariables;  // rhoBar.
    modified[2] = modif::staticVariables;  // j.
    modified[3] = modif::staticVariables;  // Mass.
    modified[4] = modif::staticVariables;  // Volume-fraction.
    modified[5] = modif::staticVariables;  // Flag-status.
    modified[6] = modif::nothing;          // Normal.
    modified[7] = modif::nothing;          // Interface-lists.
    modified[8] = modif::nothing;          // Curvature.
    modified[9] = modif::staticVariables;  // Outside density.
  }

 private:
  bool isInsideSphere(T x, T y) {
    return (normSqr(Array<T, 2>(x, y) - center) < radius * radius);
  }
  void subDomainVolumeFraction(plint globalX, plint globalY, int &flag,
                               T &volumeFraction);

 private:
  Array<T, 2> center;
  T radius;
  T rho0;
  plint subDivision;
};

template <typename T, template <typename U> class Descriptor>
class CalculateAverageSphereDensity2D
    : public PlainReductiveBoxProcessingFunctional2D {
 public:
  CalculateAverageSphereDensity2D(Array<T, 2> const &center_, T radius_)
      : center(center_),
        radius(radius_),
        averageDensityId(this->getStatistics().subscribeAverage()) {}
  virtual void processGenericBlocks(Box2D domain,
                                    std::vector<AtomicBlock2D *> atomicBlocks);
  virtual CalculateAverageSphereDensity2D<T, Descriptor> *clone() const {
    return new CalculateAverageSphereDensity2D<T, Descriptor>(*this);
  }
  T getAverageDensity() const {
    return this->getStatistics().getAverage(averageDensityId);
  }
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    std::fill(modified.begin(), modified.end(), modif::nothing);
  }

 private:
  Array<T, 2> center;
  T radius;
  plint averageDensityId;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_INITIALIZER_2D_H
