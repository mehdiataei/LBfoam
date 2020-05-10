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

#ifndef FREE_SURFACE_INITIALIZER_2D_HH
#define FREE_SURFACE_INITIALIZER_2D_HH

#include <cstdlib>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"
#include "lbfoam/models/freeSurfaceInitializer2D.h"

namespace plb {

namespace lbfoam {

/* *************** Class DefaultInitializeFreeSurface2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void DefaultInitializeFreeSurface2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  typedef Descriptor<T> D;

  // In the following, spot the interface cells and tag them.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isEmpty(param.flag(iX, iY))) {
        for (plint iPop = 1; iPop < D::q; iPop++) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];
          // Note: there is no conflict of concurrent read/write on param.flag,
          // because read tests is a cell is fluid, and write only converts
          // empty cells to interface.
          if (isFullWet(param.flag(nextX, nextY))) {
            param.flag(iX, iY) = interface;
          }
        }
      }
    }
  }
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T rho = 0.;
      if (useRhoIni) {
        rho = rhoIni;
      } else {
        rho = param.getDensity(iX, iY);
      }
      Array<T, 2> j((T)0., (T)0.);
      if (!useZeroMomentum) {
        j = param.getMomentum(iX, iY);
      }

      if (initializeCell) {
        T rhoBar = Descriptor<T>::rhoBar(rho);
        T jSqr = normSqr(j);
        T thetaBar = 0.;
        param.cell(iX, iY).getDynamics().computeEquilibria(
            param.cell(iX, iY).getRawPopulations(), rhoBar, j, jSqr, thetaBar);
      }

      param.setDensity(iX, iY, rho);
      param.setMomentum(iX, iY, j);
      param.setForce(iX, iY, force);
      param.outsideDensity(iX, iY) = rhoIni;
      switch (param.flag(iX, iY)) {
        case fluid:
        case protect:
        case temporarilyProtect:
          param.attributeDynamics(iX, iY, dynamicsTemplate->clone());
          param.mass(iX, iY) = rho;
          param.volumeFraction(iX, iY) = (T)1.;
          break;
        case interface:
          param.attributeDynamics(iX, iY, dynamicsTemplate->clone());
          param.mass(iX, iY) = 0.5 * rho;
          param.volumeFraction(iX, iY) = (T)0.5;
          break;
        case empty:
        case protectEmpty:
          param.attributeDynamics(iX, iY,
                                  new NoDynamics<T, Descriptor>(rhoIni));
          param.setForce(iX, iY,
                         Array<T, D::ExternalField::sizeOfForce>::zero());
          param.mass(iX, iY) = (T)0.;
          param.volumeFraction(iX, iY) = (T)0.;
          break;
        case wall:
          param.attributeDynamics(iX, iY,
                                  new BounceBack<T, Descriptor>(rhoIni));
          param.mass(iX, iY) = (T)0.;
          param.volumeFraction(iX, iY) = (T)0.;
          break;
        case slipWall:
          // If no reflection plane is given, SpecularReflection behaves like
          // BounceBack.
          param.attributeDynamics(
              iX, iY, new SpecularReflection<T, Descriptor>(rhoIni));
          param.mass(iX, iY) = (T)0.;
          param.volumeFraction(iX, iY) = (T)0.;
          break;
        default:
          // Invalid free-surface flag.
          PLB_ASSERT(false);
      }
    }
  }
}

/* *************** Class PartiallyDefaultInitializeFreeSurface2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void PartiallyDefaultInitializeFreeSurface2D<T, Descriptor>::
    processGenericBlocks(Box2D domain,
                         std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  typedef Descriptor<T> D;

  // In the following, spot the interface cells and tag them.
  // This time set the volume fraction to 0.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isEmpty(param.flag(iX, iY))) {
        for (plint iPop = 1; iPop < D::q; iPop++) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];
          // Note: there is no conflict of concurrent read/write on param.flag,
          // because read tests is a cell is fluid, and write only converts
          // empty cells to interface.
          if (isFullWet(param.flag(nextX, nextY))) {
            param.flag(iX, iY) = interface;
            param.volumeFraction(iX, iY) = (T)0;
          }
        }
      }
    }
  }
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T rho = 0.;
      if (useRhoIni) {
        rho = rhoIni;
      } else {
        rho = param.getDensity(iX, iY);
      }
      Array<T, 2> j((T)0., (T)0.);
      if (!useZeroMomentum) {
        j = param.getMomentum(iX, iY);
      }

      if (initializeCell) {
        T rhoBar = Descriptor<T>::rhoBar(rho);
        T jSqr = normSqr(j);
        T thetaBar = 0.;
        param.cell(iX, iY).getDynamics().computeEquilibria(
            param.cell(iX, iY).getRawPopulations(), rhoBar, j, jSqr, thetaBar);
      }

      param.setDensity(iX, iY, rho);
      param.setMomentum(iX, iY, j);
      param.setForce(iX, iY, force);
      param.outsideDensity(iX, iY) = rhoIni;
      switch (param.flag(iX, iY)) {
        case fluid:
        case protect:
        case temporarilyProtect:
          param.attributeDynamics(iX, iY, dynamicsTemplate->clone());
          param.mass(iX, iY) = rho;
          break;
        case interface:
          param.attributeDynamics(iX, iY, dynamicsTemplate->clone());
          param.mass(iX, iY) = rho * param.volumeFraction(iX, iY);
          break;
        case empty:
        case protectEmpty:
          param.attributeDynamics(iX, iY,
                                  new NoDynamics<T, Descriptor>(rhoIni));
          param.setForce(iX, iY,
                         Array<T, D::ExternalField::sizeOfForce>::zero());
          param.mass(iX, iY) = (T)0.;
          break;
        case wall:
          param.attributeDynamics(iX, iY,
                                  new BounceBack<T, Descriptor>(rhoIni));
          param.mass(iX, iY) = (T)0.;
          break;
        case slipWall:
          // If no reflection plane is given, SpecularReflection behaves like
          // BounceBack.
          param.attributeDynamics(
              iX, iY, new SpecularReflection<T, Descriptor>(rhoIni));
          param.mass(iX, iY) = (T)0.;
          break;
        default:
          // Invalid free-surface flag.
          PLB_ASSERT(false);
      }
    }
  }
}

/* *************** Class ConstantIniVelocityFreeSurface2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ConstantIniVelocityFreeSurface2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  T rho = rhoIni;
  Array<T, 2> j(velocity * rho);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isWet(param.flag(iX, iY))) {
        iniCellAtEquilibrium(param.cell(iX, iY), rho, velocity);
        param.setMomentum(iX, iY, j);
      }
    }
  }
}

/* *************** Class InletConstVolumeFraction2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void InletConstVolumeFraction2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      param.mass(iX, iY) = param.getDensity(iX, iY) * volumeFraction;
    }
  }
}
/* *************** Class MaskedInletConstVolumeFraction2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void MaskedInletConstVolumeFraction2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 3);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[0]);
  PLB_ASSERT(rhoBar);
  ScalarField2D<T> *mass = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
  PLB_ASSERT(mass);
  ScalarField2D<int> *mask =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[2]);
  PLB_ASSERT(mask);

  Dot2D ofsMass = computeRelativeDisplacement(*rhoBar, *mass);
  Dot2D ofsMask = computeRelativeDisplacement(*rhoBar, *mask);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (mask->get(iX + ofsMask.x, iY + ofsMask.y) == whichFlag) {
        mass->get(iX + ofsMass.x, iY + ofsMass.y) =
            volumeFraction * Descriptor<T>::fullRho(rhoBar->get(iX, iY));
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
void maskedInletConstVolumeFraction2D(MultiScalarField2D<T> &rhoBar,
                                      MultiScalarField2D<T> &mass,
                                      MultiScalarField2D<int> &mask,
                                      T volumeFraction, int whichFlag,
                                      Box2D domain) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&rhoBar);
  args.push_back(&mass);
  args.push_back(&mask);
  applyProcessingFunctional(new MaskedInletConstVolumeFraction2D<T, Descriptor>(
                                volumeFraction, whichFlag),
                            domain, args);
}

/* *************** Class OutletMaximumVolumeFraction2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void OutletMaximumVolumeFraction2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T maximumMass = param.getDensity(iX, iY) * volumeFraction;
      if (param.mass(iX, iY) > maximumMass) {
        param.mass(iX, iY) = maximumMass;
      }
    }
  }
}

/* *************** Class OutletVolumeFractionInRange2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void OutletVolumeFractionInRange2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T maximumMass = param.getDensity(iX, iY) * maxFraction;
      T minimumMass = param.getDensity(iX, iY) * minFraction;
      if (param.mass(iX, iY) > maximumMass) {
        param.mass(iX, iY) = maximumMass;
      } else if (param.mass(iX, iY) < minimumMass) {
        param.mass(iX, iY) = minimumMass;
      }
    }
  }
}

/* *************** Class OutletMaximumVolumeFraction2_2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void OutletMaximumVolumeFraction2_2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T maximumMass = param.getDensity(iX, iY) * volumeFraction;
      if (param.mass(iX, iY) > maximumMass) {
        param.mass(iX, iY) = maximumMass;
      }

      Cell<T, Descriptor> &cell = param.cell(iX, iY);
      T oldRhoBar;
      Array<T, 2> oldJ;
      momentTemplates<T, Descriptor>::get_rhoBar_j(cell, oldRhoBar, oldJ);
      T oldJsqr = normSqr(oldJ);
      T rhoBar = Descriptor<T>::rhoBar(param.getDensity(iX, iY));
      for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T oldEq = cell.getDynamics().computeEquilibrium(iPop, oldRhoBar, oldJ,
                                                        oldJsqr);
        T newEq =
            cell.getDynamics().computeEquilibrium(iPop, rhoBar, oldJ, oldJsqr);
        cell[iPop] += newEq - oldEq;
      }
    }
  }
}

/* *************** Class NoSlipMaximumVolumeFraction2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void NoSlipMaximumVolumeFraction2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T maximumMass = param.getDensity(iX, iY) * volumeFraction;
      if (param.mass(iX, iY) > maximumMass) {
        param.setDensity(iX, iY, param.mass(iX, iY) / volumeFraction);
      }
    }
  }
}

/* *************** Class PunchSphere2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void PunchSphere2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  typedef Descriptor<T> D;

  Dot2D offset = param.absOffset();
  Array<T, 2> localCenter(center - Array<T, 2>(offset.x, offset.y));

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (normSqr(Array<T, 2>(iX, iY) - localCenter) < radius * radius) {
        bool isBoundary = false;
        for (plint iPop = 0; iPop < D::q; ++iPop) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];
          if (normSqr((Array<T, 2>(nextX, nextY) - localCenter)) >=
              radius * radius) {
            isBoundary = true;
          }
        }
        if (isBoundary) {
          param.flag(iX, iY) = interface;
          param.volumeFraction(iX, iY) = (T)0.5;
          param.mass(iX, iY) = 0.5 * rho0;
          param.setDensity(iX, iY, rho0);
          param.outsideDensity(iX, iY) = rho0;
        } else {
          param.flag(iX, iY) = empty;
          param.attributeDynamics(iX, iY, new NoDynamics<T, Descriptor>(rho0));
          param.mass(iX, iY) = T();
          param.volumeFraction(iX, iY) = T();
          param.setDensity(iX, iY, rho0);
          param.setForce(iX, iY,
                         Array<T, D::ExternalField::sizeOfForce>::zero());
          param.setMomentum(iX, iY, Array<T, 2>(T(), T()));
          param.outsideDensity(iX, iY) = rho0;
        }
      }
    }
  }
}

/* *************** Class AnalyticalPunchSphere2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void AnalyticalPunchSphere2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  Dot2D offset = param.absOffset();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    plint globalX = iX + offset.x;
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint globalY = iY + offset.y;

      if (param.flag(iX, iY) == empty || isAnyWall(param.flag(iX, iY))) {
        continue;
      }

      int nextFlag = 0;
      T nextVolumeFraction = 0.;
      subDomainVolumeFraction(globalX, globalY, nextFlag, nextVolumeFraction);

      if (nextFlag == fluid) {
        // Do nothing.
      } else if (nextFlag == empty) {
        param.flag(iX, iY) = empty;
        param.attributeDynamics(iX, iY, new NoDynamics<T, Descriptor>(rho0));
        param.mass(iX, iY) = T();
        param.volumeFraction(iX, iY) = T();
        param.setDensity(iX, iY, rho0);
        // param.setForce(iX,iY,
        // Array<T,Descriptor<T>::ExternalField::sizeOfForce>::zero());
        param.setMomentum(iX, iY, Array<T, 2>(T(), T()));
        param.outsideDensity(iX, iY) = rho0;
      } else {
        param.flag(iX, iY) = interface;
        param.volumeFraction(iX, iY) = nextVolumeFraction;
        param.mass(iX, iY) = nextVolumeFraction * rho0;
        param.setDensity(iX, iY, rho0);
        param.outsideDensity(iX, iY) = rho0;
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
void AnalyticalPunchSphere2D<T, Descriptor>::subDomainVolumeFraction(
    plint globalX, plint globalY, int &flag, T &volumeFraction) {
  plint numInside = 0;
  plint numOutside = 0;

  srand(1.0);
  T xi = (T)globalX - 0.5;
  T yi = (T)globalY - 0.5;
  for (plint xSub = 0; xSub < subDivision; ++xSub) {
    T xPos = xi + (T)rand() / (T)RAND_MAX;
    for (plint ySub = 0; ySub < subDivision; ++ySub) {
      T yPos = yi + (T)rand() / (T)RAND_MAX;
      if (isInsideSphere(xPos, yPos)) {
        ++numInside;
      } else {
        ++numOutside;
      }
    }
  }

  if (numInside == 0) {
    flag = freeSurfaceFlag2D::fluid;
    volumeFraction = (T)1;
  } else if (numOutside == 0) {
    flag = freeSurfaceFlag2D::empty;
    volumeFraction = (T)0;
  } else {
    flag = freeSurfaceFlag2D::interface;
    volumeFraction = (T)numOutside / ((T)numInside + (T)numOutside);
  }
}

/* *************** Class CalculateAverageSphereDensity2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void CalculateAverageSphereDensity2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  Dot2D offset = param.absOffset();
  Array<T, 2> localCenter(center - Array<T, 2>(offset.x, offset.y));
  BlockStatistics &statistics = this->getStatistics();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (normSqr(Array<T, 2>(iX, iY) - localCenter) < radius * radius) {
        statistics.gatherAverage(averageDensityId, param.getDensity(iX, iY));
        statistics.incrementStats();
      }
    }
  }
}

/* *************** Class AnalyticalIniVolumeFraction2D
 * ******************************************* */

template <typename T, class InsideFunction>
void AnalyticalIniVolumeFraction2D<T, InsideFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  PLB_ASSERT(atomicBlocks.size() == 2);
  ScalarField2D<T> *volumeFraction =
      dynamic_cast<ScalarField2D<T> *>(atomicBlocks[0]);
  PLB_ASSERT(volumeFraction);
  ScalarField2D<int> *flag =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[1]);
  PLB_ASSERT(flag);

  Dot2D offset = computeRelativeDisplacement(*volumeFraction, *flag);
  Dot2D absOfs = volumeFraction->getLocation();

  // In the following, spot the interface cells and tag them.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      int nextFlag = 0;
      T nextVolumeFraction = 0.;
      subDomainVolumeFraction(iX + absOfs.x, iY + absOfs.y, nextFlag,
                              nextVolumeFraction);
      volumeFraction->get(iX, iY) = nextVolumeFraction;
      flag->get(iX + offset.x, iY + offset.y) = nextFlag;
    }
  }
}

template <typename T, class InsideFunction>
void AnalyticalIniVolumeFraction2D<T, InsideFunction>::subDomainVolumeFraction(
    plint iX, plint iY, int &flag, T &volumeFraction) {
  plint numInside = 0;
  plint numOutside = 0;

  srand(1.0);
  T xi = (T)iX - 0.5;
  T yi = (T)iY - 0.5;
  for (plint xSub = 0; xSub < subDivision; ++xSub) {
    T xPos = xi + (T)rand() / (T)RAND_MAX;
    for (plint ySub = 0; ySub < subDivision; ++ySub) {
      T yPos = yi + (T)rand() / (T)RAND_MAX;
      if (insideFunction(xPos, yPos)) {
        ++numInside;
      } else {
        ++numOutside;
      }
    }
  }

  if (numInside == 0) {
    flag = freeSurfaceFlag2D::empty;
    volumeFraction = (T)0;
  } else if (numOutside == 0) {
    flag = freeSurfaceFlag2D::fluid;
    volumeFraction = (T)1;
  } else {
    flag = freeSurfaceFlag2D::interface;
    volumeFraction = (T)numInside / ((T)numInside + (T)numOutside);
  }
}

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(MultiScalarField2D<T> &volumeFraction,
                                 MultiScalarField2D<int> &flagStatus,
                                 InsideFunction const &insideFunction,
                                 Box2D domain, plint subDivision) {
  std::vector<MultiBlock2D *> args;
  args.push_back(&volumeFraction);
  args.push_back(&flagStatus);
  applyProcessingFunctional(
      new AnalyticalIniVolumeFraction2D<T, InsideFunction>(insideFunction,
                                                           subDivision),
      domain, args);
}

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(MultiScalarField2D<T> &volumeFraction,
                                 MultiScalarField2D<int> &flagStatus,
                                 InsideFunction const &insideFunction,
                                 plint subDivision) {
  analyticalIniVolumeFraction(volumeFraction, flagStatus, insideFunction,
                              volumeFraction.getBoundingBox(), subDivision);
}

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_INITIALIZER_2D_HH
