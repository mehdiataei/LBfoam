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

#ifndef FREE_SURFACE_BOUNDARY_CONDITION_2D_HH
#define FREE_SURFACE_BOUNDARY_CONDITION_2D_HH

#include <cmath>
#include <iostream>

#include "lbfoam/models/freeSurfaceBoundaryCondition2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
FreeSurfaceFadingArea2D<T, Descriptor>::FreeSurfaceFadingArea2D(T factor_)
    : factor(factor_) {}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceFadingArea2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice) {
  std::vector<T> decomposedVariables;

  enum {
    momentumStoredOffset = Descriptor<T>::ExternalField::momentumBeginsAt,
    densityStoredOffset = Descriptor<T>::ExternalField::densityBeginsAt,
  };

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      Cell<T, Descriptor> &cell = lattice.get(iX, iY);
      plint order = 0;
      cell.getDynamics().decompose(cell, decomposedVariables, order);

      T density = Descriptor<T>::fullRho(decomposedVariables[0]);
      if (density > T(0)) density *= factor;
      decomposedVariables[0] = Descriptor<T>::rhoBar(density);
      cell.getDynamics().recompose(cell, decomposedVariables, order);

      *cell.getExternal(densityStoredOffset) = density;

      Array<T, Descriptor<T>::d> j;
      j.resetToZero();
      T rhoBar;
      momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

      // TODO: What about mass, volumeFraction, flagStatus?
      j.to_cArray(cell.getExternal(momentumStoredOffset));
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceFadingArea2D<T, Descriptor>
    *FreeSurfaceFadingArea2D<T, Descriptor>::clone() const {
  return new FreeSurfaceFadingArea2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void RemoveMass2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      // param.attributeDynamics(iX,iY, new NoDynamics<T,Descriptor>((T)1.));
      param.setDensity(iX, iY, (T)1.);
      param.setMomentum(iX, iY, Array<T, 2>((T)0., (T)0.));
      param.mass(iX, iY) = (T)0;
      param.volumeFraction(iX, iY) = (T)0;
      // param.flag(iX,iY) = empty;
    }
  }
}

template <typename T, template <typename U> class Descriptor>
RemoveMass2D<T, Descriptor> *RemoveMass2D<T, Descriptor>::clone() const {
  return new RemoveMass2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PouringLiquid2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T iniRho = T(1);
      param.attributeDynamics(iX, iY, dynamicsTemplate->clone());
      iniCellAtEquilibrium(param.cell(iX, iY), iniRho, injectionVelocity);
      param.setDensity(iX, iY, iniRho);
      param.setMomentum(iX, iY, iniRho * injectionVelocity);
      param.mass(iX, iY) = iniRho;
      param.volumeFraction(iX, iY) = (T)1;
      param.flag(iX, iY) = fluid;
    }
  }
}

template <typename T, template <typename U> class Descriptor>
PouringLiquid2D<T, Descriptor> *PouringLiquid2D<T, Descriptor>::clone() const {
  return new PouringLiquid2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ShortenBounceBack2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  typedef Descriptor<T> D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  Box2D extDomain = domain.enlarge(1);

  for (plint iX = extDomain.x0; iX <= extDomain.x1; ++iX) {  // mehdi: not sure
    bool xBoundary = iX == extDomain.x0 || iX == extDomain.x1;
    for (plint iY = extDomain.y0; iY <= extDomain.y1; ++iY) {
      if (isWall(param.flag(iX, iY))) {
        bool yBoundary = xBoundary || iY == extDomain.y0 || iY == extDomain.y1;

        for (plint iNeighbor = 1; iNeighbor < D::q; ++iNeighbor) {
          plint nextX = iX + D::c[iNeighbor][0];
          plint nextY = iY + D::c[iNeighbor][1];
          if (!yBoundary || contained(nextX, nextY, domain)) {
            if (isWet(param.flag(nextX, nextY))) {
              plint opp = indexTemplates::opposite<D>(iNeighbor);
              param.cell(nextX, nextY)[iNeighbor] = param.cell(iX, iY)[opp];
            }
          }
        }
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
ShortenBounceBack2D<T, Descriptor> *ShortenBounceBack2D<T, Descriptor>::clone()
    const {
  return new ShortenBounceBack2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone2D<T, Descriptor>::FreeSurfaceSpongeZone2D(
    plint nx_, plint ny_, Array<plint, 4> const &numSpongeCells_,
    Array<T, 4> const &translationParameters_,
    Array<T, 4> const &scaleParameters_, bool incompressibleModel_)
    : nx(nx_),
      ny(ny_),
      numSpongeCells(numSpongeCells_),
      translationParameters(translationParameters_),
      scaleParameters(scaleParameters_),
      incompressibleModel(incompressibleModel_),
      useTanhSpongeFunction(true) {
  for (int iZone = 0; iZone < 4; iZone++) {
    if (numSpongeCells[iZone] > 0) {
      PLB_ASSERT(translationParameters[iZone] > (T)0 &&
                 translationParameters[iZone] < (T)1);
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone2D<T, Descriptor>::FreeSurfaceSpongeZone2D(
    plint nx_, plint ny_, Array<plint, 4> const &numSpongeCells_,
    bool incompressibleModel_)
    : nx(nx_),
      ny(ny_),
      numSpongeCells(numSpongeCells_),
      incompressibleModel(incompressibleModel_),
      useTanhSpongeFunction(false) {}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceSpongeZone2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks) {
  PLB_ASSERT(blocks.size() == 6);

  ScalarField2D<T> *volumeFraction =
      dynamic_cast<ScalarField2D<T> *>(blocks[0]);
  PLB_ASSERT(volumeFraction);
  ScalarField2D<T> *targetVolumeFraction =
      dynamic_cast<ScalarField2D<T> *>(blocks[1]);
  PLB_ASSERT(targetVolumeFraction);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[2]);
  PLB_ASSERT(rhoBar);
  ScalarField2D<T> *mass = dynamic_cast<ScalarField2D<T> *>(blocks[3]);
  PLB_ASSERT(mass);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[4]);
  PLB_ASSERT(j);
  TensorField2D<T, 2> *targetVelocity =
      dynamic_cast<TensorField2D<T, 2> *>(blocks[5]);
  PLB_ASSERT(targetVelocity);

  Dot2D location = volumeFraction->getLocation();

  Dot2D ofsTVF =
      computeRelativeDisplacement(*volumeFraction, *targetVolumeFraction);
  Dot2D ofsRB = computeRelativeDisplacement(*volumeFraction, *rhoBar);
  Dot2D ofsM = computeRelativeDisplacement(*volumeFraction, *mass);
  Dot2D ofsJ = computeRelativeDisplacement(*volumeFraction, *j);
  Dot2D ofsTV = computeRelativeDisplacement(*volumeFraction, *targetVelocity);

  plint spongePositions[4] = {0};
  spongePositions[0] = numSpongeCells[0];
  spongePositions[1] = nx - 1 - numSpongeCells[1];
  spongePositions[2] = numSpongeCells[2];
  spongePositions[3] = ny - 1 - numSpongeCells[3];

  plint spongeCenters[4] = {0};
  T invSigma[4];
  T invNumSpongeCells[4];
  if (useTanhSpongeFunction) {
    spongeCenters[0] =
        util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
    spongeCenters[1] = util::roundToInt(
        spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
    spongeCenters[2] =
        util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
    spongeCenters[3] = util::roundToInt(
        spongePositions[3] + translationParameters[3] * numSpongeCells[3]);

    for (int iZone = 0; iZone < 4; iZone++) {
      invSigma[iZone] =
          (T)1 / (scaleParameters[iZone] * (T)numSpongeCells[iZone]);
    }
  } else {
    for (int iZone = 0; iZone < 4; iZone++) {
      invNumSpongeCells[iZone] = (T)1 / (T)numSpongeCells[iZone];
    }
  }
  T pi = std::acos((T)-1);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    plint x = iX + location.x;
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint y = iY + location.y;

      T spongeFunction = (T)1;
      bool modify = false;
      if (numSpongeCells[0] > 0 && x <= spongePositions[0]) {
        if (useTanhSpongeFunction) {
          spongeFunction *=
              (T)0.5 *
              ((T)1 + std::tanh((T)(x - spongeCenters[0]) * invSigma[0]));
        } else {
          spongeFunction *=
              (T)0.5 * ((T)1 - std::cos(pi * (T)x * invNumSpongeCells[0]));
        }
        modify = true;
      }
      if (numSpongeCells[1] > 0 && x >= spongePositions[1]) {
        if (useTanhSpongeFunction) {
          spongeFunction *=
              (T)0.5 *
              ((T)1 - std::tanh((T)(x - spongeCenters[1]) * invSigma[1]));
        } else {
          spongeFunction *=
              (T)0.5 * ((T)1 + std::cos(pi * (T)(x - spongePositions[1]) *
                                        invNumSpongeCells[1]));
        }
        modify = true;
      }

      if (numSpongeCells[2] > 0 && y <= spongePositions[2]) {
        if (useTanhSpongeFunction) {
          spongeFunction *=
              (T)0.5 *
              ((T)1 + std::tanh((T)(y - spongeCenters[2]) * invSigma[2]));
        } else {
          spongeFunction *=
              (T)0.5 * ((T)1 - std::cos(pi * (T)y * invNumSpongeCells[2]));
        }
        modify = true;
      }
      if (numSpongeCells[3] > 0 && y >= spongePositions[3]) {
        if (useTanhSpongeFunction) {
          spongeFunction *=
              (T)0.5 *
              ((T)1 - std::tanh((T)(y - spongeCenters[3]) * invSigma[3]));
        } else {
          spongeFunction *=
              (T)0.5 * ((T)1 + std::cos(pi * (T)(y - spongePositions[3]) *
                                        invNumSpongeCells[3]));
        }
        modify = true;
      }

      if (modify) {
        T rate = (T)1 - spongeFunction;

        T currentVF = volumeFraction->get(iX, iY);
        T targetVF = targetVolumeFraction->get(iX + ofsTVF.x, iY + ofsTVF.y);
        T newVolumeFraction = currentVF + rate * (targetVF - currentVF);
        volumeFraction->get(iX, iY) = newVolumeFraction;

        T rho = Descriptor<T>::fullRho(rhoBar->get(iX + ofsRB.x, iY + ofsRB.y));
        mass->get(iX + ofsM.x, iY + ofsM.y) = rho * newVolumeFraction;

        Array<T, 2> const &targetV =
            targetVelocity->get(iX + ofsTV.x, iY + ofsTV.y);
        Array<T, 2> const &currentJ = j->get(iX + ofsJ.x, iY + ofsJ.y);
        if (incompressibleModel) {
          j->get(iX + ofsJ.x, iY + ofsJ.y) =
              currentJ + rate * (targetV - currentJ);
        } else {
          j->get(iX + ofsJ.x, iY + ofsJ.y) =
              currentJ + rate * (rho * targetV - currentJ);
        }
      }
    }
  }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone2D<T, Descriptor>
    *FreeSurfaceSpongeZone2D<T, Descriptor>::clone() const {
  return new FreeSurfaceSpongeZone2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceSpongeZone2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::staticVariables;  // Volume fraction.
  modified[1] = modif::nothing;          // Target volume fraction.
  modified[2] = modif::nothing;          // rhoBar.
  modified[3] = modif::staticVariables;  // mass.
  modified[4] = modif::staticVariables;  // j.
  modified[5] = modif::nothing;          // Target j.
}

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_BOUNDARY_CONDITION_2D_HH
