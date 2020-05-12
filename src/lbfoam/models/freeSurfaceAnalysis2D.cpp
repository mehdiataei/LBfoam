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

#include "lbfoam/models/freeSurfaceAnalysis2D.h"

#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.hh"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.hh"

namespace plb {
namespace lbfoam {

/* ************** class FreeSurfaceCountWallCells2D
 * ********************************** */

FreeSurfaceCountWallCells2D::FreeSurfaceCountWallCells2D()
    : sumScalarId(this->getStatistics().subscribeIntSum()) {}

void FreeSurfaceCountWallCells2D::process(Box2D domain,
                                          ScalarField2D<int> &flag) {
  BlockStatistics &statistics = this->getStatistics();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (freeSurfaceFlag2D::isAnyWall(flag.get(iX, iY))) {
        statistics.gatherIntSum(sumScalarId, (plint)1);
      }
    }
  }
}

FreeSurfaceCountWallCells2D *FreeSurfaceCountWallCells2D::clone() const {
  return new FreeSurfaceCountWallCells2D(*this);
}

plint FreeSurfaceCountWallCells2D::getNumWallCells() const {
  return this->getStatistics().getIntSum(sumScalarId);
}

plint freeSurfaceCountWallCells(MultiScalarField2D<int> &flag, Box2D domain) {
  FreeSurfaceCountWallCells2D functional;
  applyProcessingFunctional(functional, domain, flag);
  return functional.getNumWallCells();
}

/* ************** class FreeSurfaceComputeFluidBoundingBox2D
 * ********************************** */

FreeSurfaceComputeFluidBoundingBox2D::FreeSurfaceComputeFluidBoundingBox2D()
    : minIdX(this->getStatistics().subscribeMax()),
      minIdY(this->getStatistics().subscribeMax()),
      maxIdX(this->getStatistics().subscribeMax()),
      maxIdY(this->getStatistics().subscribeMax())

{}

void FreeSurfaceComputeFluidBoundingBox2D::process(Box2D domain,
                                                   ScalarField2D<int> &flag) {
  Dot2D absOfs = flag.getLocation();
  BlockStatistics &statistics = this->getStatistics();
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    plint x = iX + absOfs.x;
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint y = iY + absOfs.y;

      if (freeSurfaceFlag2D::isWet(flag.get(iX, iY))) {
        // BlockStatistics computes only maximum, no minimum. Therefore,
        //   the relation min(x) = -max(-x) is used.
        statistics.gatherMax(minIdX, (double)-x);
        statistics.gatherMax(minIdY, (double)-y);
        statistics.gatherMax(maxIdX, (double)x);
        statistics.gatherMax(maxIdY, (double)y);
      }
    }
  }
}

FreeSurfaceComputeFluidBoundingBox2D *
FreeSurfaceComputeFluidBoundingBox2D::clone() const {
  return new FreeSurfaceComputeFluidBoundingBox2D(*this);
}

Box2D FreeSurfaceComputeFluidBoundingBox2D::getFluidBoundingBox() const {
  // The minus sign accounts for the relation min(x) = -max(-x).
  double doubleMinX = -this->getStatistics().getMax(minIdX);
  double doubleMinY = -this->getStatistics().getMax(minIdY);
  double doubleMaxX = this->getStatistics().getMax(maxIdX);
  double doubleMaxY = this->getStatistics().getMax(maxIdY);

  plint minX = util::roundToInt(doubleMinX);
  plint minY = util::roundToInt(doubleMinY);
  plint maxX = util::roundToInt(doubleMaxX);
  plint maxY = util::roundToInt(doubleMaxY);

  return Box2D(minX, maxX, minY, maxY);
}

Box2D freeSurfaceComputeFluidBoundingBox(MultiScalarField2D<int> &flag,
                                         Box2D domain) {
  FreeSurfaceComputeFluidBoundingBox2D functional;
  applyProcessingFunctional(functional, domain, flag);
  return functional.getFluidBoundingBox();
}

Box2D freeSurfaceComputeFluidBoundingBox(MultiScalarField2D<int> &flag) {
  return freeSurfaceComputeFluidBoundingBox(flag, flag.getBoundingBox());
}

}  // namespace lbfoam
}  // namespace plb
