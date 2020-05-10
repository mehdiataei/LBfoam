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

#ifndef FREE_SURFACE_TEMPLATES_2D_HH
#define FREE_SURFACE_TEMPLATES_2D_HH

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/models/freeSurfaceModel2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
struct freeSurfaceTemplates {
  static void massExchangeFluidCell(
      FreeSurfaceProcessorParam2D<T, Descriptor> &param, plint iX, plint iY) {
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag2D;
    // Calculate mass at time t+1 --> eq 6 Thurey's paper.
    for (plint iPop = 0; iPop < D::q; ++iPop) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];
      int nextFlag = param.flag(nextX, nextY);
      if (nextFlag == fluid || nextFlag == interface) {
        // In Thuerey's paper, the mass balance is computed locally on one cell,
        // but N. Thuerey uses outgoing populations. Palabos works with incoming
        // populations and uses the relation f_out_i(x,t) =
        // f_in_opp(i)(x+c_i,t+1).
        plint opp = indexTemplates::opposite<D>(iPop);
        param.mass(iX, iY) +=
            param.cell(iX, iY)[opp] - param.cell(nextX, nextY)[iPop];
      }
    }
  }
};

template <typename T>
struct freeSurfaceTemplates<T, descriptors::ForcedD2Q9Descriptor> {
  static void massExchangeFluidCell(
      FreeSurfaceProcessorParam2D<T, descriptors::ForcedD2Q9Descriptor> &param,
      plint iX, plint iY) {
    using namespace freeSurfaceFlag2D;
    // Calculate mass at time t+1 --> eq 6 Thurey's paper.
    plint nextX, nextY, nextFlag;

    nextX = iX + -1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[6] - param.cell(nextX, nextY)[2];
    }

    nextX = iX + 0;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[8] - param.cell(nextX, nextY)[4];
    }

    nextX = iX + -1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[5] - param.cell(nextX, nextY)[1];
    }
    nextX = iX + -1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[7] - param.cell(nextX, nextY)[3];
    }

    // other side

    nextX = iX + 1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[2] - param.cell(nextX, nextY)[6];
    }

    nextX = iX + 0;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[4] - param.cell(nextX, nextY)[8];
    }

    nextX = iX + 1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[3] - param.cell(nextX, nextY)[7];
    }

    nextX = iX + 1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[1] - param.cell(nextX, nextY)[5];
    }
  }
};

template <typename T>
struct freeSurfaceTemplates<T, descriptors::D2Q9Descriptor> {
  static void massExchangeFluidCell(
      FreeSurfaceProcessorParam2D<T, descriptors::D2Q9Descriptor> &param,
      plint iX, plint iY) {
    using namespace freeSurfaceFlag2D;
    // Calculate mass at time t+1 --> eq 6 Thurey's paper.
    plint nextX, nextY, nextFlag;

    nextX = iX + -1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[6] - param.cell(nextX, nextY)[2];
    }

    nextX = iX + 0;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[8] - param.cell(nextX, nextY)[4];
    }

    nextX = iX + -1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[5] - param.cell(nextX, nextY)[1];
    }
    nextX = iX + -1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[7] - param.cell(nextX, nextY)[3];
    }

    // other side

    nextX = iX + 1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[2] - param.cell(nextX, nextY)[6];
    }

    nextX = iX + 0;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[4] - param.cell(nextX, nextY)[8];
    }

    nextX = iX + 1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[3] - param.cell(nextX, nextY)[7];
    }

    nextX = iX + 1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[1] - param.cell(nextX, nextY)[5];
    }
  }
};

template <typename T>
struct freeSurfaceTemplates<T, descriptors::ExtendedForcedD2Q9Descriptor> {
  static void massExchangeFluidCell(
      FreeSurfaceProcessorParam2D<T, descriptors::ExtendedForcedD2Q9Descriptor>
          &param,
      plint iX, plint iY) {
    using namespace freeSurfaceFlag2D;
    // Calculate mass at time t+1 --> eq 6 Thurey's paper.
    plint nextX, nextY, nextFlag;

    nextX = iX + -1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[6] - param.cell(nextX, nextY)[2];
    }

    nextX = iX + 0;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[8] - param.cell(nextX, nextY)[4];
    }

    nextX = iX + -1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[5] - param.cell(nextX, nextY)[1];
    }
    nextX = iX + -1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[7] - param.cell(nextX, nextY)[3];
    }

    // other side

    nextX = iX + 1;
    nextY = iY + 0;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[2] - param.cell(nextX, nextY)[6];
    }

    nextX = iX + 0;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[4] - param.cell(nextX, nextY)[8];
    }

    nextX = iX + 1;
    nextY = iY + 1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[3] - param.cell(nextX, nextY)[7];
    }

    nextX = iX + 1;
    nextY = iY + -1;
    nextFlag = param.flag(nextX, nextY);
    if (nextFlag == fluid || nextFlag == interface) {
      param.mass(iX, iY) += param.cell(iX, iY)[1] - param.cell(nextX, nextY)[5];
    }
  }
};

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_TEMPLATES_2D_HH
