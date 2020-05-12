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

/*
C++ implementaiton of VOFTOOLS v3.1
Copyright (C) 2016 J. Lopez and J. Hernandez released under GNU General Public
License Version 3

Reference:
López, J.; Hernández, J.; Gómez, P.; Faura, F. (2017),
“VOFTools - A software package of calculation tools for volume of
fluid methods using general convex grids”, Mendeley Data, v1"
https://data.mendeley.com/datasets/brrgt645bh/1
*/

#include "core/array.h"

namespace plb {
namespace lbfoam {

template <typename T>
void squareMesh(std::vector<std::vector<T>> &vertp) {
  vertp[0][0] = 0.0;
  vertp[0][1] = 0.0;
  vertp[1][0] = 1.0;
  vertp[1][1] = 0.0;
  vertp[2][0] = 1.0;
  vertp[2][1] = 1.0;
  vertp[3][0] = 0.0;
  vertp[3][1] = 1.0;
}

template <typename T>
void enforv2dsz(T &c, T &dx, T &dy, T &v, std::vector<std::vector<T>> &vertp,
                plb::Array<T, 2> normal) {
  T cMin = 1.0e+14;
  T cMax = -1.0e+14;
  T vt = dx * dy;
  T vBack = v;
  v = v / vt;

  T ci, sn, xm, ym, xmi, ymi, m, m1, v1, alpha;

  for (int i = 0; i < 4; i++) {
    ci = -1. * (vertp[i][0] * normal[0] + vertp[i][1] * normal[1]);

    if (ci <= cMin) {
      cMin = ci;
    }
    if (ci >= cMax) {
      cMax = ci;
    }
  }

  if ((vBack / vt) <= 0.5) {
    ci = cMin;
  } else {
    ci = cMax;
    v = 1. - v;
  }

  sn = std::abs(normal[0]) + std::abs(normal[1]);
  xm = normal[0] / sn;
  ym = normal[1] / sn;
  xmi = xm * dx;
  ymi = ym * dy;
  sn = std::abs(xmi) + std::abs(ymi);
  xm = std::abs(xmi) / sn;
  ym = std::abs(ymi) / sn;

  m1 = std::min(xm, ym);
  m = m1;
  v1 = m / (2. * (1. - m));

  if (v >= 0 && v < v1) {
    alpha = std::sqrt(2. * m * (1. - m) * v);
  } else {
    alpha = v * (1. - m) + m / 2.;
  }

  if ((vBack / vt) <= 0.5) {
    c = cMin + alpha * std::abs(cMax - cMin);
  } else {
    c = cMax - alpha * std::abs(cMax - cMin);
  }

  v = vBack;
}

}  // namespace lbfoam
}  // namespace plb