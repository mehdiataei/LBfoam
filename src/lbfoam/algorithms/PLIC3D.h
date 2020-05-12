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
void cubicmesh(std::vector<std::vector<T>> &vertp) {
  T d0 = 0.;
  T d1 = 1.;

  vertp[0][0] = d1;
  vertp[0][1] = d0;
  vertp[0][2] = d1;
  vertp[1][0] = d1;
  vertp[1][1] = d0;
  vertp[1][2] = d0;
  vertp[2][0] = d1;
  vertp[2][1] = d1;
  vertp[2][2] = d0;
  vertp[3][0] = d1;
  vertp[3][1] = d1;
  vertp[3][2] = d1;
  vertp[4][0] = d0;
  vertp[4][1] = d0;
  vertp[4][2] = d1;
  vertp[5][0] = d0;
  vertp[5][1] = d0;
  vertp[5][2] = d0;
  vertp[6][0] = d0;
  vertp[6][1] = d1;
  vertp[6][2] = d0;
  vertp[7][0] = d0;
  vertp[7][1] = d1;
  vertp[7][2] = d1;
}

template <typename T>
void enforv3dsz(T &c, T &dx, T &dy, T &dz, T &v,
                std::vector<std::vector<T>> &vertp, Array<T, 3> normal) {
  T tol = 1.e-9;
  T cMin = 1.e+14;
  T cMax = -1.0e+14;
  T vt = dx * dy * dz;
  T vBack = v;
  v = v / vt;
  T ci, sn, xm, ym, zm, xmi, ymi, zmi, m1, m2, m3, m12, v1, v2, v3, a0, a1, a2,
      a3, alpha, p0, q0, theta;
  a0 = 0;
  a1 = 0;
  a2 = 0;

  for (int i = 0; i < 8; i++) {
    ci = -1.0 * (vertp[i][0] * normal[0] + vertp[i][1] * normal[1] +
                 vertp[i][2] * normal[2]);

    if (ci <= cMin) {
      cMin = ci;
    }

    if (ci >= cMax) {
      cMax = ci;
    }
  }

  if (((T)vBack / vt) <= 0.5) {
    ci = cMin;
  } else {
    ci = cMax;
    v = 1. - v;
  }

  sn = std::abs(normal[0]) + std::abs(normal[1]) + std::abs(normal[2]);
  xm = normal[0] / sn;
  ym = normal[1] / sn;
  zm = normal[2] / sn;
  xmi = xm * dx;
  ymi = ym * dy;
  zmi = zm * dz;

  sn = std::abs(xmi) + std::abs(ymi) + std::abs(zmi);
  xm = std::abs(xmi) / sn;
  ym = std::abs(ymi) / sn;
  zm = std::abs(zmi) / sn;

  // Defining the binary function
  T myComps[] = {xm, ym, zm};
  std::vector<T> comp(myComps, myComps + sizeof(myComps) / sizeof(T));

  typename std::vector<T>::iterator minE =
      std::min_element(comp.begin(), comp.end());
  typename std::vector<T>::iterator maxE =
      std::max_element(comp.begin(), comp.end());

  m1 = *minE;
  m3 = *maxE;

  if ((m1 == xm && m3 == zm) || (m1 == zm && m3 == xm)) {
    m2 = ym;
  } else {
    if ((m1 == xm && m3 == ym) || (m1 == ym && m3 == xm)) {
      m2 = zm;
    } else {
      m2 = xm;
    }
  }

  m12 = m1 + m2;
  v1 = (std::pow(m1, 2)) / std::max(6 * m2 * m3, tol);
  v2 = v1 + (m2 - m1) / (2. * m3);

  if (m3 < m12) {
    v3 = ((3.0 * m12 - m3) * (std::pow(m3, 2)) +
          (m1 - 3.0 * m3) * (std::pow(m1, 2)) +
          (m2 - 3.0 * m3) * (std::pow(m2, 2))) /
         (6.0 * m1 * m2 * m3);
  } else {
    v3 = m12 / (2 * m3);
  }

  if (v >= v2 && v < v3) {
    a3 = -1.0;
    a2 = 3.0 * m12 / a3;
    a1 = -3.0 * (std::pow(m1, 2) + std::pow(m2, 2)) / a3;
    a0 = (std::pow(m1, 3) + std::pow(m2, 3) - (6.0 * m1 * m2 * m3 * v)) / a3;
    a3 = 1.0;
  } else {
    if (v >= v3 && v <= 0.5 && m3 < m12) {
      a3 = -2.0;
      a2 = 3.0 / a3;
      a1 = -3.0 * (std::pow(m1, 2) + std::pow(m2, 2) + std::pow(m3, 2)) / a3;
      a0 = (std::pow(m1, 3) + std::pow(m2, 3) + std::pow(m3, 3) -
            (6.0 * m1 * m2 * m3 * v)) /
           a3;
      a3 = 1.0;
    }
  }

  if (v >= 0 && v < v1) {
    alpha = std::pow((6 * m1 * m2 * m3 * v), (T)1. / 3.);
  } else if (v >= v1 && v < v2) {
    alpha =
        0.5 * (m1 + std::pow((std::pow(m1, 2) + 8 * m2 * m3 * (v - v1)), 0.5));
  } else if ((v >= v2 && v < v3) || (v >= v3 && v <= 0.5 && m3 < m12)) {
    p0 = (a1 / 3.) - ((std::pow(a2, 2)) / 9.);
    q0 = ((a1 * a2 - 3.0 * a0) / 6.0) - (std::pow(a2, 3)) / 27.0;
    T arg = q0 / (std::pow((-1.0 * std::pow(p0, 3)), 0.5));
    theta = std::acos(arg) / 3.;
    alpha = (std::pow((-1.0 * p0), 0.5)) *
                (std::sin(theta) * (std::pow(3, 0.5)) - std::cos(theta)) -
            (a2 / 3.0);
  } else {
    alpha = m3 * v + m12 / 2.;
  }

  if (((T)vBack / vt) <= 0.5) {
    c = cMin + alpha * std::abs(cMax - cMin);
  } else {
    c = cMax - alpha * std::abs(cMax - cMin);
  }

  v = vBack;
}
}  // namespace lbfoam
}  // namespace plb