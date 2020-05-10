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

#ifndef CREATE_BUBBLES_2D_HH
#define CREATE_BUBBLES_2D_HH

#include <limits>

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "createBubbles2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "offLattice/makeSparse2D.h"

namespace plb {

namespace lbfoam {

template <typename T, template <typename U> class Descriptor>
void punchSphere(FreeSurfaceFields2D<T, Descriptor> &fields,
                 Array<T, 2> const &center, T radius, T rhoEmpty, T rho0,
                 Dynamics<T, Descriptor> &dynamics) {
  applyProcessingFunctional(
      new PunchSphere2D<T, Descriptor>(center, radius, rho0),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoEmpty),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
          dynamics.clone(), Array<T, 2>((T)0., (T)0.)),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoEmpty),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(new FreeSurfaceComputeStatistics2D<T, Descriptor>,
                            fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(FreeSurfaceFields2D<T, Descriptor> &fields,
                           Array<T, 2> const &center, T radius, T rhoEmpty,
                           T rho0, plint subDivision,
                           Dynamics<T, Descriptor> &dynamics) {
  applyProcessingFunctional(new AnalyticalPunchSphere2D<T, Descriptor>(
                                center, radius, rho0, subDivision),
                            fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoEmpty),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
          dynamics.clone(), Array<T, 2>((T)0., (T)0.)),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoEmpty),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
      fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

  applyProcessingFunctional(new FreeSurfaceComputeStatistics2D<T, Descriptor>,
                            fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(FreeSurfaceFields2D<T, Descriptor> &fields,
                              Array<T, 2> const &center, T radius) {
  CalculateAverageSphereDensity2D<T, Descriptor> functional(center, radius);
  applyProcessingFunctional(functional, fields.lattice.getBoundingBox(),
                            fields.freeSurfaceArgs);
  return functional.getAverageDensity();
}

template <typename T, template <typename U> class Descriptor>
void punchSphere(FreeSurfaceSetup2D<T, Descriptor> &setup,
                 Array<T, 2> const &center, T radius, T rhoEmpty, T rho0,
                 Dynamics<T, Descriptor> &dynamics) {
  applyProcessingFunctional(
      new PunchSphere2D<T, Descriptor>(center, radius, rho0),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoEmpty),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
          dynamics.clone(), Array<T, 2>((T)0., (T)0.)),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoEmpty),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(new FreeSurfaceComputeStatistics2D<T, Descriptor>,
                            setup.getGroup().getBoundingBox(),
                            setup.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(FreeSurfaceSetup2D<T, Descriptor> &setup,
                           Array<T, 2> const &center, T radius, T rhoEmpty,
                           T rho0, plint subDivision,
                           Dynamics<T, Descriptor> &dynamics) {
  applyProcessingFunctional(new AnalyticalPunchSphere2D<T, Descriptor>(
                                center, radius, rho0, subDivision),
                            setup.getGroup().getBoundingBox(),
                            setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceComputeInterfaceLists2D<T, Descriptor>(),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>(rhoEmpty),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>(
          dynamics.clone(), Array<T, 2>((T)0., (T)0.)),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>(rhoEmpty),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(
      new FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>(),
      setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

  applyProcessingFunctional(new FreeSurfaceComputeStatistics2D<T, Descriptor>,
                            setup.getGroup().getBoundingBox(),
                            setup.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(FreeSurfaceSetup2D<T, Descriptor> &setup,
                              Array<T, 2> const &center, T radius) {
  CalculateAverageSphereDensity2D<T, Descriptor> functional(center, radius);
  applyProcessingFunctional(functional, setup.getGroup().getBoundingBox(),
                            setup.freeSurfaceArgs);
  return functional.getAverageDensity();
}

}  // namespace lbfoam
}  // namespace plb

#endif  // CREATE_BUBBLES_2D_HH
