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

#ifndef ADVECTION_DIFFUSION_BUBBLE_GROWTH_3D_H
#define ADVECTION_DIFFUSION_BUBBLE_GROWTH_3D_H

#include <memory>

#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "lbfoam/bubble/bubbleGrowth3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {
namespace lbfoam {

/* Mehdi: Coupling between the free-surface model and an advection-diffusion
 * equation for growth. The "empty space" of the free-surface simulation is
 * approximated by an empty medium with empty dynamics.
 */
template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
class GrowthCoupling3D : public BoxProcessingFunctional3D {
 public:
  GrowthCoupling3D(Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics_,
                   Dynamics<T, AD_Descriptor> *emptyDynamics_, T Kh_,
                   std::map<plint, BubbleInfo3D> const &bubbles_,
                   bool surfaceDiffusion_, T source_);
  ~GrowthCoupling3D();
  GrowthCoupling3D(
      GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
  GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> &operator=(
      GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs);
  void swap(GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs);
  virtual void processGenericBlocks(Box3D domain,
                                    std::vector<AtomicBlock3D *> blocks);
  virtual GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> *clone() const;
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const;

 private:
  Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics;
  Dynamics<T, AD_Descriptor> *emptyDynamics;
  T Kh;
  std::map<plint, lbfoam::BubbleInfo3D> bubbles;
  bool surfaceDiffusion;
  T source;
};

}  // namespace lbfoam
}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_H
