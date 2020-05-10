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
#ifndef ADVECTION_DIFFUSION_BUBBLE_GROWTH_3D_HH
#define ADVECTION_DIFFUSION_BUBBLE_GROWTH_3D_HH

#include <cmath>

#include "atomicBlock/blockLattice3D.h"
#include "core/dynamics.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/dynamics/advectionDiffusionBubbleGrowth3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

namespace plb {

namespace lbfoam {

/* Mehdi: Coupling between advection-diffusion for growth and free-surface. */

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::GrowthCoupling3D(
    Dynamics<T, AD_Descriptor> *advectionDiffusionDynamics_,
    Dynamics<T, AD_Descriptor> *emptyDynamics_, T Kh_,
    std::map<plint, lbfoam::BubbleInfo3D> const &bubbles_,
    bool surfaceDiffusion_, T source_)
    : advectionDiffusionDynamics(advectionDiffusionDynamics_),
      emptyDynamics(emptyDynamics_),
      Kh(Kh_),
      bubbles(bubbles_),
      surfaceDiffusion(surfaceDiffusion_),
      source(source_)

{}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::~GrowthCoupling3D() {
  delete advectionDiffusionDynamics;
  delete emptyDynamics;
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks) {
  using namespace freeSurfaceFlag3D;
  typedef typename InterfaceLists3D<T, AD_Descriptor>::Node Node;
  PLB_ASSERT(blocks.size() == 6);

  BlockLattice3D<T, AD_Descriptor> *gasLattice =
      dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(blocks[0]);
  BlockLattice3D<T, FS_Descriptor> *fluidLattice =
      dynamic_cast<BlockLattice3D<T, FS_Descriptor> *>(blocks[1]);
  ScalarField3D<int> *flags = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
  TensorField3D<T, 3> *j = dynamic_cast<TensorField3D<T, 3> *>(blocks[3]);
  ScalarField3D<T> *outsideDensity =
      dynamic_cast<ScalarField3D<T> *>(blocks[4]);
  ScalarField3D<plint> *tags = dynamic_cast<ScalarField3D<plint> *>(blocks[5]);

  plint EmptyId = emptyDynamics->getId();
  plint adID = advectionDiffusionDynamics->getId();

  BlockLattice3D<T, AD_Descriptor> gasBak(*gasLattice);

  PLB_ASSERT(gasLattice);
  PLB_ASSERT(fluidLattice);
  PLB_ASSERT(flags);
  PLB_ASSERT(j);
  PLB_ASSERT(tags);

  std::map<Node, Array<T, AD_Descriptor<T>::q>> neighborOppositePop;
  Dot3D ofsTags = computeRelativeDisplacement(*gasLattice, *tags);
  Dot3D ofsFluid = computeRelativeDisplacement(*gasLattice, *fluidLattice);
  Dot3D ofsFlag = computeRelativeDisplacement(*gasLattice, *flags);
  Dot3D relativeOffsetJ = computeRelativeDisplacement(*gasLattice, *j);
  Dot3D ofsOutsideDensity =
      computeRelativeDisplacement(*gasLattice, *outsideDensity);

  const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
  Array<T, 3> zero;
  zero.resetToZero();

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
        Cell<T, AD_Descriptor> &gasCell = gasLattice->get(iX, iY, iZ);
        Cell<T, FS_Descriptor> &fluidCell = fluidLattice->get(
            iX + ofsFluid.x, iY + ofsFluid.y, iZ + ofsFluid.z);
        plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y, iZ + ofsFlag.z);
        T *u = gasCell.getExternal(velOffset);
        if (isFullWet(flag)) {
          gasLattice->attributeDynamics(iX, iY, iZ,
                                        advectionDiffusionDynamics->clone());
          gasCell.setExternalField(
              AD_Descriptor<T>::ExternalField::scalarBeginsAt,
              AD_Descriptor<T>::ExternalField::sizeOfScalar, &source);
          Array<T, 3> velocity;
          fluidCell.computeVelocity(velocity);
          velocity.to_cArray(u);
        } else {
          gasLattice->attributeDynamics(iX, iY, iZ, emptyDynamics->clone());
          zero.to_cArray(u);
        }

        if (flag == interface) {
          bool needsModification = false;
          Array<T, AD_Descriptor<T>::q> savedPop;
          savedPop[0] = T() - 2.;

          for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop) {
            plint prevX = iX - AD_Descriptor<T>::c[iPop][0];
            plint prevY = iY - AD_Descriptor<T>::c[iPop][1];
            plint prevZ = iZ - AD_Descriptor<T>::c[iPop][2];

            plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);

            bool gasIsEmpty =
                gasBak.get(prevX, prevY, prevZ).getDynamics().getId() ==
                EmptyId;
            // If the f_i[iPop] would be streamed from an empty cell
            if (gasIsEmpty)  // No complete sealing of the melt by the mold.
            {
              savedPop[iPop] = gasLattice->get(
                  iX, iY, iZ)[opp];  // So that g_\bar{i}^in=g(x,t+dt)
              needsModification = true;
            } else {
              savedPop[iPop] = (T)-2.;
            }
          }

          if (surfaceDiffusion == false) {
            plint tag =
                tags->get(iX + ofsTags.x, iY + ofsTags.y, iZ + ofsTags.z);

            std::map<plint, BubbleInfo3D>::const_iterator it =
                bubbles.find(tag);

            if (it->second.isFrozen()) {
              needsModification = false;

              gasLattice->attributeDynamics(
                  iX, iY, iZ, advectionDiffusionDynamics->clone());
              Cell<T, AD_Descriptor> &cell = gasLattice->get(iX, iY, iZ);
              for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                cell[iPop] = T();
              }
              plint numNeighbors = 0;
              plint d = 1;
              for (plint dx = -d; dx <= d; dx++) {
                plint i = iX + dx;
                for (plint dy = -d; dy <= d; dy++) {
                  plint j = iY + dy;
                  for (plint dz = -d; dz <= d; dz++) {
                    plint k = iZ + dz;

                    Cell<T, AD_Descriptor> const &nextCell =
                        gasBak.get(i, j, k);
                    // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y);
                    // if(!(dx==0 && dy==0) &&
                    // (nextCell.getDynamics().getId()==conductiveId) &&
                    //        (nextFlag!=wall)) {
                    if (!(dx == 0 && dy == 0 && dz == 0) &&
                        (nextCell.getDynamics().getId() == adID)) {
                      for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        cell[iPop] += nextCell[iPop];
                      }
                      ++numNeighbors;
                    }
                  }
                }
              }
              for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                if (numNeighbors == 0) {
                  cell[iPop] = AD_Descriptor<T>::t[iPop];
                } else {
                  cell[iPop] /= numNeighbors;
                }
              }
            }
          }

          if (needsModification) {
            neighborOppositePop.insert(
                std::pair<Node, Array<T, AD_Descriptor<T>::q>>(Node(iX, iY, iZ),
                                                               savedPop));
          }

          gasLattice->attributeDynamics(iX, iY, iZ,
                                        advectionDiffusionDynamics->clone());
        }
      }

      typename std::map<Node, Array<T, AD_Descriptor<T>::q>>::const_iterator
          nodes = neighborOppositePop.begin();
      for (; nodes != neighborOppositePop.end(); ++nodes) {
        Node node = nodes->first;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];
        Cell<T, AD_Descriptor> &gasCell = gasLattice->get(iX, iY, iZ);
        Array<T, AD_Descriptor<T>::q> neighborOppPop = nodes->second;
        for (plint iPop = 1; iPop < AD_Descriptor<T>::q; ++iPop) {
          if (neighborOppPop[iPop] > (T)-1.) {
            // Velocity is simply taken from the previous time step.
            Array<T, 3> j_FS =
                j->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y,
                       iZ + relativeOffsetJ.z);
            T jSqr = VectorTemplate<T, FS_Descriptor>::normSqr(j_FS);
            // Remember: the value of pressure on an interface node has been set
            // in FreeSurfaceMacroscopic3D, and is equal to the ambient pressure
            // for a single free-surface fluid, or in the case of a binary
            // pressure, an averaged value.
            T rho_FS_adj = outsideDensity->get(iX + ofsOutsideDensity.x,
                                               iY + ofsOutsideDensity.y,
                                               iZ + ofsOutsideDensity.z) *
                           Kh * AD_Descriptor<T>::cs2;  // Apply Henry's law
            T rhoBar_FS_adj = AD_Descriptor<T>::rhoBar(rho_FS_adj);
            T feq_i =
                gasCell.computeEquilibrium(iPop, rhoBar_FS_adj, j_FS, jSqr);
            plint opp = indexTemplates::opposite<AD_Descriptor<T>>(iPop);
            T feq_opp_i =
                gasCell.computeEquilibrium(opp, rhoBar_FS_adj, j_FS, jSqr);

            // Condition that the opposite iPops are unkown
            if (neighborOppPop[opp] > (T)-1.) {
              gasCell[iPop] = feq_i;
            } else {
              gasCell[iPop] = feq_i + feq_opp_i - neighborOppPop[iPop];
            }
          }
        }
      }
    }
  }
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::GrowthCoupling3D(
    GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs)
    : advectionDiffusionDynamics(rhs.advectionDiffusionDynamics->clone()),
      emptyDynamics(rhs.emptyDynamics->clone()),
      Kh(rhs.Kh),
      bubbles(rhs.bubbles),
      surfaceDiffusion(rhs.surfaceDiffusion)

{}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>
    *GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::clone() const {
  return new GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>(*this);
}

template <typename T, template <typename U1> class AD_Descriptor,
          template <typename U2> class FS_Descriptor>
void GrowthCoupling3D<T, AD_Descriptor, FS_Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const {
  modified[0] = modif::dataStructure;  // Gas
  modified[1] = modif::nothing;        // FreeSurface Fluid
  modified[2] = modif::nothing;        // FreeSurface Flags
  modified[3] = modif::nothing;        // FreeSurface j
  modified[4] = modif::nothing;        // OutsideDensity
  modified[5] = modif::nothing;        // tags
}

}  // namespace lbfoam
}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_HH
