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

#ifndef FREE_SURFACE_MODEL_2D_HH
#define FREE_SURFACE_MODEL_2D_HH

#include <limits>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "lbfoam/models/freeSurfaceModel2D.h"
#include "lbfoam/models/freeSurfaceTemplates2D.h"

namespace plb {

namespace lbfoam {

/* *************** Class FreeSurfaceComputeNormals2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeNormals2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // Smooth the volume fraction twice. (At the end include also a 1-cell layer
  // around "domain".)
  plint nx = domain.getNx() + 4;
  plint ny = domain.getNy() + 4;
  ScalarField2D<T> smoothVolumeFractionTmp(nx, ny);
  Dot2D ofsSVFT(-domain.x0 + 2, -domain.y0 + 2);
  for (plint iX = domain.x0 - 2; iX <= domain.x1 + 2; ++iX) {
    plint i = iX + ofsSVFT.x;
    for (plint iY = domain.y0 - 2; iY <= domain.y1 + 2; ++iY) {
      plint j = iY + ofsSVFT.y;

      // smoothVolumeFractionTmp.get(i, j) =
      // param.smooth(*param.volumeFractionP(),
      //        param.volumeFractionOffset(), iX, iY);
      // smoothVolumeFractionTmp.get(i, j) =
      //    param.template lbmSmooth<T, Descriptor>(
      //            *param.volumeFractionP(), param.volumeFractionOffset(), iX,
      //            iY);
      smoothVolumeFractionTmp.get(i, j) = param.template lbmSmooth<
          T, descriptors::AdvectionDiffusionD2Q5Descriptor>(
          *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY);
    }
  }

  nx = domain.getNx() + 2;
  ny = domain.getNy() + 2;
  ScalarField2D<T> smoothVolumeFraction(nx, ny);
  Dot2D ofsSVF(-domain.x0 + 1, -domain.y0 + 1);
  for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
    plint i = iX + ofsSVF.x;
    for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
      plint j = iY + ofsSVF.y;

      // smoothVolumeFraction.get(i, j) = param.smooth(smoothVolumeFractionTmp,
      // ofsSVFT, iX, iY); smoothVolumeFraction.get(i, j) =
      //    param.template lbmSmooth<T, Descriptor>(
      //            smoothVolumeFractionTmp, ofsSVFT, iX, iY);
      smoothVolumeFraction.get(i, j) = param.template lbmSmooth<
          T, descriptors::AdvectionDiffusionD2Q5Descriptor>(
          smoothVolumeFractionTmp, ofsSVFT, iX, iY);
    }
  }

  // The outward pointing unit normal is: n = - grad(VOF) / ||grad(VOF)||.

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      Array<T, 2> normal((T)0, (T)0);

      if (isAnyWall(param.flag(iX, iY))) {
        param.setNormal(iX, iY, normal);
        continue;
      }

      /*
          int useLB = 1;
          for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
              plint nextX = iX + Descriptor<T>::c[iPop][0];
              plint nextY = iY + Descriptor<T>::c[iPop][1];
              if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                  useLB = 0;
                  break;
              }
          }

          if (useLB) {
          */
      // Compute the gradient of the smoothed volume fraction "the lattice
      // Boltzmann way". normal = param.template lbmComputeGradient<T,
      // Descriptor>(smoothVolumeFraction, ofsSVF, iX, iY);

      //     } else {
      // Compute the gradient of the smoothed volume fraction with finite
      // differences excluding the wall cells.
      plint h = 1;  // This must be 1 because above we included only a 1-cell
                    // layer around "domain".
      normal = param.computeGradient(smoothVolumeFraction, ofsSVF, h, iX, iY);
      //   }

      T nn = norm(normal);
      if (util::isZero(nn)) {
        normal.resetToZero();
      } else {
        normal /= -nn;
      }
      param.setNormal(iX, iY, normal);
    }
  }
}

/* *************** Class FreeSurfaceGeometry2D ********************************
 */

template <typename T, template <typename U> class Descriptor>
ScalarField2D<int> *FreeSurfaceGeometry2D<T, Descriptor>::getInterfaceFlags(
    Box2D domain, FreeSurfaceProcessorParam2D<T, Descriptor> &param) {
  using namespace freeSurfaceFlag2D;

  // Define a temporary scalar field for local use in this function. This scalar
  // field will contain 1 extra layer of cells around "domain".
  plint nx = domain.x1 - domain.x0 + 1;
  plint ny = domain.y1 - domain.y0 + 1;
  ScalarField2D<int> *tmp =
      new ScalarField2D<int>(nx + 2, ny + 2, (int)unTagged);
  PLB_ASSERT(tmp);

  // First tag all regular and contact line interface cells. (Loop along 1
  // envelope cell as well). All interface tags are stored in the temporary
  // storage.
  for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
    plint indX = iX - domain.x0 + 1;
    for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
      plint indY = iY - domain.y0 + 1;

      if (param.flag(iX, iY) != interface) {
        tmp->get(indX, indY) = notInterface;
        continue;
      }

      // Find all wall neighbors and store their indices.
      int numWallNeighbors = 0;
      std::vector<Array<plint, 2>> wallNeighborIndex;
      for (int dx = -1; dx < 2; dx++) {
        plint i = iX + dx;
        for (int dy = -1; dy < 2; dy++) {
          plint j = iY + dy;

          if (!(dx == 0 && dy == 0)) {
            if (isAnyWall(param.flag(i, j))) {
              numWallNeighbors++;
              wallNeighborIndex.push_back(Array<plint, 2>(i, j));
            }
          }
        }
      }

      if (numWallNeighbors == 0) {
        tmp->get(indX, indY) = regular;
        continue;
      }

      for (int dx = -1; dx < 2; dx++) {
        plint i = iX + dx;
        for (int dy = -1; dy < 2; dy++) {
          plint j = iY + dy;

          if (!(dx == 0 && dy == 0)) {
            if ((contactAngle > 90.0 && isFullWet(param.flag(i, j))) ||
                (contactAngle <= 90.0 && isEmpty(param.flag(i, j)))) {
              for (int dxx = -1; dxx < 2; dxx++) {
                plint ii = i + dxx;
                for (int dyy = -1; dyy < 2; dyy++) {
                  plint jj = j + dyy;

                  if (!(dxx == 0 && dyy == 0)) {
                    if (isAnyWall(param.flag(ii, jj))) {
                      for (int iWall = 0; iWall < numWallNeighbors; iWall++) {
                        if (ii == wallNeighborIndex[iWall][0] &&
                            jj == wallNeighborIndex[iWall][1]) {
                          tmp->get(indX, indY) = contactLine;
                          goto label0;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    label0:
      continue;
    }
  }

  // Define a scalar field with the interface flags that will be returned from
  // this function.
  ScalarField2D<int> *interfaceFlag =
      new ScalarField2D<int>(nx, ny, (int)unTagged);
  PLB_ASSERT(interfaceFlag);

  // Now tag all adjacent interface cells and copy all information to the scalar
  // field to be returned. At this point all cells that have the flag "unTagged"
  // are non-contact-line interface cells with wall neighbors, so they are
  // either "regular" or "adjacent" interface cells.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    plint indXtmp = iX - domain.x0 + 1;
    plint indX = iX - domain.x0;
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      plint indYtmp = iY - domain.y0 + 1;
      plint indY = iY - domain.y0;

      if (tmp->get(indXtmp, indYtmp) != unTagged) {
        interfaceFlag->get(indX, indY) = tmp->get(indXtmp, indYtmp);
      } else {
        int isAdjacent = 0;
        for (int dx = -1; dx < 2; dx++) {
          plint i = indXtmp + dx;
          for (int dy = -1; dy < 2; dy++) {
            plint j = indYtmp + dy;

            if (!(dx == 0 && dy == 0)) {
              if (tmp->get(i, j) == contactLine) {
                isAdjacent = 1;
                interfaceFlag->get(indX, indY) = adjacent;
                goto label1;
              }
            }
          }
        }
      label1:
        if (!isAdjacent) {
          interfaceFlag->get(indX, indY) = regular;
        }
      }
    }
  }

  // Check for untagged cells
#ifdef PLB_DEBUG
  for (plint i = 0; i < nx; i++) {
    for (plint j = 0; j < ny; j++) {
      PLB_ASSERT(interfaceFlag->get(i, j) != unTagged);
    }
  }
#endif

  delete tmp;

  return interfaceFlag;
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceGeometry2D<T, Descriptor>::computeHeights2D(
    FreeSurfaceProcessorParam2D<T, Descriptor> &param, int integrationDirection,
    plint iX, plint iY, T h[3]) {
  using namespace freeSurfaceFlag2D;

  // Compute the vector parallel to the integration direction.
  Array<int, 2> integrationVector;
  integrationVector[0] = integrationDirection == 0 ? 1 : 0;
  integrationVector[1] = integrationDirection == 1 ? 1 : 0;

  // Compute the vectors tangent to the plane which is normal to the integration
  // vector.
  int iTangentDirection0 = integrationDirection == 0 ? 1 : 0;
  Array<int, 2> tangent0;
  tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
  tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;

  // Calculate the integration stencil width.
  int maxLim = 3;
  for (int d0 = -1; d0 <= 1; d0++) {
    plint posX = iX + d0 * tangent0[0];
    plint posY = iY + d0 * tangent0[1];
    if (isAnyWall(param.flag(posX, posY))) {
      continue;
    }
    for (int d = 1; d <= maxLim; d++) {
      plint nextX = posX + d * integrationVector[0];
      plint nextY = posY + d * integrationVector[1];
      if (isAnyWall(param.flag(nextX, nextY))) {
        maxLim = std::min(maxLim, d - 1);
        break;
      }
    }
  }

  int minLim = 3;
  for (int d0 = -1; d0 <= 1; d0++) {
    plint posX = iX + d0 * tangent0[0];
    plint posY = iY + d0 * tangent0[1];
    if (isAnyWall(param.flag(posX, posY))) {
      continue;
    }
    for (int d = 1; d <= minLim; d++) {
      plint nextX = posX - d * integrationVector[0];
      plint nextY = posY - d * integrationVector[1];
      if (isAnyWall(param.flag(nextX, nextY))) {
        minLim = std::min(minLim, d - 1);
        break;
      }
    }
  }

  // Properly initialize heights to -1.
  for (int i = 0; i < 3; i++) {
    h[i] = -1.0;
  }

  // Integrate.
  for (int d0 = -1; d0 <= 1; d0++) {
    int i = d0 + 1;

    plint posX = iX + d0 * tangent0[0];
    plint posY = iY + d0 * tangent0[1];
    if (isAnyWall(param.flag(posX, posY))) {
      continue;
    }
    h[i] = 0.0;
    for (int d = -minLim; d <= maxLim; d++) {
      plint nextX = posX + d * integrationVector[0];
      plint nextY = posY + d * integrationVector[1];
      h[i] += param.volumeFraction(nextX, nextY);
    }
  }

  // Extrapolate on walls. (No contact angle algorithm).
  for (int i = 0; i < 3; i++) {
    if (util::isZero(h[i] + (T)1)) {
      h[i] = h[1];
    }
  }
}

// template<typename T,template<typename U> class Descriptor> // disabled for
// contact angle mehdi void
// FreeSurfaceGeometry2D<T,Descriptor>::computeHeights2DC(FreeSurfaceProcessorParam2D<T,Descriptor>&
// param,
//        Array<int,2>& wallTangent0, int integrationDirection, plint iX, plint
//        iY, T h[3])
//{
//    using namespace freeSurfaceFlag2D;

//    // Compute the vector parallel to the integration direction.
//    Array<int,2> integrationVector;
//    integrationVector[0] = integrationDirection == 0 ? wallTangent0[0] :
//    wallTangent1[0]; //mehdi: not sure integrationVector[1] =
//    integrationDirection == 0 ? wallTangent0[1] : wallTangent1[1];
//    // Compute the vector tangent to the line which is normal to the
//    integration vector. Array<int,2> tangent; tangent[0] =
//    integrationDirection == 0 ? wallTangent1[0] : wallTangent0[0]; tangent[1]
//    = integrationDirection == 0 ? wallTangent1[1] : wallTangent0[1];

//    // Calculate the integration stencil width.
//    int maxLim = 3;
//    for (int d0 = -1; d0 <= 1; d0++) {
//        plint posX = iX + d0 * tangent[0];
//        plint posY = iY + d0 * tangent[1];
//        if (isAnyWall(param.flag(posX, posY))) {
//            continue;
//        }
//        for (int d = 1; d <= maxLim; d++) {
//            plint nextX = posX + d * integrationVector[0];
//            plint nextY = posY + d * integrationVector[1];
//            if (isAnyWall(param.flag(nextX, nextY))) {
//                maxLim = std::min(maxLim, d - 1);
//                break;
//            }
//        }
//    }

//    int minLim = 3;
//    for (int d0 = -1; d0 <= 1; d0++) {
//        plint posX = iX + d0 * tangent[0];
//        plint posY = iY + d0 * tangent[1];
//        if (isAnyWall(param.flag(posX, posY))) {
//            continue;
//        }
//        for (int d = 1; d <= minLim; d++) {
//            plint nextX = posX - d * integrationVector[0];
//            plint nextY = posY - d * integrationVector[1];
//            if (isAnyWall(param.flag(nextX, nextY))) {
//                minLim = std::min(minLim, d - 1);
//                break;
//            }
//        }
//    }

//    // Properly initialize heights to -1.
//    h[0] = h[1] = h[2] = -1.0;

//    // Integrate.
//    for (int d0 = -1; d0 <= 1; d0++) {
//        int i = d0 + 1;
//        plint posX = iX + d0 * tangent[0];
//        plint posY = iY + d0 * tangent[1];
//        if (isAnyWall(param.flag(posX, posY))) {
//            continue;
//        }
//        h[i] = 0.0;
//        for (int d = -minLim; d <= maxLim; d++) {
//            plint nextX = posX + d * integrationVector[0];
//            plint nextY = posY + d * integrationVector[1];
//            h[i] += param.volumeFraction(nextX, nextY);
//        }
//    }

//    // Extrapolate on walls. (No contact angle algorithm).
//    for (int i = 0; i < 3; i++) {
//        if (util::isZero(h[i] + (T) 1)) {
//            h[i] = h[1];
//        }
//    }
//}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceGeometry2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  static T degToRad = (T)3.14159265358979323844L / (T)180;

  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  Array<T, 2> zeroVector((T)0, (T)0);

  if (!useContactAngle) {
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
        if (param.flag(iX, iY) == interface) {
          if (param.isBoundary(iX, iY)) {
            param.curvature(iX, iY) = 0.0;
            param.setNormal(iX, iY, zeroVector);
            continue;
          }
          // Locally smooth the volume fraction to compute an estimate of the
          // normal.
          T svfcp = param.smooth(*param.volumeFractionP(),
                                 param.volumeFractionOffset(), iX, iY);
          T svfx0 = !isAnyWall(param.flag(iX - 1, iY))
                        ? param.smooth(*param.volumeFractionP(),
                                       param.volumeFractionOffset(), iX - 1, iY)
                        : svfcp;
          T svfx1 = !isAnyWall(param.flag(iX + 1, iY))
                        ? param.smooth(*param.volumeFractionP(),
                                       param.volumeFractionOffset(), iX + 1, iY)
                        : svfcp;
          T svfy0 = !isAnyWall(param.flag(iX, iY - 1))
                        ? param.smooth(*param.volumeFractionP(),
                                       param.volumeFractionOffset(), iX, iY - 1)
                        : svfcp;
          T svfy1 = !isAnyWall(param.flag(iX, iY + 1))
                        ? param.smooth(*param.volumeFractionP(),
                                       param.volumeFractionOffset(), iX, iY + 1)
                        : svfcp;

          // Compute a normalized grad(VF) (inward-pointing normal).
          Array<T, 2> gradVF;
          gradVF[0] = 0.5 * (svfx1 - svfx0);
          gradVF[1] = 0.5 * (svfy1 - svfy0);
          T norm_gradVF = norm(gradVF);
          if (util::isZero(norm_gradVF)) {
            param.curvature(iX, iY) = 0.0;
            param.setNormal(iX, iY, zeroVector);
            continue;
          }
          gradVF /= norm_gradVF;

          T abs0 = std::fabs(gradVF[0]);
          T abs1 = std::fabs(gradVF[1]);

          int integrationDirection = 1;

          if (abs0 > abs1) {
            integrationDirection = 0;
          }

          T h[3];
          computeHeights2D(param, integrationDirection, iX, iY, h);

          T dh0 = 0.5 * (h[2] - h[0]);

          T dh00 = h[2] - 2.0 * h[1] + h[0];

          T value = -(dh00) / std::pow((T)1.0 + dh0 * dh0, (T)1.5);

          param.curvature(iX, iY) = value;

          T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
          Array<T, 2> normal;
          if (integrationDirection == 0) {
            normal = Array<T, 2>(sgn, -dh0);
          } else {
            normal = Array<T, 2>(-dh0, sgn);
          }
          T norm_normal = norm(normal);
          if (util::isZero(norm_normal)) {
            param.setNormal(iX, iY, zeroVector);
          } else {
            param.setNormal(iX, iY, normal / norm_normal);
          }
        } else {
          param.curvature(iX, iY) = 0.0;
          param.setNormal(iX, iY, zeroVector);
        }
      }
    }
  } /*else { // Use contact angles. Mehdi: disabled contact angle for now
      // First compute the flags.
      ScalarField2D<int> *interfaceFlag = getInterfaceFlags(domain, param);

      // New contact angle algorithm. This algorithm does not treat properly
either the adjacent cells
       //                              or the free-slip wall cells.

      // First loop over all the regular and adjacent interface cells and
calculate the curvature and the normal vectors.
      // When the appropriate algorithm for the adjacent cells is implemented,
they must be removed from these loops. for (plint iX=domain.x0; iX<=domain.x1;
++iX) { plint i = iX - domain.x0; for (plint iY=domain.y0; iY<=domain.y1; ++iY)
{ plint j = iY - domain.y0;

                  if (interfaceFlag->get(i, j) == regular ||
interfaceFlag->get(i, j) == adjacent) { if (param.isBoundary(iX, iY)) {
                          param.curvature(iX, iY) = 0.0;
                          param.setNormal(iX, iY, zeroVector);
                          continue;
                      }
                      // Locally smooth the volume fraction to compute an
estimate of the normal. T svfcp = param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), iX, iY); T svfx0 = !isAnyWall(param.flag(iX - 1,
iY)) ? param.smooth(*param.volumeFractionP(), param.volumeFractionOffset(), iX -
1, iY) : svfcp; T svfx1 = !isAnyWall(param.flag(iX + 1, iY)) ?
                          param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), iX + 1, iY) : svfcp; T svfy0 =
!isAnyWall(param.flag(iX, iY - 1)) ? param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), iX, iY - 1) : svfcp; T svfy1 =
!isAnyWall(param.flag(iX, iY + 1)) ? param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), iX, iY + 1) : svfcp;


                      // Compute a normalized grad(VF) (inward-pointing normal).
                      Array<T,2> gradVF;
                      gradVF[0] = 0.5 * (svfx1 - svfx0);
                      gradVF[1] = 0.5 * (svfy1 - svfy0);
                      T norm_gradVF = norm(gradVF);
                      if (util::isZero(norm_gradVF)) {
                          param.curvature(iX, iY) = 0.0;
                          param.setNormal(iX, iY, zeroVector);
                          continue;
                      }
                      gradVF /= norm_gradVF;

                      T abs0 = std::fabs(gradVF[0]);
                      T abs1 = std::fabs(gradVF[1]);

                      int integrationDirection=2;

                      if (abs0 > abs1) {

                              integrationDirection = 0;
                          }
                          else {
                              integrationDirection = 1;
                          }


                      T h[3];
                      computeHeights2D(param, integrationDirection, iX, iY, h);

                      T dh0 = 0.5 * (h[2] - h[0]);

                      T dh00 = h[2] - 2.0 * h[1] + h[0];

                      T value = -(dh00) /
                          std::pow((T)1.0 + dh0 * dh0, (T)1.5);

                      param.curvature(iX, iY) = value;

                      T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
                      Array<T,2> normal;
                      if (integrationDirection == 0) {
                          normal = Array<T,2>(sgn, -dh0);
                      } else  {
                          normal = Array<T,2>(-dh0, sgn);
                      }

                      T norm_normal = norm(normal);
                      if (util::isZero(norm_normal)) {
                          param.setNormal(iX, iY, zeroVector);
                      } else {
                          param.setNormal(iX, iY, normal / norm_normal);
                      }
                  } else {
                      param.curvature(iX, iY) = 0.0;
                      param.setNormal(iX, iY, zeroVector);
                  }

          }
      }

      // Then loop over all the contact-line interface cells and calculate the
curvature and the
      // normal vectors according to the specified contact angle.
      T contactAngleRad = contactAngle * degToRad;
      T tanContactAngle = std::tan(contactAngleRad);

      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
          plint i = iX - domain.x0;
          for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
              plint j = iY - domain.y0;

                  if (interfaceFlag->get(i, j) == contactLine) {
                      if (param.isBoundary(iX, iY)) {
                          param.curvature(iX, iY) = 0.0;
                          param.setNormal(iX, iY, zeroVector);
                          continue;
                      }

                      // First decide where is the wall.
                      int numWallCells = 0;
                      // Computation of the inward-pointing wall normal (not
unitary). Array<int,2> inwardWallNormal(0, 0); for (int dx = -1; dx < 2; dx++) {
                          for (int dy = -1; dy < 2; dy++) {
                                      if (isAnyWall(param.flag(iX + dx, iY +
dy))) { inwardWallNormal += Array<int,2>(-dx, -dy); numWallCells++;
                                  }

                          }
                      }
                      PLB_ASSERT(numWallCells != 0);
#ifdef PLB_DEBUG
                      int norm2_inwardWallNormal = inwardWallNormal[0] *
inwardWallNormal[0] + inwardWallNormal[1] * inwardWallNormal[1];

#endif
                      PLB_ASSERT(norm2_inwardWallNormal != 0);

                      int iWallNormalDirection;
                      // The inwardWallNormal is aligned with one axis.
                      if (inwardWallNormal[0] != 0 && inwardWallNormal[1] == 0)
{ iWallNormalDirection = 0; } else if (inwardWallNormal[1] != 0 &&
inwardWallNormal[0] == 0) { iWallNormalDirection = 1;
                      }
                       else {
                          // The inwardWallNormal is not aligned with one axis.
                          Array<int,2> sumDirection[2];
                          sumDirection[0] = inwardWallNormal[0] == 0 ?
Array<int,2>( 0,  0) : (inwardWallNormal[0] >  0 ? Array<int,2>( 1,  0) :
                                                                       Array<int,2>(-1,
0)); sumDirection[1] = inwardWallNormal[1] == 0 ? Array<int,2>( 0,  0) :
                                           (inwardWallNormal[1] >  0 ?
Array<int,2>( 0,  1) : Array<int,2>( 0, -1));

                          T sum[2] = { std::numeric_limits<T>::max(),
                                       std::numeric_limits<T>::max()
                                        };
                          for (int iSum = 0; iSum < 2; iSum++) {
                              if (sumDirection[iSum][0] + sumDirection[iSum][1]
!= 0) { sum[iSum] = 0.0; for (int d = 0; d <= 3; d++) { plint posX = iX + d *
sumDirection[iSum][0]; plint posY = iY + d * sumDirection[iSum][1]; if
(!isAnyWall(param.flag(posX, posY))) { sum[iSum] += param.volumeFraction(posX,
posY);
                                      }
                                  }
                              }
                          }

                          // The wall normal direction is the direction of the
smallest sum. if (sum[0] < sum[1]) {

                                  iWallNormalDirection = 0;

                              else {
                                  iWallNormalDirection = 1;
                              }
                          }
                          // sum[0] >= sum[1]

                      }

                      // Reset the inward wall normal to be unitary and to
contain information on the direction. inwardWallNormal[0] = iWallNormalDirection
!= 0 ? 0 : (inwardWallNormal[0] > 0 ? 1 : -1); inwardWallNormal[1] =
iWallNormalDirection != 1 ? 0 : (inwardWallNormal[1] > 0 ? 1 : -1);

                      // Define a wall normal that shows only the wall normal
axis. Array<int,2> wallNormal; wallNormal[0] = iWallNormalDirection == 0 ? 1 :
0; wallNormal[1] = iWallNormalDirection == 1 ? 1 : 0;

                      // Compute the wall tangent vectors.
                      int iWallTangentDirection0 = iWallNormalDirection == 0 ? 1
: 0; Array<int,2> wallTangent0; wallTangent0[0] = iWallTangentDirection0 == 0 ?
1 : 0; wallTangent0[1] = iWallTangentDirection0 == 1 ? 1 : 0;


                      // Locally smooth the volume fraction to compute an
estimate of the 2D normal. T svfcp = param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), iX, iY); plint posX, posY; posX = iX -
wallTangent0[0]; posY = iY - wallTangent0[1]; T svf00 =
!isAnyWall(param.flag(posX, posY)) ? param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), posX, posY) : svfcp; posX = iX + wallTangent0[0];
                      posY = iY + wallTangent0[1];
                      T svf01 = !isAnyWall(param.flag(posX, posY)) ?
                          param.smooth(*param.volumeFractionP(),
param.volumeFractionOffset(), posX, posY) : svfcp;


                      // Compute a normalized 2D grad(VF) (inward-pointing 2D
normal). Array<T,2> gradVF2D; gradVF2D[0] = 0.5 * (svf01 - svf00); gradVF2D[1] =
0.5 * (svf11 - svf10); T norm_gradVF2D = norm(gradVF2D); if
(util::isZero(norm_gradVF2D)) { param.curvature(iX, iY) = 0.0;
                          param.setNormal(iX, iY, zeroVector);
                          continue;
                      }
                      gradVF2D /= norm_gradVF2D;


                      int integrationDirection2D = 0; // wallTangent0.


                      T h2D[3];
                      computeHeights2DC(param, wallTangent0, wallTangent1,
integrationDirection2D, iX, iY, h2D);

                      T dh2D = 0.5 * (h2D[2] - h2D[0]);

                      T sgn2D = -gradVF2D[integrationDirection2D] < 0.0 ? -1.0
: 1.0; Array<T,2> normal2D; normal2D = Array<T,2>(sgn2D, -dh2D); T norm_normal2D
= norm(normal2D); if (util::isZero(norm_normal2D)) { param.curvature(iX, iY) =
0.0; param.setNormal(iX, iY, zeroVector); continue;
                      }

                      Array<T,2> normal; // 2D outward unit normal.
                      T wallNormalComponent = norm_normal2D / tanContactAngle;
                      normal[0] = normal2D[0] * wallTangent0[0] +
wallNormalComponent * wallNormal[0]; normal[1] = normal2D[0] * wallTangent0[1] +
wallNormalComponent * wallNormal[1]; T norm_normal = norm(normal); if
(util::isZero(norm_normal)) { param.setNormal(iX, iY, zeroVector); } else {
                          param.setNormal(iX, iY, normal / norm_normal);
                      }

                      // Now compute the curvature.
                      // First compute the 2D height functions.
                      int integrationDirection;

                      integrationDirection = wallTangent0[0] != 0 ? 0 : 1;

                      T h[3];
                      computeHeights2D(param, integrationDirection, iX, iY, h);

                      // Determine the orientation of the elements of h.
                      int iTangentDirection0 = integrationDirection == 0 ? 1 :
0; Array<int,2> tangent0; tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
                      tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;
                      int i0 = -1;
                      if (inwardWallNormal[0] == tangent0[0] &&
                          inwardWallNormal[1] == tangent0[1]
                          ) {
                          i0 = 0;
                      } else if (inwardWallNormal[0] == -tangent0[0] &&
                                 inwardWallNormal[1] == -tangent0[1]
                                 ) {
                          i0 = 2;

                      } else {
                          PLB_ASSERT(false);
                      }

                      Array<T,2> v1, v2; // In the wallTangent0, wallTangent1
base. v1[0] = std::fabs(normal2D[0]); v1[1] = std::fabs(normal2D[1]); if
(integrationDirection2D == 0) { v2[0] = 1.0; v2[1] = 0.0; } else { v2[0] = 0.0;
                          v2[1] = 1.0;
                      }
                      T cosAlpha = std::cos(angleBetweenVectors(v1, v2));
                      T correction = 1.0 / (tanContactAngle * cosAlpha);

                      if (i0 != -1) {

                              h[i0] = h[1] + correction;

                      }
                       else {
                          PLB_ASSERT(false);
                      }

                      T dh0 = 0.5 * (h[2] - h[0]);

                      T dh00 = h[2] - 2.0 * h[1] + h[0];


                      T value = -(dh00) /
                          std::pow((T)1.0 + dh0 * dh0, (T)1.5);

                      param.curvature(iX, iY) = value;
                  }

          }
      }

      delete interfaceFlag;
  }*/
}

/* *************** Class FreeSurfaceComputeCurvature2D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeCurvature2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  static T degToRad = (T)3.14159265358979323844L / (T)180;

  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  Dot2D location = param.absOffset();

  int smooth = 0;  // 0 for no normal vector field smoothing.
  // 1 for smoothing the normal vector field once before imposing the contact
  // angle. 2 for smoothing the normal vector field once after imposing the
  // contact angle.

  TensorField2D<T, 2> *normalForContactAngle = 0;
  Dot2D ofsNCA(0, 0);
  plint w = 0;

  if (smooth == 0 || smooth == 2) {
    // Tensor field to hold a temporary vector field of unit normals. (Include
    // also a w-cell layer around "domain".)
    w = smooth == 0 ? 1 : 2;
    plint nx = domain.getNx() + 2 * w;
    plint ny = domain.getNy() + 2 * w;
    normalForContactAngle = new TensorField2D<T, 2>(nx, ny);
    ofsNCA = Dot2D(-domain.x0 + w, -domain.y0 + w);
    for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
      plint i = iX + ofsNCA.x;
      for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
        plint j = iY + ofsNCA.y;

        normalForContactAngle->get(i, j) = param.getNormal(iX, iY);

        param.curvature(iX, iY) = 0.0;
      }
    }
  } else {
    // Tensor field to hold a temporary vector field of unit normals. (Include
    // also a w-cell layer around "domain".) This normal field is smoothed once
    // (early tests show that if the normal field is smoothed twice, then the
    // results are pretty bad).
    w = 1;
    plint nx = domain.getNx() + 2 * w;
    plint ny = domain.getNy() + 2 * w;
    normalForContactAngle = new TensorField2D<T, 2>(nx, ny);
    ofsNCA = Dot2D(-domain.x0 + w, -domain.y0 + w);
    for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
      plint i = iX + ofsNCA.x;
      for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
        plint j = iY + ofsNCA.y;

        normalForContactAngle->get(i, j) = param.template lbmSmooth<
            T, descriptors::AdvectionDiffusionD2Q5Descriptor>(
            *param.normalP(), param.normalOffset(), iX, iY);

        param.curvature(iX, iY) = 0.0;
      }
    }
  }

  // Enforce contact angles.
  if (useContactAngle) {
    for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
      plint i = iX + ofsNCA.x;
      for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
        plint j = iY + ofsNCA.y;

        T curv = 0.0;

        // If we want to operate on the "sharp free surface", which is defined
        // by the lattice cells with the "interface" flag, then the next clause
        // must be used.
        // if (param.flag(iX, iY) == interface) {
        // If we want a "diffuse interface view", then instead of the above
        // clause, we need to have one of the following, which operate on a
        // "thick free surface" that is an outcome of the smoothing proceedure
        // (diffusion) performed when computing the field of unit normals. There
        // are two ways to do this: one is to use the original normal
        // ("getNormal"), which is already "thick", and the other way is to use
        // the "normalForContactAngle" which is localy smoothed one more time if
        // "smooth = 1".
        if (!isAnyWall(param.flag(iX, iY)) &&
            !util::isZero(normSqr(param.getNormal(iX, iY)))) {
          // if (!isAnyWall(param.flag(iX, iY)) &&
          // !util::isZero(normSqr(normalForContactAngle->get(i, j)))) {
          int isaContactAngleCell = 0;
          int numWallCells = 0;
          // Computation of the inward-pointing wall normal (towards the fluid).
          Array<int, 2> tmpWallNormal(0, 0);
          plint h = 1;  // This cannot be greater than 1 when w = 2!
          for (plint dx = -h; dx <= h; dx++) {
            for (plint dy = -h; dy <= h; dy++) {
              int flg = param.flag(iX + dx, iY + dy);
              if (isWall(flg)) {  // Here we want only the no-slip walls, not
                                  // the free-slip ones.
                tmpWallNormal += Array<int, 2>(-dx, -dy);
                numWallCells++;
              }
            }
          }
          Array<T, 2> wallNormal;
          if (numWallCells != 0) {
            int norm2tmpWallNormal = tmpWallNormal[0] * tmpWallNormal[0] +
                                     tmpWallNormal[1] * tmpWallNormal[1];
            if (norm2tmpWallNormal != 0) {
              T tmpNormWallNormal = std::sqrt((T)norm2tmpWallNormal);
              wallNormal[0] = (T)tmpWallNormal[0] / tmpNormWallNormal;
              wallNormal[1] = (T)tmpWallNormal[1] / tmpNormWallNormal;

              isaContactAngleCell = 1;
            }
          }
          if (isaContactAngleCell) {
            // Here again, there are two possible normals one can use: the
            // original normal ("getNormal"), which is already "thick", and the
            // "normalForContactAngle" which is localy smoothed one more time if
            // "smooth = 1".
            Array<T, 2> normal = param.getNormal(iX, iY);
            // Array<T,2> normal = normalForContactAngle->get(i, j);
            T cosPhi = dot(normal, wallNormal);
            Array<T, 2> wallTangent = normal - cosPhi * wallNormal;
            T normWallTangent = norm(wallTangent);
            if (!util::isZero(normWallTangent)) {
              wallTangent /= normWallTangent;
            } else {
              // This option to not enforce a contact angle when the interface
              // normal is parallel to the wall normal must be revised in the
              // framework of long thin films.
              continue;
            }

            // Enforce the contact angle by taking at most two actions:
            //
            // 1) Add a penalty-like contribution to the curvature, to
            // implicitly correct
            //    towards the equilibrium contact angle. For the form of the
            //    penalty-like term, see: Attar et al, "Lattice Boltzmann method
            //    for dynamic wetting problems", Journal of Colloid and
            //    Interface Science, 335 (2009) 84-93.
            //
            //    When operating in very low resolution, one should consider not
            //    using this penalty term, for better accuracy.
            //
            // 2) Modify the free surface normal at the local interface cell, so
            // that
            //    the next time the curvature is computed, the equilibrium
            //    contact angle is considered.
            //
            // Both of these actions impose the contact angle indirectly,
            // through Young's law which relates the curvature, the surface
            // tension and the pressure drop on the free surface.

            T thetaEq =
                (contactAngleFunction == 0
                     ? contactAngle
                     : contactAngleFunction(iX + location.x, iY + location.y) *
                           degToRad);  // Equilibrium contact angle in radians.

            // curv = cosPhi - std::cos(thetaEq); // At convergence, this term
            // is zero.
            normal = std::cos(thetaEq) * wallNormal +
                     std::sin(thetaEq) * wallTangent;

            // This is yet another way, that does not try to impose the
            // equilibrium contact angle, but an extrapolation of it, taking
            // under consideration the current value of the dynamic contact
            // angle (phi): theta = (T) 2 * thetaEq - phi. After convergence is
            // achieved, then: phi = thetaEq, and theta = thetaEq.
            //
            // Another way is that one can impose the curvature as:
            // curv = cosPhi - std::cos(thetaEq), but the normal as:
            // normal = std::cos(theta) * wallNormal + std::sin(theta) *
            // wallTangent, with: theta = (T) 2 * thetaEq - phi. This is a
            // mixture of the two previously discussed methods. Sometimes it
            // gives nice results (it enforces the contact angle more strongly).

            // T phi = angleBetweenVectors(normal, wallNormal);
            // T theta = (T) 2 * thetaEq - phi;

            // curv = std::cos(theta) - std::cos(thetaEq); // At convergence,
            // this term is zero. normal = std::cos(theta) * wallNormal +
            // std::sin(theta) * wallTangent;

            normalForContactAngle->get(i, j) = normal;
          }
        }

        param.curvature(iX, iY) = curv;
      }
    }
  }

  TensorField2D<T, 2> *normalForCurvature = 0;
  Dot2D ofsNC(0, 0);

  if (smooth == 2) {
    // Tensor field to hold a temporary vector field of unit normals. (Include
    // also a t-cell layer around "domain".) This normal field is smoothed once
    // (early tests show that if the normal field is smoothed twice, then the
    // results are pretty bad).
    plint t = 1;
    plint nx = domain.getNx() + 2 * t;
    plint ny = domain.getNy() + 2 * t;
    normalForCurvature = new TensorField2D<T, 2>(nx, ny);
    ofsNC = Dot2D(-domain.x0 + t, -domain.y0 + t);
    for (plint iX = domain.x0 - t; iX <= domain.x1 + t; ++iX) {
      plint i = iX + ofsNC.x;
      for (plint iY = domain.y0 - t; iY <= domain.y1 + t; ++iY) {
        plint j = iY + ofsNC.y;

        normalForCurvature->get(i, j) = param.template lbmSmooth<
            T, descriptors::AdvectionDiffusionD2Q5Descriptor>(
            *normalForContactAngle, ofsNCA, iX, iY);
      }
    }
  } else {
    normalForCurvature = normalForContactAngle;
    ofsNC = ofsNCA;
  }

  // Compute the curvature as the divergence of the vector field of unit
  // normals.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      T curv = 0.0;

      if (param.flag(iX, iY) != interface) {
        param.curvature(iX, iY) = curv;
        continue;
      }

      /*
          int useLB = 1;
          for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
              plint nextX = iX + Descriptor<T>::c[iPop][0];
              plint nextY = iY + Descriptor<T>::c[iPop][1];
              if (isAnyWall(param.flag(nextX, nextY))) {
                  useLB = 0;
                  break;
              }
          }

          if (useLB) {
          */
      // Compute the divergence of the normal vector field "the lattice
      // Boltzmann way".
      curv = param.template lbmComputeDivergence<T, Descriptor>(
          *normalForCurvature, ofsNC, iX, iY);
      /*
          } else {
              // Compute the divergence of the normal vector field with finite
         differences on the interface cells
              // excluding wall cells.
              plint h = 1;    // This must be 1 because above we included only a
         1-cell layer around "domain". curv =
         param.computeDivergence(*normalForCurvature, ofsNC, h, iX, iY);
          }
          */

      // We restrict the radius of curvature to be more always >= 0.5, in
      // lattice units. A smaller radius makes no sense anyway, numerically
      // speaking, and in this way we avoid problems of the "division by zero"
      // kind. (radius = 2/curvature)
      // if (std::fabs(curv) > (T) 4) {
      //    if (curv < (T) 0) {
      //        curv = -4.0;
      //    }
      //    else {
      //        curv = 4.0;
      //    }
      //}

      param.curvature(iX, iY) +=
          curv;  // We add in order to include the "penalty-like" term if any.
    }
  }

  delete normalForContactAngle;
  if (smooth == 2) {
    delete normalForCurvature;
  }
}

/* *************** Class FreeSurfaceMassChange2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceMassChange2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // This loop updates the mass, summarizing  Eq. 6/7, and Eq.8, in
  // the N. Thuerey e.a. technical report "Interactive Free Surface Fluids
  // with the Lattice Boltzmann Method".
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      Cell<T, Descriptor> &cell = param.cell(iX, iY);
      int flag = param.flag(iX, iY);
      if (isFullWet(flag)) {
        freeSurfaceTemplates<T, Descriptor>::massExchangeFluidCell(param, iX,
                                                                   iY);
      } else if (flag == interface) {
        for (plint iPop = 0; iPop < D::q; ++iPop) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];
          int nextFlag = param.flag(nextX, nextY);
          plint opp = indexTemplates::opposite<D>(iPop);
          // Calculate mass at time t+1 on interface cell --> eq 7 Thurey's
          // paper.
          if (isFullWet(nextFlag)) {
            param.mass(iX, iY) += (cell[opp] - param.cell(nextX, nextY)[iPop]);
          } else if (nextFlag == interface) {
            param.mass(iX, iY) += (cell[opp] - param.cell(nextX, nextY)[iPop]) *
                                  0.5 *
                                  (param.volumeFraction(nextX, nextY) +
                                   param.volumeFraction(iX, iY));
          }
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceCompletion2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceCompletion2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  using namespace freeSurfaceFlag2D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;

  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // In this data processor, populations are both written locally and read
  // non-locally. To guarantee data consistency, a first loop makes only read
  // accesses and stores the necessary information into the list
  // neighborOppositePop. A second loop reads from this list and assigns values
  // to populations.
  std::map<Node, Array<T, D::q>> neighborOppositePop;
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      // This is the old form of the completion scheme. There is this extra
      // condition mentioned by Thurey which has to do with the normal to the
      // interface. We found that this condition is responsible for an
      // instability when one increases both the spatial and temporal resolution
      // while respecting the diffusive limit in the presence of surface
      // tension. We also found that it causes an instability at the simple test
      // case of a fluid sphere which is subject to surface tension but not to
      // any other force. This sphere should remain still, but in the presence
      // of this condition it starts moving.
      /*
          if (param.flag(iX,iY) == interface) {
              // Here we are on an interface node. The entire set of fi's is
         reconstructed.
              // The normal is recomputed as in eq. 10 of Thurey's paper.
              Array<T,2> normalToInterface;
              normalToInterface = param.getNormal(iX, iY);

              bool needsModification = false;
              Array<T,D::q> savedPop;
              savedPop[0] = -2.;
              for(plint iPop=1; iPop < D::q; ++iPop )
              {
                  // This is one of the tricky points of the code
                  // we have to decide if the f_is from the neighborhood
                  // have to be re-update by using the Thurey's rule, which
                  // states that f_i's coming from nearest neighs. that are
         empty cells,
                  // have to be re-updated.
                  // I like the eq.   f^{in}_i(x,t+dt) = f^{out}_i(x-e_i,t);
                  // This eq. makes me think that the neigh. that I have to
         check
                  // (to control is status e.g. empty or fluid ?) has to be
         pos-c_i plint prevX = iX-D::c[iPop][0]; plint prevY = iY-D::c[iPop][1];

                  plint opp = indexTemplates::opposite<D>(iPop);
                  T scalarProduct = D::c[opp][0]*normalToInterface[0] +
                                    D::c[opp][1]*normalToInterface[1];

                  // Should I also change particle distribution function coming
         from
                  // bounceBack nodes? Well ideally no ... but there is for sure
         some
                  // cell configuration where these f_is are not well defined
         because
                  // they are probably coming from empty cells

                  // If the f_i[iPop] would be streamed from an empty cell, or
         whenever the scalar product is positive. if ( scalarProduct > 0 ||
         param.flag(prevX,prevY) == empty ||
                       isAnyWall(param.flag(prevX,prevY,prevZ)) )
                  {
                      savedPop[iPop] = param.cell(prevX,prevY)[opp];
                      needsModification = true;
                  }
                  else {
                      savedPop[iPop] = (T)-2.;
                  }
              }
              if (needsModification) {
                  neighborOppositePop.insert(std::pair<Node,Array<T,D::q>
         >(Node(iX,iY), savedPop));
              }
          }
          */

      if (param.flag(iX, iY) == interface) {
        // Here we are on an interface node. The entire set of fi's is
        // reconstructed.
        bool needsModification = false;
        Array<T, D::q> savedPop;
        savedPop[0] = -2.;
        for (plint iPop = 1; iPop < D::q; ++iPop) {
          // This is one of the tricky points of the code
          // we have to decide if the f_is from the neighborhood
          // have to be re-update by using the Thurey's rule, which
          // states that f_i's coming from nearest neighs. that are empty cells,
          // have to be re-updated.
          // I like the eq.   f^{in}_i(x,t+dt) = f^{out}_i(x-e_i,t);
          // This eq. makes me think that the neigh. that I have to check
          // (to control is status e.g. empty or fluid ?) has to be pos-c_i
          plint prevX = iX - D::c[iPop][0];
          plint prevY = iY - D::c[iPop][1];

          plint opp = indexTemplates::opposite<D>(iPop);

          // Should I also change particle distribution function coming from
          // bounceBack nodes? Well ideally no ... but there is for sure some
          // cell configuration where these f_is are not well defined because
          // they are probably coming from empty cells

          // If the f_i[iPop] would be streamed from an empty cell
          if (isEmpty(param.flag(prevX, prevY)) ||
              isAnyWall(param.flag(prevX, prevY))) {
            savedPop[iPop] = param.cell(prevX, prevY)[opp];
            needsModification = true;
          } else {
            savedPop[iPop] = (T)-2.;
          }
        }
        if (needsModification) {
          neighborOppositePop.insert(
              std::pair<Node, Array<T, D::q>>(Node(iX, iY), savedPop));
        }
      }
    }
  }

  typename std::map<Node, Array<T, D::q>>::const_iterator nodes =
      neighborOppositePop.begin();
  for (; nodes != neighborOppositePop.end(); ++nodes) {
    Node node = nodes->first;
    plint iX = node[0];
    plint iY = node[1];
    Array<T, D::q> neighborOppPop = nodes->second;
    for (plint iPop = 1; iPop < D::q; ++iPop) {
      if (neighborOppPop[iPop] > (T)-1.) {
        // Velocity is simply taken from the previous time step.
        Array<T, 2> j = param.getMomentum(iX, iY);
        T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
        // Remember: the value of pressure on an interface node has been set in
        // F, and is equal to the ambient pressure for a
        // single free-surface fluid, or in the case of a binary pressure, an
        // averaged value.
        T rhoBar = Descriptor<T>::rhoBar(param.getDensity(iX, iY));
        T feq_i = param.cell(iX, iY).computeEquilibrium(iPop, rhoBar, j, jSqr);
        plint opp = indexTemplates::opposite<D>(iPop);
        T feq_opp_i =
            param.cell(iX, iY).computeEquilibrium(opp, rhoBar, j, jSqr);
        param.cell(iX, iY)[iPop] = feq_i + feq_opp_i - neighborOppPop[iPop];
      }
    }
  }
}

/* *************** Class FreeSurfaceMacroscopic2D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceMacroscopic2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  T lostMass = param.getSumLostMass();
  plint numInterfaceCells = param.getNumInterfaceCells();
  T massPerCell = T();
  if (numInterfaceCells > 0) {
    massPerCell = lostMass / (T)numInterfaceCells;
  }

  // Save macroscopic fields in external scalars and update the mass and the
  // volume-fraction.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isWet(param.flag(iX, iY))) {
        T rhoBar;
        Array<T, 2> j;
        momentTemplates<T, Descriptor>::get_rhoBar_j(param.cell(iX, iY), rhoBar,
                                                     j);
        T density = Descriptor<T>::fullRho(rhoBar);
        param.setDensity(iX, iY, density);

        if (param.flag(iX, iY) == interface) {
          param.mass(iX, iY) += massPerCell;
          T newDensity = param.outsideDensity(iX, iY);
          param.volumeFraction(iX, iY) = param.mass(iX, iY) / newDensity;
          // On interface cells, adjust the pressure to the ambient pressure.
          param.setDensity(iX, iY, newDensity);
          if (!incompressibleModel) {
            j *= newDensity / density;
          }
        } else if (isFullWet(param.flag(iX, iY))) {
          param.volumeFraction(iX, iY) = T(1);
        }

        param.setMomentum(iX, iY, j);
      }
    }
  }
}

/* *************** Class FreeSurfaceAddSurfaceTension2D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddSurfaceTension2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // Save macroscopic fields in external scalars and add the surface tension
  // effect.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.flag(iX, iY) == interface) {
        T density = param.getDensity(iX, iY);
        T newDensity = density;
        // pcout << "old dens surface: " << newDensity << std::endl;
        // Stored curvature is computed to be twice the mean curvature.
        // pcout << surfaceTension * param.curvature(iX,iY) * D::invCs2 <<
        // std::endl;
        newDensity += surfaceTension * param.curvature(iX, iY) * D::invCs2;
        // pcout << "new dens surface: " << newDensity << std::endl;

        param.volumeFraction(iX, iY) = param.mass(iX, iY) / newDensity;
        // On interface cells, adjust the pressure to incorporate surface
        // tension.
        param.setDensity(iX, iY, newDensity);
        if (!incompressibleModel) {
          Array<T, 2> j = param.getMomentum(iX, iY);
          Array<T, 2> newJ = j * newDensity / density;
          param.setMomentum(iX, iY, newJ);
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceAddSurfaceTensionFromScalarField2D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddSurfaceTensionFromScalarField2D<T, Descriptor>::
    processGenericBlocks(Box2D domain,
                         std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  using namespace freeSurfaceFlag2D;

  std::vector<AtomicBlock2D *> freeSurfaceBlocks(atomicBlocks.begin(),
                                                 atomicBlocks.end() - 1);
  FreeSurfaceProcessorParam2D<T, Descriptor> param(freeSurfaceBlocks);

  ScalarField2D<T> *surfaceTensionField =
      dynamic_cast<ScalarField2D<T> *>(*(atomicBlocks.end() - 1));
  PLB_ASSERT(surfaceTensionField);

  Dot2D ofsSTF =
      computeRelativeDisplacement(*atomicBlocks[0], *surfaceTensionField);

  // Save macroscopic fields in external scalars and add the surface tension
  // effect.
  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.flag(iX, iY) == interface) {
        T density = param.getDensity(iX, iY);
        T newDensity = density;
        T surfaceTension =
            surfaceTensionField->get(iX + ofsSTF.x, iY + ofsSTF.y);
        // Stored curvature is computed to be twice the mean curvature.
        newDensity += surfaceTension * param.curvature(iX, iY) * D::invCs2;
        param.volumeFraction(iX, iY) = param.mass(iX, iY) / newDensity;
        // On interface cells, adjust the pressure to incorporate surface
        // tension.
        param.setDensity(iX, iY, newDensity);
        if (!incompressibleModel) {
          Array<T, 2> j = param.getMomentum(iX, iY);
          Array<T, 2> newJ = j * newDensity / density;
          param.setMomentum(iX, iY, newJ);
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceStabilize2D ********************************
 */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceStabilize2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (param.flag(iX, iY) == interface) {
        Cell<T, Descriptor> &cell = param.cell(iX, iY);
        T oldRhoBar;
        Array<T, 2> oldJ;
        momentTemplates<T, Descriptor>::get_rhoBar_j(cell, oldRhoBar, oldJ);
        T oldJsqr = normSqr(oldJ);
        T newDensity = param.getDensity(iX, iY);
        T newRhoBar = Descriptor<T>::rhoBar(newDensity);
        Array<T, 2> newJ = param.getMomentum(iX, iY);
        T newJsqr = normSqr(newJ);
        for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
          T oldEq = cell.getDynamics().computeEquilibrium(iPop, oldRhoBar, oldJ,
                                                          oldJsqr);
          T newEq = cell.getDynamics().computeEquilibrium(iPop, newRhoBar, newJ,
                                                          newJsqr);
          cell[iPop] += newEq - oldEq;
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceComputeInterfaceLists2D
 * ******************************************* */

template <typename T, template <typename> class Descriptor>
T FreeSurfaceComputeInterfaceLists2D<T, Descriptor>::kappa = 1.e-3;

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeInterfaceLists2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);
  using namespace freeSurfaceFlag2D;

  param.massExcess().clear();
  param.interfaceToFluid().clear();
  param.interfaceToEmpty().clear();
  param.emptyToInterface().clear();

  // interfaceToFluid needs to be computed in bulk+2.
  for (plint iX = domain.x0 - 2; iX <= domain.x1 + 2; ++iX) {
    for (plint iY = domain.y0 - 2; iY <= domain.y1 + 2; ++iY) {
      Node node(iX, iY);
      // Eq. 11 in Thuerey's technical report.
      if (param.flag(iX, iY) == interface) {  // Interface cell.
        if (param.volumeFraction(iX, iY) >
            T(1) + kappa) {  // Interface cell is filled.
          // Elements are added even if they belong to the envelope, because
          // they may be
          //   needed further down in the same data processor.
          param.interfaceToFluid().insert(node);
        } else if (param.volumeFraction(iX, iY) <
                   kappa) {  // Interface cell is empty.
          // Elements are added even if they belong to the envelope, because
          // they may be
          //   needed further down in the same data processor.
          param.interfaceToEmpty().insert(node);
        }
      }
    }
  }

  // Where interface cells have become fluid, neighboring cells must be
  // prevented from
  //   being empty, because otherwise there's no interface cell between empty
  //   and fluid.
  typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
  for (; iEle != param.interfaceToFluid().end(); ++iEle) {
    // The node here may belong to the 1st envelope.
    Node node = *iEle;
    plint iX = node[0];
    plint iY = node[1];

    for (plint iPop = 1; iPop < D::q; ++iPop) {
      plint nextX = iX + D::c[iPop][0];
      plint nextY = iY + D::c[iPop][1];
      Node nextNode(nextX, nextY);

      // If one of my neighbors switches interface->fluid, then I shall be
      // prevented
      //     from switching interface->empty at the same time step.
      if (contained(nextX, nextY, domain.enlarge(1)) &&
          param.flag(nextX, nextY) == interface) {
        param.interfaceToEmpty().erase(nextNode);
      }
      // If one of my neighbors switches interface->fluid and I am empty I shall
      // become
      //   interface.
      else if (contained(nextX, nextY, domain.enlarge(1)) &&
               isEmpty(param.flag(nextX, nextY))) {
        param.emptyToInterface().insert(nextNode);
      }
    }
  }
}

/* *************** Class FreeSurfaceIniInterfaceToAnyNodes2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
FreeSurfaceIniInterfaceToAnyNodes2D<
    T, Descriptor>::FreeSurfaceIniInterfaceToAnyNodes2D(T rhoDefault_)
    : rhoDefault(rhoDefault_) {}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceIniInterfaceToAnyNodes2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  using namespace freeSurfaceFlag2D;

  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // 1. For interface->fluid nodes, update in the flag matrix,
  //   and compute and store mass excess from these cells.
  typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
  for (; iEle != param.interfaceToFluid().end(); ++iEle) {
    Node node = *iEle;

    plint iX = node[0];
    plint iY = node[1];

    if (contained(iX, iY, domain.enlarge(1))) {
      T saveMass = param.mass(iX, iY);
      param.mass(iX, iY) = param.getDensity(iX, iY);
      param.volumeFraction(iX, iY) = (T)1;
      param.flag(iX, iY) = fluid;

      T massExcess = saveMass - param.getDensity(iX, iY);
      param.massExcess().insert(std::pair<Node, T>(node, massExcess));
    }
  }

  // 2. For interface->empty nodes, update in the flag matrix,
  //   and compute and store mass excess from these cells.
  iEle = param.interfaceToEmpty().begin();
  for (; iEle != param.interfaceToEmpty().end(); ++iEle) {
    Node node = *iEle;
    plint iX = node[0];
    plint iY = node[1];

    if (contained(iX, iY, domain.enlarge(1))) {
      // Avoid the case where an empty cell has a fluid neighbor without
      // interface cell between them.
      bool isAdjacentToProtected = false;
      for (plint iPop = 1; iPop < D::q; ++iPop) {
        plint nextX = iX + D::c[iPop][0];
        plint nextY = iY + D::c[iPop][1];
        if (isProtected(param.flag(nextX, nextY))) {
          isAdjacentToProtected = true;
          break;
        }
      }
      if (!isAdjacentToProtected) {
        param.flag(iX, iY) = empty;
        param.attributeDynamics(iX, iY,
                                new NoDynamics<T, Descriptor>(rhoDefault));

        T massExcess = param.mass(iX, iY);
        param.massExcess().insert(std::pair<Node, T>(node, massExcess));

        param.mass(iX, iY) = T();
        param.volumeFraction(iX, iY) = T();
        param.setDensity(iX, iY, rhoDefault);
        param.setForce(iX, iY, Array<T, D::ExternalField::sizeOfForce>::zero());
        param.setMomentum(iX, iY, Array<T, 2>(T(), T()));
        for (plint iPop = 1; iPop < D::q; ++iPop) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];

          // The concurrent read/write on param.flag is not an issue here,
          // because the result in any case is that all adjacent fluid cells
          // have become interface.
          if (param.flag(nextX, nextY) == fluid) {
            param.flag(nextX, nextY) = interface;
          }
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceIniEmptyToInterfaceNodes2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceIniEmptyToInterfaceNodes2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  // In this data processor, density and momentum are potentially read and
  // written
  //   from the same node, because nodes can switch state. The following two
  //   vectors store temporary variables to avoid read/write in undefined order.
  std::vector<T> newDensity(param.emptyToInterface().size());
  std::vector<Array<T, 2>> newMomentum(param.emptyToInterface().size());
  std::fill(newDensity.begin(), newDensity.end(), T());
  std::fill(newMomentum.begin(), newMomentum.end(), Array<T, 2>(T(), T()));

  // Compute density and momentum for cells that will switch state
  // empty->interface.
  //   It is sufficient to do this is bulk+0.
  //   This loop performs read-only access to the lattice.
  plint i = 0;
  typename std::set<Node>::iterator iEle = param.emptyToInterface().begin();
  for (; iEle != param.emptyToInterface().end(); ++iEle, ++i) {
    Node node = *iEle;
    plint iX = node[0];
    plint iY = node[1];

    // If non-bulk elements are left in the list, disregard to avoid accessing
    // undefined neighbors.
    if (contained(iX, iY, domain)) {
      // For initialization of the new cell, compute average density
      //   and momentum on neighbors.
      T averageDensity = T(0);
      Array<T, 2> averageMomentum(T(0), T(0));
      T sumWeights = (T)0;
      for (plint iPop = 1; iPop < D::q; ++iPop) {
        plint nextX = iX + D::c[iPop][0];
        plint nextY = iY + D::c[iPop][1];

        // Warning: it is not accounted for the fact that neighbors can have
        // excess mass. It
        //   might be good to account for this in the future.
        if (isWet(param.flag(nextX, nextY))) {
          T weight = D::t[iPop];
          sumWeights += weight;
          averageDensity += weight * param.getDensity(nextX, nextY);
          averageMomentum += weight * param.getMomentum(nextX, nextY);
        }
      }
      T invSum = T(1) / sumWeights;
      averageDensity *= invSum;
      averageMomentum *= invSum;
      newDensity[i] = averageDensity;
      newMomentum[i] = averageMomentum;
    }
  }

  // Elements that have switched state empty->interface are initialized at
  // equilibrium.
  //   It is sufficient to initialize them in bulk+0.
  //   This loop performs write-only access on the lattice.
  i = 0;
  iEle = param.emptyToInterface().begin();
  for (; iEle != param.emptyToInterface().end(); ++iEle, ++i) {
    Node node = *iEle;

    plint iX = node[0];
    plint iY = node[1];

    // If non-bulk elements are left in the list, disregard to avoid accessing
    // undefined neighbors.
    if (contained(iX, iY, domain)) {
      T averageDensity = newDensity[i];
      Array<T, 2> averageMomentum = newMomentum[i];

      param.attributeDynamics(iX, iY, dynamicsTemplate->clone());

      iniCellAtEquilibrium(param.cell(iX, iY), averageDensity,
                           averageMomentum / averageDensity);
      param.setForce(iX, iY, force);
      // Change density, but leave mass and volumeFraction at 0, as they are
      // later
      //   recomputed (Warning: this is probably correct, but there remains a
      //   small doubt).
      param.setMomentum(iX, iY, averageMomentum);
      param.setDensity(iX, iY, averageDensity);
      param.mass(iX, iY) = T();
      param.volumeFraction(iX, iY) = T();
      param.flag(iX, iY) = interface;
    }
  }
}

/* *************** Class FreeSurfaceRemoveFalseInterfaceCells2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceRemoveFalseInterfaceCells2D<T, Descriptor>::
    processGenericBlocks(Box2D domain,
                         std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  using namespace freeSurfaceFlag2D;

  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  /// In the following, the flag status of cells is read (non-locally) and
  /// modified (locally). To avoid conflict, two loops are made, the first
  /// of which reads only, and the second writes. The vectors
  /// "interfaceToFluidNodes" and "interfaceToEmptyNodes" store coordinates of
  /// nodes that will switch status.
  std::vector<Node> interfaceToFluidNodes, interfaceToEmptyNodes;
  for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
    for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
      Node node(iX, iY);
      if (param.flag(iX, iY) == interface) {
        bool noEmptyNeighbor = true;
        bool noFluidNeighbor = true;

        for (plint iPop = 1; iPop < D::q; iPop++) {
          plint nextX = iX + D::c[iPop][0];
          plint nextY = iY + D::c[iPop][1];

          if (isEmpty(param.flag(nextX, nextY))) noEmptyNeighbor = false;
          if (isFullWet(param.flag(nextX, nextY))) noFluidNeighbor = false;
        }

        if (noEmptyNeighbor) {
          interfaceToFluidNodes.push_back(Node(iX, iY));
          // Store the coordinates, so flag on this node
          // can be changed in a loop outside the current one.

          T massExcess = param.mass(iX, iY) - param.getDensity(iX, iY);
          param.massExcess().insert(std::pair<Node, T>(node, massExcess));
          param.mass(iX, iY) = param.getDensity(iX, iY);
          param.volumeFraction(iX, iY) = T(1);
        } else if (noFluidNeighbor) {
          interfaceToEmptyNodes.push_back(Node(iX, iY));
          // Store the coordinates, so flag on this node
          // can be changed in a loop outside the current one.

          T massExcess = param.mass(iX, iY);
          param.massExcess().insert(std::pair<Node, T>(node, massExcess));

          param.attributeDynamics(iX, iY,
                                  new NoDynamics<T, Descriptor>(rhoDefault));
          param.mass(iX, iY) = T();
          param.setDensity(iX, iY, rhoDefault);
          param.volumeFraction(iX, iY) = T();
          param.setForce(iX, iY,
                         Array<T, D::ExternalField::sizeOfForce>::zero());
          param.setMomentum(iX, iY, Array<T, 2>(T(), T()));
        }
      }
    }
  }

  for (pluint i = 0; i < interfaceToFluidNodes.size(); ++i) {
    Node const &pos = interfaceToFluidNodes[i];
    param.flag(pos[0], pos[1]) = fluid;
  }
  for (pluint i = 0; i < interfaceToEmptyNodes.size(); ++i) {
    Node const &pos = interfaceToEmptyNodes[i];
    param.flag(pos[0], pos[1]) = empty;
  }
}

/* *************** Class FreeSurfaceEqualMassExcessReDistribution2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceEqualMassExcessReDistribution2D<T, Descriptor>::
    processGenericBlocks(Box2D domain,
                         std::vector<AtomicBlock2D *> atomicBlocks) {
  typedef Descriptor<T> D;
  typedef typename InterfaceLists2D<T, Descriptor>::Node Node;
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  Box2D originalDomain(domain);

  typename std::map<Node, T>::iterator iEle = param.massExcess().begin();
  for (; iEle != param.massExcess().end(); ++iEle) {
    Array<plint, 2> node = iEle->first;
    plint iX = node[0];
    plint iY = node[1];

    // Check for valid interface neighbors to re-distribute mass
    if (contained(iX, iY, domain.enlarge(1))) {
      std::vector<int> indX, indY;
      plint numValidNeighbors = 0;

      // Check for interface neighbors in the LB directions.
      for (plint iPop = 1; iPop < D::q; iPop++) {
        plint nextX = iX + D::c[iPop][0];
        plint nextY = iY + D::c[iPop][1];
        if (param.flag(nextX, nextY) == interface) {
          if (contained(nextX, nextY, domain)) {
            indX.push_back(nextX);
            indY.push_back(nextY);
          }
          numValidNeighbors++;
        }
      }

      // Mass re-distribution
      if (numValidNeighbors != 0) {
        int indSize = (int)indX.size();
        T massToRedistribute = iEle->second / (T)numValidNeighbors;

        for (int i = 0; i < indSize; i++) {
          int nextX = indX[i];
          int nextY = indY[i];

          param.mass(nextX, nextY) += massToRedistribute;
          param.volumeFraction(nextX, nextY) =
              param.mass(nextX, nextY) / param.getDensity(nextX, nextY);
        }
      } else {
        if (contained(iX, iY, originalDomain)) {
          param.addToLostMass(iEle->second);
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceComputeStatistics2D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeStatistics2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isWet(param.flag(iX, iY))) {
        param.addToTotalMass(param.mass(iX, iY));
        if (param.flag(iX, iY) == interface) {
          param.addToInterfaceCells(1);
        }
      }
    }
  }
}

/* *************** Class FreeSurfaceAddExternalForce2D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddExternalForce2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;
  FreeSurfaceProcessorParam2D<T, Descriptor> param(atomicBlocks);

  if (Descriptor<T>::ExternalField::sizeOfForce == 0) {
    return;
  }
  PLB_ASSERT(Descriptor<T>::ExternalField::sizeOfForce == 2);  // mehdi: or 3?

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (isWet(param.flag(iX, iY))) {
        Array<T, 2> newJ = param.getMomentum(iX, iY);
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force =
            param.getForce(iX, iY);

        T dynamicOmega = param.cell(iX, iY).getDynamics().getDynamicParameter(
            dynamicParams::dynamicOmega,
            param.cell(iX, iY));  // In case of a Smagorinsky model.
        T tau = 0.0;
        if (!util::isZero(dynamicOmega)) {
          tau = (T)1 / dynamicOmega;
        } else {
          tau = (T)1 / param.cell(iX, iY).getDynamics().getOmega();
        }

        // Two comments:
        // - Here the force is multiplied by rho0 and not rho so that, under
        //   gravity, a linear pressure profile is obtained.
        // - The force is not multiplied by the volume fraction (some authors
        //   do multiply it by the volumeFraction), because there is a
        //   point-wise interpretation of quantities like momentum.
        for (plint i = 0; i < Descriptor<T>::ExternalField::sizeOfForce; i++) {
          newJ[i] += rhoDefault * tau * force[i];
        }
        param.setMomentum(iX, iY, newJ);
      }
    }
  }
}

/* *************** Class RepelInterfaceFromImmersedWalls2D
 * ******************************** */

template <typename T, class VelFunction>
void RepelInterfaceFromImmersedWalls2D<T, VelFunction>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  PLB_PRECONDITION(atomicBlocks.size() == 4);
  ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(atomicBlocks[0]);
  TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(atomicBlocks[1]);
  ScalarField2D<int> *flag =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[2]);
  AtomicContainerBlock2D *container =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[3]);
  PLB_ASSERT(rhoBar);
  PLB_ASSERT(j);
  PLB_ASSERT(flag);
  PLB_ASSERT(container);

  Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
  Dot2D ofsF = computeRelativeDisplacement(*rhoBar, *flag);

  ImmersedWallData2D<T> *wallData =
      dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2>> const &vertices = wallData->vertices;
  std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
  PLB_ASSERT(vertices.size() == globalVertexIds.size());

  if (strongRepelling) {
    // Define temporary matrices to facilitate the computation of the
    // mean position and velocity of the immersed boundary vertices
    // in a vicinity of each interface node.
    plint nx = domain.getNx();
    plint ny = domain.getNy();
    TensorField2D<T, 2> sumVelocity(nx, ny);
    ScalarField2D<plint> numVertices(nx, ny);
    Dot2D ofs(-domain.x0, -domain.y0);

    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const &vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          if (contained(pos[0], pos[1], domain)) {
            if (flag->get(pos[0] + ofsF.x, pos[1] + ofsF.y) == interface) {
              sumVelocity.get(pos[0] + ofs.x, pos[1] + ofs.y) +=
                  velFunction(globalVertexIds[i]);
              numVertices.get(pos[0] + ofs.x, pos[1] + ofs.y) += 1;
            }
          }
        }
      }
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
        if (flag->get(iX + ofsF.x, iY + ofsF.y) == interface &&
            numVertices.get(iX + ofs.x, iY + ofs.y) != 0) {
          Array<T, 2> meanVelocityWall =
              sumVelocity.get(iX + ofs.x, iY + ofs.y) /
              (T)numVertices.get(iX + ofs.x, iY + ofs.y);

          rhoBar->get(iX, iY) = rhoDefault - (T)1.0;
          j->get(iX + ofsJ.x, iY + ofsJ.y) = rhoDefault * meanVelocityWall;
        }
      }
    }
  } else {
    // Define temporary matrices to facilitate the computation of the
    // mean position and velocity of the immersed boundary vertices
    // in a vicinity of each interface node.
    plint nx = domain.getNx();
    plint ny = domain.getNy();
    TensorField2D<T, 2> sumPosition(nx, ny);
    TensorField2D<T, 2> sumVelocity(nx, ny);
    ScalarField2D<plint> numVertices(nx, ny);
    Dot2D ofs(-domain.x0, -domain.y0);

    for (pluint i = 0; i < vertices.size(); ++i) {
      Array<T, 2> const &vertex = vertices[i];
      Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
      // x   x . x   x
      for (plint dx = -1; dx <= +2; ++dx) {
        for (plint dy = -1; dy <= +2; ++dy) {
          Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
          if (contained(pos[0], pos[1], domain)) {
            if (flag->get(pos[0] + ofsF.x, pos[1] + ofsF.y) == interface) {
              sumPosition.get(pos[0] + ofs.x, pos[1] + ofs.y) += vertex;
              sumVelocity.get(pos[0] + ofs.x, pos[1] + ofs.y) +=
                  velFunction(globalVertexIds[i]);
              numVertices.get(pos[0] + ofs.x, pos[1] + ofs.y) += 1;
            }
          }
        }
      }
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
        if (flag->get(iX + ofsF.x, iY + ofsF.y) == interface &&
            numVertices.get(iX + ofs.x, iY + ofs.y) != 0) {
          // For each interface node which is close to an immersed boundary,
          // compute the mean position vector of its neighbor empty nodes.
          // This way we can calculate on which side of the immersed wall
          // the specific bubble (or droplet) is situated.
          Array<T, 2> meanPositionEmpty((T)0, (T)0);
          plint numEmpty = 0;
          for (plint dx = -2; dx <= +2; ++dx) {
            plint x = iX + dx;
            for (plint dy = -2; dy <= +2; ++dy) {
              plint y = iY + dy;
              if (isEmpty(flag->get(x + ofsF.x, y + ofsF.y))) {
                meanPositionEmpty += Array<T, 2>((T)x, (T)y);
                numEmpty++;
              }
            }
          }
          PLB_ASSERT(numEmpty);
          meanPositionEmpty /= (T)numEmpty;

          // Decide if the immersed boundary and the bubble (or droplet) are
          // approaching one another. In such a case change the density and the
          // momentum of the bubble (or droplet) interface nodes.
          Array<T, 2> meanPositionWall =
              sumPosition.get(iX + ofs.x, iY + ofs.y) /
              (T)numVertices.get(iX + ofs.x, iY + ofs.y);
          Array<T, 2> meanVelocityWall =
              sumVelocity.get(iX + ofs.x, iY + ofs.y) /
              (T)numVertices.get(iX + ofs.x, iY + ofs.y);
          T dotProd =
              dot(meanVelocityWall, meanPositionEmpty - meanPositionWall);

          if (util::greaterEqual(dotProd, (T)0)) {
            rhoBar->get(iX, iY) = rhoDefault - (T)1.0;
            j->get(iX + ofsJ.x, iY + ofsJ.y) = rhoDefault * meanVelocityWall;
          }
        }
      }
    }
  }
}

/* *************** Class TemporarilyProtectImmersedWalls2D
 * ******************************** */

template <typename T>
void TemporarilyProtectImmersedWalls2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  PLB_PRECONDITION(atomicBlocks.size() == 2);
  ScalarField2D<int> *flag =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[0]);
  AtomicContainerBlock2D *container =
      dynamic_cast<AtomicContainerBlock2D *>(atomicBlocks[1]);
  PLB_ASSERT(flag);
  PLB_ASSERT(container);

  ImmersedWallData2D<T> *wallData =
      dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
  PLB_ASSERT(wallData);

  std::vector<Array<T, 2>> const &vertices = wallData->vertices;

  for (pluint i = 0; i < vertices.size(); ++i) {
    Array<T, 2> const &vertex = vertices[i];
    Array<plint, 2> intPos((plint)vertex[0], (plint)vertex[1]);
    // x   x . x   x
    for (plint dx = -1; dx <= +2; ++dx) {
      for (plint dy = -1; dy <= +2; ++dy) {
        Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
        if (contained(pos[0], pos[1], domain)) {
          if (flag->get(pos[0], pos[1]) == fluid) {
            flag->get(pos[0], pos[1]) = temporarilyProtect;
          }
        }
      }
    }
  }
}

/* *************** Class RemoveProtectionFromImmersedWalls2D
 * ******************************** */

template <typename T>
void RemoveProtectionFromImmersedWalls2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks) {
  using namespace freeSurfaceFlag2D;

  PLB_PRECONDITION(atomicBlocks.size() == 1);
  ScalarField2D<int> *flag =
      dynamic_cast<ScalarField2D<int> *>(atomicBlocks[0]);
  PLB_ASSERT(flag);

  for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
      if (flag->get(iX, iY) == temporarilyProtect) {
        flag->get(iX, iY) = fluid;
      }
    }
  }
}

}  // namespace lbfoam
}  // namespace plb

#endif  // FREE_SURFACE_MODEL_2D_HH
