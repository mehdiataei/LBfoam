/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* The breaking dam free surface problem. This code demonstrates the basic usage
 * of the free surface module in Palabos. Surface tension and contact angles are
 * optional.
 */

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace lbfoam;

#define DESCRIPTOR descriptors::ForcedD2Q9Descriptor

typedef double T;

// Smagorinsky constant for LES model.
const T cSmago = 0.14;

// Physical dimensions of the system (in meters).
const T lx = 4.0;
const T ly = 1.0;

const T rhoEmpty = T(1);

plint writeImagesIter = 100;
plint getStatisticsIter = 200;

plint maxIter;
plint N;
plint nx, ny;
T delta_t, delta_x;
Array<T, 2> externalForce;
T nuPhys, nuLB, tau, omega, Bo, surfaceTensionLB, contactAngle;

std::string outDir;
plint obstacleCenterXYplane, obstacleLength, obstacleWidth, obstacleHeight,
    beginWaterReservoir, waterReservoirHeight;
plint waterLevelOne, waterLevelTwo, waterLevelThree, waterLevelFour;

void setupParameters() {
  delta_x = ly / N;
  nx = util::roundToInt(lx / delta_x);
  ny = util::roundToInt(ly / delta_x);

  // Gravity in lattice units.
  T gLB = 9.8 * delta_t * delta_t / delta_x;
  externalForce = Array<T, 2>(0., -gLB);
  tau = (nuPhys * DESCRIPTOR<T>::invCs2 * delta_t) / (delta_x * delta_x) + 0.5;
  omega = 1. / tau;
  nuLB = (tau - 0.5) * DESCRIPTOR<T>::cs2;  // Viscosity in lattice units.

  surfaceTensionLB = rhoEmpty * gLB * delta_x * delta_x / Bo;

  obstacleCenterXYplane = util::roundToInt(1.1 * N);
  obstacleWidth = util::roundToInt(0.3 * N);
  obstacleHeight = util::roundToInt(0.3 * N);
  beginWaterReservoir = util::roundToInt((0.744 + 1.248) * N);
  waterReservoirHeight = util::roundToInt(0.7 * N);

  waterLevelOne = util::roundToInt(0.496 * N);
  waterLevelTwo = util::roundToInt(2. * 0.496 * N);
  waterLevelThree = util::roundToInt(3. * 0.496 * N);
  waterLevelFour = util::roundToInt((3. * 0.496 + 1.150) * N);
}

// Specifies the initial condition for the fluid (each cell is assigned the
// flag "fluid", "empty", or "wall").
int initialFluidFlags(plint iX, plint iY) {
  // Place an obstacle on the left end, which is hit by the fluid.
  bool insideObstacle = iX >= obstacleCenterXYplane - obstacleWidth / 2 &&
                        iX <= obstacleCenterXYplane + obstacleWidth / 2 &&
                        iY <= obstacleHeight + 1;

  if (insideObstacle) {
    return freeSurfaceFlag2D::wall;
  } else if (iX >= beginWaterReservoir && iY <= waterReservoirHeight) {
    return freeSurfaceFlag2D::fluid;
  } else {
    return freeSurfaceFlag2D::empty;
  }
}

void writeResults(MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
                  MultiScalarField2D<T>& volumeFraction, plint iT) {
  static const plint nx = lattice.getNx();
  static const plint ny = lattice.getNy();
  Box2D slice(0, nx - 1, 0, ny - 1);
  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeScaledGif(createFileName("u", iT, 6),
                             *computeVelocityNorm(lattice, slice));

  imageWriter.writeScaledGif(createFileName("rho", iT, 6),
                             *computeDensity(lattice, slice));

  imageWriter.writeScaledGif(createFileName("volumeFraction", iT, 6),
                             *extractSubDomain(volumeFraction, slice));

  VtkImageOutput2D<T> vtkOut(createFileName("volumeFraction", iT, 6), 1.);
  vtkOut.writeData<float>(volumeFraction, "vf", 1.);
}

void writeStatistics(FreeSurfaceFields2D<T, DESCRIPTOR>& fields) {
  pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- "
        << std::endl;
  T averageMass = freeSurfaceAverageMass<T, DESCRIPTOR>(
      fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
  pcout << "Average Mass: " << averageMass << std::endl;
  T averageDensity = freeSurfaceAverageDensity<T, DESCRIPTOR>(
      fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
  pcout << "Average Density: " << std::setprecision(12) << averageDensity
        << std::endl;

  T averageVolumeFraction = freeSurfaceAverageVolumeFraction<T, DESCRIPTOR>(
      fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
  pcout << "Average Volume-Fraction: " << std::setprecision(12)
        << averageVolumeFraction << std::endl;

  pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- "
        << std::endl;
}

int main(int argc, char** argv) {
  plbInit(&argc, &argv);
  global::directories().setInputDir("./");

  if (global::argc() != 8) {
    pcout << "Error missing some input parameter\n";
  }

  try {
    global::argv(1).read(outDir);
    global::directories().setOutputDir(outDir + "/");
    global::argv(2).read(nuPhys);
    global::argv(3).read(Bo);
    global::argv(4).read(contactAngle);
    global::argv(5).read(N);
    global::argv(6).read(delta_t);
    global::argv(7).read(maxIter);
  } catch (PlbIOException& except) {
    pcout << except.what() << std::endl;
    pcout << "The parameters for this program are :\n";
    pcout << "1. Output directory name.\n";
    pcout << "2. kinematic viscosity in physical Units (m^2/s) .\n";
    pcout << "3. Bond number (Bo = rho * g * L^2 / gamma).\n";
    pcout << "4. Contact angle (in degrees) not working now so set -1.\n";
    pcout << "5. number of lattice nodes for lz .\n";
    pcout << "6. delta_t .\n";
    pcout << "7. maxIter .\n";
    pcout << "Reasonable parameters on a desktop computer are: "
          << (std::string)global::argv(0)
          << " tmp 1.e-6 100 -1 80 1.e-4 80000\n\n";
    pcout << "Reasonable parameters on a parallel machine are: "
          << (std::string)global::argv(0)
          << " tmp 1.e-6 100 -1 250 1.e-4 80000\n";
    exit(EXIT_FAILURE);
  }

  setupParameters();

  pcout << "delta_t= " << delta_t << std::endl;
  pcout << "delta_x= " << delta_x << std::endl;
  pcout << "delta_t*delta_t/delta_x= " << delta_t * delta_t / delta_x
        << std::endl;
  pcout << "externalForce= " << externalForce[1] << std::endl;
  pcout << "relaxation time= " << tau << std::endl;
  pcout << "omega= " << omega << std::endl;
  pcout << "kinematic viscosity physical units = " << nuPhys << std::endl;
  pcout << "kinematic viscosity lattice units= " << nuLB << std::endl;

  global::timer("initialization").start();

  SparseBlockStructure2D blockStructure(createRegularDistribution2D(nx, ny));

  Dynamics<T, DESCRIPTOR>* dynamics =
      new SmagorinskyBGKdynamics<T, DESCRIPTOR>(omega, cSmago);

  // If surfaceTensionLB is 0, then the surface tension algorithm is
  // deactivated. If contactAngle is less than 0, then the contact angle
  // algorithm is deactivated.
  FreeSurfaceFields2D<T, DESCRIPTOR> fields(blockStructure, dynamics->clone(),
                                            rhoEmpty, surfaceTensionLB,
                                            contactAngle, externalForce);
  // integrateProcessingFunctional(new ShortenBounceBack3D<T,DESCRIPTOR>,
  // fields.lattice.getBoundingBox(), fields.freeSurfaceArgs, 0);

  // Set all outer-wall cells to "wall" (here, bulk-cells are also set to
  // "wall", but it doesn't matter, because they are overwritten on the next
  // line).
  setToConstant(fields.flag, fields.flag.getBoundingBox(),
                (int)freeSurfaceFlag2D::wall);
  // In the bulk (all except outer wall layer), initialize the flags as
  // specified by the function "initialFluidFlags".
  setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1),
                initialFluidFlags);

  fields.defaultInitialize();

  pcout << "Time spent for setting up lattices: "
        << global::timer("initialization").stop() << std::endl;
  T lastIterationTime = T();

  for (plint iT = 0; iT <= maxIter; ++iT) {
    global::timer("iteration").restart();

    T sum_of_mass_matrix = T();
    T lost_mass = T();
    if (iT % getStatisticsIter == 0) {
      pcout << std::endl;
      pcout << "ITERATION = " << iT << std::endl;
      pcout << "Time of last iteration is " << lastIterationTime << " seconds"
            << std::endl;
      writeStatistics(fields);
      sum_of_mass_matrix = fields.lattice.getInternalStatistics().getSum(0);
      pcout << "Sum of mass matrix: " << sum_of_mass_matrix << std::endl;
      lost_mass = fields.lattice.getInternalStatistics().getSum(1);
      pcout << "Lost mass: " << lost_mass << std::endl;
      pcout << "Total mass: " << sum_of_mass_matrix + lost_mass << std::endl;
      pcout << "Interface cells: "
            << fields.lattice.getInternalStatistics().getIntSum(0) << std::endl;
    }

    if (iT % writeImagesIter == 0) {
      global::timer("images").start();
      writeResults(fields.lattice, fields.volumeFraction, iT);
      pcout << "Total time spent for writing images: "
            << global::timer("images").stop() << std::endl;
    }

    // This includes the collision-streaming cycle, plus all free-surface
    // operations.
    fields.lattice.executeInternalProcessors();
    fields.lattice.evaluateStatistics();
    fields.lattice.incrementTime();

    lastIterationTime = global::timer("iteration").stop();
  }
}
