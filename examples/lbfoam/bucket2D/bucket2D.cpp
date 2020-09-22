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

#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "fenv.h"
#include "palabos2D.h"
#include "palabos2D.hh"
using namespace plb;
using namespace lbfoam;

// Define type of lattice
#define DESCRIPTOR descriptors::ForcedD2Q9Descriptor
//#define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionWithSourceD2Q5Descriptor

// Define dynamics model (BGK, RLB, etc.)

#define ADYNAMICS AdvectionDiffusionBGKdynamics
#define ADYNAMICSWS AdvectionDiffusionWithSourceRLBdynamics
#define PADDING 8

typedef double T;

std::string outDir;

struct SimulationParameters {
  // Geometric parameters.
  std::map<int, Array<T, 2>> nucleiCenters;
  int numNuclei;
  int numRows;
  int numberOfBubbles;
  T shift;
  T radius;
  T bdy;
  T bdx;
  std::string distribution;
  plint packingOffset;

  T contactAngle;
  T gasIni_LB;
  T temperature;
  T p_ini;
  T R_s;
  plint bucketThickness_LB;
  plint bucketOffset_LB;
  plint bucketHeight_LB;
  plint fluidPoolHeight_LB;

  plint maxIter;
  T cSmago;
  bool freezeLargestBubble;
  bool surfaceDiffusion;
  bool gravity;
  bool entrapBubbles;

  T bubbleVolumeRatio;
  T alpha, beta;

  plint statIter;  // Output parameters.
  plint outIter;

  plint nx, ny;
  Array<T, 2> gVector_LB;
  T g_LB;
  T rho_LB;
  T tau_LB;
  T sigma_LB;
  T tauD_LB;
  T omega;
  T adOmega;
  T kh_LB;
  T pi_LB;
  T source_LB;
};

SimulationParameters param;

// This function reads the simulation parameters from the XML file
void readUserDefinedSimulationParameters(std::string xmlInputFileName,
                                         SimulationParameters &param) {
  XMLreader document(xmlInputFileName);

  document["geometry"]["simulationDomain"]["nx"].read(param.nx);
  document["geometry"]["simulationDomain"]["ny"].read(param.ny);
  document["geometry"]["fluidPoolHeight_LB"].read(param.fluidPoolHeight_LB);
  document["geometry"]["bucketOffset_LB"].read(param.bucketOffset_LB);
  document["geometry"]["bucketHeight_LB"].read(param.bucketHeight_LB);
  document["geometry"]["bucketThickness_LB"].read(param.bucketThickness_LB);

  document["Nucleation"]["numNuclei"].read(param.numNuclei);
  document["Nucleation"]["numRows"].read(param.numRows);
  document["Nucleation"]["radius"].read(param.radius);
  document["Nucleation"]["shift"].read(param.shift);
  document["Nucleation"]["distribution"].read(param.distribution);
  document["Nucleation"]["packingOffset"].read(param.packingOffset);
  document["Nucleation"]["numberOfBubbles"].read(param.numberOfBubbles);

  document["fluid"]["rho_LB"].read(param.rho_LB);
  document["fluid"]["R_s"].read(param.R_s);
  document["fluid"]["p_ini"].read(param.p_ini);
  document["fluid"]["temperature"].read(param.temperature);
  document["fluid"]["tau_LB"].read(param.tau_LB);
  document["fluid"]["tauD_LB"].read(param.tauD_LB);
  document["fluid"]["sigma_LB"].read(param.sigma_LB);
  document["fluid"]["kh_LB"].read(param.kh_LB);
  document["fluid"]["gasIni_LB"].read(param.gasIni_LB);
  document["fluid"]["pi_LB"].read(param.pi_LB);
  document["fluid"]["contactAngle"].read(param.contactAngle);
  document["fluid"]["surfaceDiffusion"].read(param.surfaceDiffusion);
  document["fluid"]["source_LB"].read(param.source_LB);

  T pi = acos((T)-1);
  param.contactAngle *= pi / 180.0;

  document["numerics"]["maxIter"].read(param.maxIter);
  document["numerics"]["cSmago"].read(param.cSmago);
  document["numerics"]["freezeLargestBubble"].read(param.freezeLargestBubble);
  document["numerics"]["gravity"].read(param.gravity);
  document["numerics"]["g_LB"].read(param.g_LB);
  document["numerics"]["entrapBubbles"].read(param.entrapBubbles);

  document["numerics"]["bubbleVolumeRatio"].read(param.bubbleVolumeRatio);
  document["numerics"]["alpha"].read(param.alpha);
  document["numerics"]["beta"].read(param.beta);

  document["output"]["statIter"].read(param.statIter);
  document["output"]["outIter"].read(param.outIter);
  document["output"]["outDir"].read(outDir);
}

// This class is used to implement the gas source term
template <typename T, template <typename U> class Descriptor>
class SourceTerm : public BoxProcessingFunctional2D_L<T, Descriptor> {
 public:
  SourceTerm(T source_) : source(source_){};
  virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice) {
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
      for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
        lattice.get(iX, iY).setExternalField(
            Descriptor<T>::ExternalField::scalarBeginsAt,
            Descriptor<T>::ExternalField::sizeOfScalar, &source);
      }
    }
  };
  virtual SourceTerm<T, Descriptor> *clone() const {
    return new SourceTerm<T, Descriptor>(*this);
  };
  virtual void getTypeOfModification(
      std::vector<modif::ModifT> &modified) const {
    modified[0] = modif::staticVariables;
  };

 private:
  T source;
};

void calculateDerivedSimulationParameters(SimulationParameters &param) {
  // Derived quantities.

  plint totBubinTwoRows = 2 * param.numNuclei - 1;

  if (!param.gravity) {
    param.g_LB = 0.;
  }

  param.gVector_LB = Array<T, 2>((T)0, (T)-param.g_LB);
  param.adOmega = 1.0 / param.tauD_LB;
  param.omega = 1.0 / param.tau_LB;

  // Calculate the offsets for bubble nucleation

  plint totalOffset =
      param.shift + param.bucketThickness_LB + param.bucketOffset_LB;
  param.bdx = (param.nx - param.bucketThickness_LB - param.bucketOffset_LB -
               2 * (param.shift)) /
              (totBubinTwoRows - 1);
  param.bdy = param.fluidPoolHeight_LB / (param.numRows + 1);

  // Distribution of bubbles
  if (param.distribution == "organized") {
    int N = 0;

    for (int j = 0; j < param.numRows; j++) {
      for (int i = 0; i <= totBubinTwoRows; i++) {
        T X = 0;
        T Y = 0;

        if ((i + j) % 2 == 0) {
          if (j % 2 == 1 && i == totBubinTwoRows) {
            continue;
          }

          X = i * param.bdx + totalOffset;
          Y = (j + 1) * param.bdy;

          Array<T, 2> center(X, Y);

          param.nucleiCenters.insert(std::pair<int, Array<T, 2>>(N, center));

          N++;
        } else {
          continue;
        }
      }
    }
  } else if (param.distribution == "random") {
    T lowX = totalOffset;
    T highX = param.nx - param.shift - param.bucketOffset_LB -
              param.bucketThickness_LB;

    T lowY = param.shift;
    T highY = param.fluidPoolHeight_LB - param.shift;

    namespace pds = thinks::poisson_disk_sampling;

    // Input parameters.

    T diameter = 2 * param.radius + param.packingOffset;
    const auto min = std::array<T, 2>{{lowX, lowY}};
    const auto max = std::array<T, 2>{{highX, highY}};

    // Samples returned as std::vector<std::array<float, 2>>.
    // Default seed and max attempts.
    auto samples = pds::PoissonDiskSampling<T>(diameter, min, max);

    auto rng = std::default_random_engine(param.numberOfBubbles);

    std::shuffle(samples.begin(), samples.end(), rng);

    int N = 0;
    for (auto sample : samples) {
      if (N == param.numberOfBubbles) {
        break;
      }

      // int yDist = std::abs(sample[1] - param.fluidPoolHeight_LB / 2);

      // int s = rand() % (param.fluidPoolHeight_LB / 2);

      Array<T, 2> center(sample[0], sample[1]);

      // Store the location of each bubble

      param.nucleiCenters.insert(std::pair<int, Array<T, 2>>(N, center));

      N++;
    }
  }

  else if (param.distribution == "one") {
    T X = param.nx / 2;
    T Y = param.fluidPoolHeight_LB / 2;

    Array<T, 2> center(X, Y);

    param.nucleiCenters.insert(std::pair<int, Array<T, 2>>(0, center));
  }
}

void printSimulationParameters(SimulationParameters const &param) {
  pcout << "fluidPoolHeight_LB = " << param.fluidPoolHeight_LB << std::endl;

  pcout << "g_LB = (" << param.gVector_LB[0] << ", " << param.gVector_LB[1]
        << " )" << std::endl;
  pcout << "rho_LB = " << param.rho_LB << std::endl;
  pcout << "sigma_LB = " << param.sigma_LB << std::endl;
  pcout << "omega = " << param.omega << std::endl;
  pcout << "tau_LB = " << param.tau_LB << std::endl;
  pcout << std::endl;
  pcout << "contractAngle = " << param.contactAngle * 180.0 / acos((T)-1)
        << std::endl;

  pcout << "gas_ini = " << param.gasIni_LB << std::endl;
  pcout << "tau_LB = " << param.tau_LB << std::endl;
  pcout << "adOmega = " << param.adOmega << std::endl;
  pcout << "Kh_LB = " << param.kh_LB << std::endl;
  pcout << "maxIter = " << param.maxIter << std::endl;
  pcout << "cSmago = " << param.cSmago << std::endl;
  pcout << "freezeLargestBubble = "
        << (param.freezeLargestBubble ? "true" : "false") << std::endl;
  pcout << "bubbleVolumeRatio = " << param.bubbleVolumeRatio << std::endl;
  pcout << "alpha = " << param.alpha << std::endl;
  pcout << "beta = " << param.beta << std::endl;

  pcout << "statIter = " << param.statIter << std::endl;
  pcout << "outIter = " << param.outIter << std::endl;

  pcout << "nx = " << param.nx << std::endl;
  pcout << "ny = " << param.ny << std::endl;
  pcout << "distribution = " << param.distribution << std::endl;

  pcout << "bucket offset = " << param.bucketOffset_LB << std::endl;
  pcout << "bucket thickness = " << param.bucketThickness_LB << std::endl;
  pcout << "bucket height = " << param.bucketHeight_LB << std::endl;
}

bool insideSphere(T x, T y) {
  bool isInside = false;

  Array<T, 2> pos(x, y);
  typename std::map<int, Array<T, 2>>::const_iterator it =
      param.nucleiCenters.begin();
  for (; it != param.nucleiCenters.end(); ++it) {
    Array<T, 2> center = it->second;
    T r = norm<T, 2>(pos - center);
    if (r <= param.radius) {
      isInside = true;
    }
  }

  return isInside;
}

void nucleateBubbles(FreeSurfaceFields2D<T, DESCRIPTOR> &fields,
                     Dynamics<T, DESCRIPTOR> *dynamics) {
  typename std::map<int, Array<T, 2>>::const_iterator it =
      param.nucleiCenters.begin();

  for (; it != param.nucleiCenters.end(); ++it) {
    Array<T, 2> center = it->second;

    punchSphere(fields, center, param.radius, param.rho_LB, param.rho_LB,
                *dynamics);
  }
}

bool insideFluidPool(T x, T y) {
  if (y <= param.fluidPoolHeight_LB &&
      x < param.nx - (param.bucketOffset_LB + param.bucketThickness_LB) &&
      x > param.bucketOffset_LB + param.bucketThickness_LB) {
    return true;
  }
  return false;
}

bool insideFluid(T x, T y) {
  if (insideFluidPool(x, y)) {
    return true;
  }
  return false;
}

int initialFluidFlags(plint iX, plint iY) {
  bool insideBucket =
      ((iX > param.bucketOffset_LB &&
        iX <= param.bucketOffset_LB + param.bucketThickness_LB) ||
       (iX < param.nx - param.bucketOffset_LB &&
        iX >= param.nx - param.bucketOffset_LB - param.bucketThickness_LB)) &&
      (iY < param.bucketHeight_LB);

  if (insideFluid(iX, iY)) {
    return freeSurfaceFlag2D::fluid;
  } else if (insideBucket) {
    return freeSurfaceFlag2D::wall;
  }
  return freeSurfaceFlag2D::empty;
}

void writeResults(FreeSurfaceFields2D<T, DESCRIPTOR> *fields,
                  MultiBlockLattice2D<T, ADESCRIPTOR> adLattice,
                  MultiScalarField2D<plint> *tagMatrix,
                  MultiScalarField2D<double> *disjoiningPressureField,
                  plint iT) {
  std::vector<T> isoLevels;
  isoLevels.push_back(0.5);

  VtkImageOutput2D<T> vtkOut(
      createFileName(outDir + "volumeData_", iT, PADDING));
  std::auto_ptr<MultiTensorField2D<T, 2>> v = computeVelocity(fields->lattice);
  std::auto_ptr<MultiScalarField2D<T>> rho = computeDensity(fields->lattice);
  std::auto_ptr<MultiScalarField2D<T>> adrho = computeDensity(adLattice);
  std::auto_ptr<MultiScalarField2D<T>> smoothVF(lbmSmoothen<T, DESCRIPTOR>(
      fields->volumeFraction, fields->volumeFraction.getBoundingBox()));
  vtkOut.writeData<2, float>(*v, "velocity");
  vtkOut.writeData<float>(*rho, "pressure", 1);
  vtkOut.writeData<float>(*adrho, "adDensity", 1);
  vtkOut.writeData<float>(fields->volumeFraction, "volumeFraction", 1.0);
  vtkOut.writeData<float>(*smoothVF, "smoothedVolumeFraction", 1.0);
  // vtkOut.writeData<float>(fields->outsideDensity, "outsidePressure",
  //        param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
  vtkOut.writeData<float>(*copyConvert<plint, T>(*tagMatrix), "bubbleTags",
                          1.0);
  vtkOut.writeData<float>(*copyConvert<double, T>(*disjoiningPressureField),
                          "disjoiningPressure", 1.0);

  ImageWriter<T> imageWriter("leeloo");

  imageWriter.writeScaledGif(createFileName(outDir + "smoothVOF", iT, 6),
                             *smoothVF);

  imageWriter.writeScaledGif(createFileName(outDir + "gas", iT, 6), *adrho);
}

plint numRememberedVolumes = 1;

int main(int argc, char **argv) {
  plbInit(&argc, &argv);

  std::cout.precision(10);
  std::scientific(std::cout);

  // Command-line arguments

  if (argc != 2) {
    pcout << "Usage: " << argv[0] << " xml-input-file-name" << std::endl;
    exit(1);
  }

  std::string xmlInputFileName;
  xmlInputFileName = std::string(argv[1]);

  // Set the simulation parameters.

  readUserDefinedSimulationParameters(xmlInputFileName, param);
  calculateDerivedSimulationParameters(param);
  printSimulationParameters(param);

  SparseBlockStructure2D blockStructure(
      createRegularDistribution2D(param.nx, param.ny));

  //    Dynamics<T, DESCRIPTOR>* dynamics =
  //            new SmagorinskyBGKdynamics<T, DESCRIPTOR>(param.omega,
  //            param.cSmago);

  Dynamics<T, DESCRIPTOR> *dynamics =
      new BGKdynamics<T, DESCRIPTOR>(param.omega);

  Dynamics<T, ADESCRIPTOR> *adynamics =
      new ADYNAMICSWS<T, ADESCRIPTOR>(param.adOmega);
  Dynamics<T, ADESCRIPTOR> *emptyDynamics = new NoDynamics<T, ADESCRIPTOR>();

  FreeSurfaceFields2D<T, DESCRIPTOR> fields(
      blockStructure, dynamics->clone(), param.rho_LB, param.sigma_LB,
      param.contactAngle, param.gVector_LB);
  MultiBlockLattice2D<T, ADESCRIPTOR> adLattice(param.nx, param.ny,
                                                adynamics->clone());
  Array<T, 2> u0((T)0, (T)0);

  // Initialization
  initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(),
                          param.gasIni_LB, u0);
  adLattice.initialize();

  pcout << "Setting up initial condition." << std::endl;

  Box2D bottom(0, param.nx - 1, 0, 0);
  Box2D top(0, param.nx - 1, param.ny - 1, param.ny - 1);
  Box2D lateral1(0, 0, 0, param.ny - 1);
  Box2D lateral2(param.nx - 1, param.nx - 1, 0, param.ny - 1);
  // Set the walls boundary conditions
  setToConstant(fields.flag, bottom, (int)freeSurfaceFlag2D::wall);
  setToConstant(fields.flag, top, (int)freeSurfaceFlag2D::wall);
  setToConstant(fields.flag, lateral1, (int)freeSurfaceFlag2D::wall);
  setToConstant(fields.flag, lateral2, (int)freeSurfaceFlag2D::wall);
  // Flag the liquid cells
  setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1),
                initialFluidFlags);
  // Nucleate bubbles
  nucleateBubbles(fields, dynamics->clone());
  // analyticalIniVolumeFraction(fields.volumeFraction, fields.flag,
  // insideFluid, 							32);

  fields.periodicityToggleAll(false);
  adLattice.periodicity().toggleAll(false);
  fields.partiallyDefaultInitialize();

  plint iniIter = 0;

  BubbleTracking2D bubbleTracking(fields.flag);
  BubbleGrowth2D<T> bubbleGrowth(fields.flag);

  std::string fname = outDir + "bubbles.log";
  FILE *fp = fopen(fname.c_str(), "w");

  pcout << std::endl;

  // MultiScalarField2D<T> newvof = fields.volumeFraction;
  MultiScalarField2D<T> oldvof = fields.volumeFraction;
  // Integrate the source term
  integrateProcessingFunctional(new SourceTerm<T, ADESCRIPTOR>(param.source_LB),
                                adLattice.getBoundingBox(), adLattice);

  // Main iteration loop.
  for (plint iT = iniIter; iT < param.maxIter; iT++) {
    if (iT % param.statIter == 0 || iT == param.maxIter - 1) {
      pcout << "At iteration " << iT << std::endl;
      T avE = computeAverageEnergy(fields.lattice);
      pcout << "Average kinetic energy: " << avE << std::endl;
      plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
      pcout << "Number of interface cells: " << numIntCells << std::endl;
      if (iT != iniIter) {
        pcout << "Time spent for each iteration: "
              << global::timer("iteration").getTime() / (T)param.statIter
              << std::endl;
        global::timer("iteration").reset();
      }
      pcout << std::endl;
    }

    if ((iT % param.outIter == 0 && iT != 0) || iT == 1 ||
        iT == param.maxIter - 1) {
      pcout << "Writing results at iteration " << iT << std::endl;
      global::timer("images").restart();
      writeResults(&fields, adLattice, bubbleTracking.getTagMatrix(),
                   bubbleTracking.getDisjoiningPressureField(), iT);
      global::timer("images").stop();
      pcout << "Time spent for writing results: "
            << global::timer("images").getTime() << std::endl;
      pcout << std::endl;
    }

    global::timer("iteration").start();
    T bubbleVolumeRatio = iT == 0 ? 1.0 : param.bubbleVolumeRatio;
    bubbleTracking.execute<T, DESCRIPTOR, ADESCRIPTOR>(
        fields.volumeFraction, fields.flag, fields.normal, fields.rhoBar,
        fields.mass, fields.j, adLattice, oldvof, false, param.pi_LB);

    oldvof = fields.volumeFraction;

    // Calculate bubble transitions and the disjoining pressure
    bubbleGrowth.transition(bubbleTracking, iT, param.temperature, param.R_s,
                            param.p_ini, 1., param.rho_LB, bubbleVolumeRatio,
                            param.entrapBubbles, numRememberedVolumes);

    // Update bubble gas content
    bubbleGrowth.updateBubbleGrowth(fields.outsideDensity, param.rho_LB,
                                    param.alpha, param.beta, (T)1.);

    // Free surface
    fields.lattice.executeInternalProcessors();
    fields.lattice.evaluateStatistics();
    fields.lattice.incrementTime();

    if (iT == 0 && param.freezeLargestBubble) {
      bubbleGrowth.freezeLargestBubble();
    }

    // order matters
    std::vector<MultiBlock2D *> couplingBlocks;
    couplingBlocks.push_back(&adLattice);
    couplingBlocks.push_back(&fields.lattice);
    couplingBlocks.push_back(&fields.flag);
    couplingBlocks.push_back(&fields.j);
    couplingBlocks.push_back(&fields.outsideDensity);
    couplingBlocks.push_back(bubbleGrowth.getOldTagMatrix());

    adLattice.collideAndStream();
    // Calculate the gas diffused into each bubble
    applyProcessingFunctional(
        new GrowthCoupling2D<T, ADESCRIPTOR, DESCRIPTOR>(
            adynamics->clone(), emptyDynamics->clone(), param.kh_LB,
            bubbleGrowth.getBubbles(), param.surfaceDiffusion),
        adLattice.getBoundingBox(), couplingBlocks);


    global::timer("iteration").stop();

    if (iT % param.statIter == 0 || iT == param.maxIter - 1) {
      bubbleGrowth.timeHistoryLog(outDir + "bubbleTimeHistory.log");
      bubbleGrowth.fullBubbleLog(outDir + "FullBubbleRecord2D.log");

      // We do not log frozen bubbles.
      std::map<plint, BubbleInfo2D>::const_iterator it =
          bubbleGrowth.getBubbles().begin();
      T totalBubbleVolume = T();
      T currentDensity = T();
      T totalDisjoining = T();
      plint numBubbles = 0;
      for (; it != bubbleGrowth.getBubbles().end(); ++it) {
        if (it->second.isFrozen()) {
          pcout << "Bubble with this ID is frozen: " << it->first << std::endl;
          continue;
        }
        numBubbles++;
        T v = it->second.getVolume();
        T d = it->second.getCurrentDensity();
        T j = it->second.getDisjoiningPressure();

        totalBubbleVolume += v;
        currentDensity += d;
        totalDisjoining += j;
      }
      pcout << "At iteration " << iT << ", the number of bubbles is "
            << numBubbles << std::endl;
      pcout << "The total volume of bubbles is: " << totalBubbleVolume
            << std::endl;
      pcout << "The total density is: " << currentDensity << std::endl;
      pcout << "The total disjoining pressure is: " << totalDisjoining
            << std::endl;

      pcout << std::endl;
      fflush(fp);
    }
  }

  fclose(fp);
  delete dynamics;
  delete emptyDynamics;
  delete adynamics;

  exit(0);
}
