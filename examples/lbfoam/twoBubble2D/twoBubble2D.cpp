/* This file is part of LBfoam library, a modified version of the Palabos library
 * Modified file: Copyright (C) 2020 Mehdi Ataei
 * E-mail contact: ataei@mie.utoronto.ca
 * 
 * Original file: Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The program is free software: you can redistribute it and/or
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


#include "palabos2D.h"
#include "palabos2D.hh"
#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>
#include "fenv.h"
#include <random>
using namespace plb;
using namespace lbfoam;

#define DESCRIPTOR descriptors::ForcedD2Q9Descriptor
//#define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionWithSourceD2Q5Descriptor

#define ADYNAMICS AdvectionDiffusionBGKdynamics
#define ADYNAMICSWS AdvectionDiffusionWithSourceRLBdynamics
#define PADDING 8

typedef double T;

std::string outDir;

struct SimulationParameters
{
	/*
	* Parameters set by the user.
	*/

	// Geometric parameters.
	std::map<int, Array<T, 2>> nucleiCenters;

	Array<T, 2> sphereCenter1_LB;
	Array<T, 2> sphereCenter2_LB;
	T sphereRadius1_LB;
	T sphereRadius2_LB;

	T gasIni_LB;
	T temperature;
	T p_ini;
	T R_s;
	plint bucketThickness_LB;
	plint bucketOffset_LB;
	plint bucketHeight_LB;

	plint maxIter;
	T cSmago;
	bool freezeLargestBubble;
	bool surfaceDiffusion;
	bool gravity;
	bool entrapBubbles;

	T bubbleVolumeRatio;
	T alpha, beta;

	plint statIter; // Output parameters.
	plint outIter;

	/*
	* Parameters NOT set by the user.
	*/

	plint nx, ny;
	T fluidPoolHeight_LB;
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

void readUserDefinedSimulationParameters(std::string xmlInputFileName,
										 SimulationParameters &param)
{
	XMLreader document(xmlInputFileName);

	document["geometry"]["simulationDomain"]["nx"].read(param.nx);
	document["geometry"]["simulationDomain"]["ny"].read(param.ny);
	document["geometry"]["fluidPoolHeight_LB"].read(param.fluidPoolHeight_LB);
	document["geometry"]["bucketOffset_LB"].read(param.bucketOffset_LB);
	document["geometry"]["bucketHeight_LB"].read(param.bucketHeight_LB);
	document["geometry"]["bucketThickness_LB"].read(param.bucketThickness_LB);

	document["nucleation"]["sphere1"]["radius1_LB"].read(param.sphereRadius1_LB);
	document["nucleation"]["sphere2"]["radius2_LB"].read(param.sphereRadius2_LB);

	std::vector<T> sphereCenter1;
	document["nucleation"]["sphere1"]["center1_LB"].read(sphereCenter1);
	PLB_ASSERT(sphereCenter1.size() == 2);
	param.sphereCenter1_LB[0] = sphereCenter1[0];
	param.sphereCenter1_LB[1] = sphereCenter1[1];

	std::vector<T> sphereCenter2;
	document["nucleation"]["sphere2"]["center2_LB"].read(sphereCenter2);
	PLB_ASSERT(sphereCenter2.size() == 2);
	param.sphereCenter2_LB[0] = sphereCenter2[0];
	param.sphereCenter2_LB[1] = sphereCenter2[1];

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
	document["fluid"]["surfaceDiffusion"].read(param.surfaceDiffusion);
	document["fluid"]["source_LB"].read(param.source_LB);

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

void calculateDerivedSimulationParameters(SimulationParameters &param)
{
	// Derived quantities.

	if (!param.gravity)
	{
		param.g_LB = 0.;
	}

	param.gVector_LB = Array<T, 2>((T)0, (T)-param.g_LB);
	param.adOmega = 1.0 / param.tauD_LB;
	param.omega = 1.0 / param.tau_LB;
}

void printSimulationParameters(SimulationParameters const &param)
{

	pcout << "fluidPoolHeight_LB = " << param.fluidPoolHeight_LB << std::endl;

	pcout << "g_LB = (" << param.gVector_LB[0] << ", " << param.gVector_LB[1] << " )"
		  << std::endl;
	pcout << "rho_LB = " << param.rho_LB << std::endl;
	pcout << "sigma_LB = " << param.sigma_LB << std::endl;
	pcout << "omega = " << param.omega << std::endl;
	pcout << "tau_LB = " << param.tau_LB << std::endl;
	pcout << std::endl;

	pcout << "gas_ini = " << param.gasIni_LB << std::endl;
	pcout << "tau_LB = " << param.tau_LB << std::endl;
	pcout << "adOmega = " << param.adOmega << std::endl;
	pcout << "Kh_LB = " << param.kh_LB << std::endl;
	pcout << "maxIter = " << param.maxIter << std::endl;
	pcout << "cSmago = " << param.cSmago << std::endl;
	pcout << "freezeLargestBubble = "
		  << (param.freezeLargestBubble ? "true" : "false") << std::endl;
}

bool insideSphere(T x, T y, Array<T, 2> sphereCenter_LB, T sphereRadius_LB)
{
	Array<T, 2> pos(x, y);
	T r = norm<T, 2>(pos - sphereCenter_LB);
	if (r <= sphereRadius_LB)
	{
		return true;
	}
	return false;
}

bool insideFluidPool(T x, T y)
{
	if (y <= param.fluidPoolHeight_LB &&
		x < param.nx - (param.bucketOffset_LB + param.bucketThickness_LB) &&
		x > param.bucketOffset_LB + param.bucketThickness_LB)
	{
		return true;
	}
	return false;
}

bool insideFluid(T x, T y)
{
	if (insideFluidPool(x, y) &&
		!insideSphere(x, y, param.sphereCenter1_LB, param.sphereRadius1_LB) &&
		!insideSphere(x, y, param.sphereCenter2_LB, param.sphereRadius2_LB))
	{
		return true;
	}
	return false;
}

// Specifies the initial condition for the fluid (each cell is assigned the
// flag "fluid", "empty", or "wall").
int initialFluidFlags(plint iX, plint iY)
{

	bool insideBucket = ((iX > param.bucketOffset_LB && iX <= param.bucketOffset_LB + param.bucketThickness_LB) ||
						 (iX < param.nx - param.bucketOffset_LB && iX >= param.nx - param.bucketOffset_LB - param.bucketThickness_LB)) &&
						(iY < param.bucketHeight_LB);

	if (insideFluid(iX, iY))
	{
		return freeSurfaceFlag2D::fluid;
	}
	else if (insideBucket)
	{
		return freeSurfaceFlag2D::wall;
	}
	return freeSurfaceFlag2D::empty;
}

void writeResults(FreeSurfaceFields2D<T, DESCRIPTOR> *fields,
				  MultiBlockLattice2D<T, ADESCRIPTOR> adLattice,
				  MultiScalarField2D<plint> *tagMatrix, MultiScalarField2D<double> *disjoiningPressureField, plint iT)
{
	std::vector<T> isoLevels;
	isoLevels.push_back(0.5);

	// T coef = 1.0 / 3.0;
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
	vtkOut.writeData<float>(*copyConvert<plint, T>(*tagMatrix), "bubbleTags",
							1.0);
	vtkOut.writeData<float>(*copyConvert<double, T>(*disjoiningPressureField), "disjoiningPressure",
							1.0);

	ImageWriter<T> imageWriter("leeloo");

	imageWriter.writeScaledGif(createFileName(outDir + "smoothVOF", iT, 6), *smoothVF);

	imageWriter.writeScaledGif(createFileName(outDir + "gas", iT, 6), *adrho);
}

plint numRememberedVolumes = 1;

int main(int argc, char **argv)
{

	//feenableexcept(FE_INVALID | FE_OVERFLOW);

	plbInit(&argc, &argv);

	std::cout.precision(10);
	std::scientific(std::cout);

	// Command-line arguments

	if (argc != 2)
	{
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
	//            new SmagorinskyBGKdynamics<T, DESCRIPTOR>(param.omega, param.cSmago);

	Dynamics<T, DESCRIPTOR> *dynamics =
		new BGKdynamics<T, DESCRIPTOR>(param.omega);

	Dynamics<T, ADESCRIPTOR> *adynamics =
		new ADYNAMICSWS<T, ADESCRIPTOR>(param.adOmega);
	Dynamics<T, ADESCRIPTOR> *emptyDynamics = new NoDynamics<T, ADESCRIPTOR>();

	FreeSurfaceFields2D<T, DESCRIPTOR> fields(blockStructure, dynamics->clone(),
											  param.rho_LB, param.sigma_LB,
											  (T)-1, param.gVector_LB);
	MultiBlockLattice2D<T, ADESCRIPTOR> adLattice(param.nx, param.ny,
												  adynamics->clone());
	Array<T, 2> u0((T)0, (T)0);

	// Initialization
	initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(), param.gasIni_LB,
							u0);
	adLattice.initialize();

	pcout << "Setting up initial condition." << std::endl;

	Box2D bottom(0, param.nx - 1, 0, 0);
	Box2D top(0, param.nx - 1, param.ny - 1, param.ny - 1);
	Box2D lateral1(0, 0, 0, param.ny - 1);
	Box2D lateral2(param.nx - 1, param.nx - 1, 0, param.ny - 1);

	setToConstant(fields.flag, bottom, (int)freeSurfaceFlag2D::wall);
	setToConstant(fields.flag, top, (int)freeSurfaceFlag2D::wall);
	setToConstant(fields.flag, lateral1, (int)freeSurfaceFlag2D::wall);
	setToConstant(fields.flag, lateral2, (int)freeSurfaceFlag2D::wall);

	setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1), initialFluidFlags);
	// analyticalIniVolumeFraction(fields.volumeFraction, fields.flag, insideFluid,
	//                             32);

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

	// Main iteration loop.
	for (plint iT = iniIter; iT < param.maxIter; iT++)
	{
		if (iT % param.statIter == 0 || iT == param.maxIter - 1)
		{
			pcout << "At iteration " << iT << std::endl;
			T avE = computeAverageEnergy(fields.lattice);
			pcout << "Average kinetic energy: " << avE << std::endl;
			plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
			pcout << "Number of interface cells: " << numIntCells << std::endl;
			if (iT != iniIter)
			{
				pcout << "Time spent for each iteration: "
					  << global::timer("iteration").getTime() / (T)param.statIter
					  << std::endl;
				global::timer("iteration").reset();
			}
			pcout << std::endl;
		}

		if (iT % param.outIter == 0 || iT == param.maxIter - 1)
		{
			pcout << "Writing results at iteration " << iT << std::endl;
			global::timer("images").restart();
			writeResults(&fields, adLattice, bubbleTracking.getTagMatrix(), bubbleTracking.getDisjoiningPressureField(), iT);
			global::timer("images").stop();
			pcout << "Time spent for writing results: "
				  << global::timer("images").getTime() << std::endl;
			pcout << std::endl;
		}

		global::timer("iteration").start();
		T bubbleVolumeRatio = iT == 0 ? 1.0 : param.bubbleVolumeRatio;
		bubbleTracking.execute<T, DESCRIPTOR, ADESCRIPTOR>(fields.volumeFraction, fields.flag, fields.normal, fields.rhoBar, fields.mass, fields.j, adLattice, oldvof, false, param.pi_LB);

		oldvof = fields.volumeFraction;

		bubbleGrowth.transition(
			bubbleTracking, iT, param.temperature, param.R_s, param.p_ini, 1.,
			param.rho_LB, bubbleVolumeRatio, param.entrapBubbles, numRememberedVolumes);
		bubbleGrowth.updateBubbleGrowth(fields.outsideDensity, param.rho_LB,
										 param.alpha, param.beta, (T)1.);

		// Free surface
		fields.lattice.executeInternalProcessors();
		fields.lattice.evaluateStatistics();
		fields.lattice.incrementTime();

		if (iT == 0 && param.freezeLargestBubble)
		{
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

		applyProcessingFunctional(
			new GrowthCoupling2D<T, ADESCRIPTOR, DESCRIPTOR>(
				adynamics->clone(), emptyDynamics->clone(), param.kh_LB, bubbleGrowth.getBubbles(), param.surfaceDiffusion, param.source_LB),
			adLattice.getBoundingBox(), couplingBlocks);

		adLattice.collideAndStream();

		global::timer("iteration").stop();

		if (iT % param.statIter == 0 || iT == param.maxIter - 1)
		{
			bubbleGrowth.timeHistoryLog(outDir + "bubbleTimeHistory.log");
			bubbleGrowth.fullBubbleLog(outDir + "FullBubbleRecord2D.log");

			// We do not log frozen bubbles.
			std::map<plint, BubbleInfo2D>::const_iterator it =
				bubbleGrowth.getBubbles().begin();
			T totalBubbleVolume = T();
			T currentDensity = T();
			T totalDisjoining = T();
			plint numBubbles = 0;
			for (; it != bubbleGrowth.getBubbles().end(); ++it)
			{
				if (it->second.isFrozen())
				{
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
