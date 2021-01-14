/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GeneratorData.h
 *
 *  Created on: Mar 12, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef GENERATORDATA_H_
#define GENERATORDATA_H_

class AbstractCostEstimator;
#include <ostream>
//#include "../structures/tree/AbstractCostEstimator.h"

using namespace std;

/**
 * Class that acts as a wrapper of the parameters associated to a tree generation process.
 */
class GeneratorData {
private:
	bool didAllocateCostEstimator;
public:
	/**
	 * Constructor with default parameters.
	 */
	GeneratorData();
	/**
	 * Constructor with user parameters.
	 * @param nLevelTest Levels for tree scaling for each new segment test.
	 * @param nTerminalTrial Number of trials before diminish dlim.
	 * @param dLimReductionFactor Factor by which the Dlim constraint diminish after N failed trials.
	 * @param perfusionAreaFactor Factor that scales the perfusion area by which Dlim is computed.
	 * @param closeNeighborhoodFactor Factor that increase the neighborhood to search nearest neighbors.
	 * @param midPointDlimFactor Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
	 * @param nBifurcationTest Number of bifurcation sites tested in the optimization process.
	 */
	GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor, int nBifurcationTest);
	/**
	 * Constructor with user parameters.
	 * @param nLevelTest Levels for tree scaling for each new segment test.
	 * @param nTerminalTrial Number of trials before diminish dlim.
	 * @param dLimReductionFactor Factor by which the Dlim constraint diminish after N failed trials.
	 * @param perfusionAreaFactor Factor that scales the perfusion area by which Dlim is computed.
	 * @param closeNeighborhoodFactor Factor that increase the neighborhood to search nearest neighbors.
	 * @param midPointDlimFactor Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
	 * @param nBifurcationTest Number of bifurcation sites tested in the optimization process.
	 * @param vesselFunction Functionality of the generated vessels.
	 */
	GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor, int nBifurcationTest, int vesselFunction);
	/**
	 * Constructor with user parameters.
	 * @param nLevelTest Levels for tree scaling for each new segment test.
	 * @param nTerminalTrial Number of trials before diminish dlim.
	 * @param dLimReductionFactor Factor by which the Dlim constraint diminish after N failed trials.
	 * @param perfusionAreaFactor Factor that scales the perfusion area by which Dlim is computed.
	 * @param closeNeighborhoodFactor Factor that increase the neighborhood to search nearest neighbors.
	 * @param midPointDlimFactor Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
	 * @param nBifurcationTest Number of bifurcation sites tested in the optimization process.
	 * @param vesselFunction Functionality of the generated vessels.
	 * @param resetDLim Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
	 */
	GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor, int nBifurcationTest, int vesselFunction, bool resetDLim);
	/**
	 * Constructor with user parameters.
	 * @param nLevelTest Levels for tree scaling for each new segment test.
	 * @param nTerminalTrial Number of trials before diminish dlim.
	 * @param dLimReductionFactor Factor by which the Dlim constraint diminish after N failed trials.
	 * @param perfusionAreaFactor Factor that scales the perfusion area by which Dlim is computed.
	 * @param closeNeighborhoodFactor Factor that increase the neighborhood to search nearest neighbors.
	 * @param midPointDlimFactor Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
	 * @param nBifurcationTest Number of bifurcation sites tested in the optimization process.
	 * @param vesselFunction Functionality of the generated vessels.
	 * @param resetDLim Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
	 * @param costEstimator Cost estimator used to compute the functional at the current stage.
	 */
	GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor, int nBifurcationTest, int vesselFunction, bool resetDLim, AbstractCostEstimator *costEstimator);
	/**
	 * Destructor
	 */
	~GeneratorData();
	/**
	 * Levels for tree scaling for each new segment test.
	 */
	int nLevelTest;
	/**
	 * Number of trials before diminish dlim.
	 */
	int nTerminalTrial;
	/**
	 * Factor by which the Dlim constraint diminish after N failed trials.
	 */
	double dLimReductionFactor;
	/**
	 * Factor by which the Dlim constraint diminish at the first trial for each terminal estimation.
	 */
	double dLimCorrectionFactor;
	/**
	 * Factor that scales the perfusion area by which Dlim is computed.
	 */
	double perfusionAreaFactor;
	/**
	 * Factor that increase the neighborhood to search nearest neighbors.
	 */
	double closeNeighborhoodFactor;
	/**
	 * Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
	 */
	double midPointDlimFactor;
	/**
	 * Number of bifurcation sites tested in the optimization process is given by @p nBifurcationTest * ( @p nBifurcationTest - 1 ). (default 8)
	 */
	int nBifurcationTest;
	/**
	 * Functionality of the vessel generated, important for Object trees.
	 */
	int vesselFunction;
	/**
	 * Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
	 */
	bool resetsDLim;
	/**
	 * Cost estimator for the given stage.
	 */
	AbstractCostEstimator *costEstimator;

};

inline ostream& operator<<(ostream& os, GeneratorData *instanceData) {
	os << "nLevelTest : " << instanceData->nLevelTest << endl;
	os << "Bifurcation tries : " << instanceData->nBifurcationTest << endl;
	os << "Tries before dlim diminution : " << instanceData->nTerminalTrial << endl;
	os << "Dlim diminution factor : " << instanceData->dLimReductionFactor << endl;
	os << "Dlim correction factor : " << instanceData->dLimCorrectionFactor << endl;
	os << "Perfusion area factor : " << instanceData->perfusionAreaFactor << endl;
	os << "Close neighborhood factor : " << instanceData->closeNeighborhoodFactor << endl;
//	os << "Cost estimator : " << typeid(instanceData->costEstimator).name() << endl;
	return os;
}

#endif /* GENERATORDATA_H_ */
