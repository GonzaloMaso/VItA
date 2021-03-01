/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GeneratorData.cpp
 *
 *  Created on: Mar 12, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "GeneratorData.h"
#include "../structures/tree/VolumetricCostEstimator.h"

GeneratorData::GeneratorData()
{
	this->nLevelTest = 128;
	this->nTerminalTrial = 500;
	this->dLimReductionFactor = 0.9;
	this->perfusionAreaFactor = 0.5;
	this->closeNeighborhoodFactor = 4.0;
	this->midPointDlimFactor = 0.25;
	this->nBifurcationTest = 7;
	this->dLimCorrectionFactor = 1.0;
	this->vesselFunction = 0;
	this->resetsDLim = false;
	this->didAllocateCostEstimator = false;
}

GeneratorData::GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor,
		int nBifurcationTest){
	this->nLevelTest = nLevelTest;
	this->nTerminalTrial = nTerminalTrial;
	this->dLimReductionFactor = dLimReductionFactor;
	this->perfusionAreaFactor = perfusionAreaFactor;
	this->closeNeighborhoodFactor = closeNeighborhoodFactor;
	this->midPointDlimFactor = midPointDlimFactor;
	this->nBifurcationTest = nBifurcationTest;

	this->dLimCorrectionFactor = 1.0;
	this->vesselFunction = 0;
	this->resetsDLim = false;
	this->costEstimator = new VolumetricCostEstimator();
	this->didAllocateCostEstimator = true;
}

GeneratorData::GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor,
		int nBifurcationTest, int vesselFunction){
	this->nLevelTest = nLevelTest;
	this->nTerminalTrial = nTerminalTrial;
	this->dLimReductionFactor = dLimReductionFactor;
	this->perfusionAreaFactor = perfusionAreaFactor;
	this->closeNeighborhoodFactor = closeNeighborhoodFactor;
	this->midPointDlimFactor = midPointDlimFactor;
	this->nBifurcationTest = nBifurcationTest;

	this->dLimCorrectionFactor = 1.0;
	this->vesselFunction = vesselFunction;
	this->resetsDLim = false;
	this->costEstimator = new VolumetricCostEstimator();
	this->didAllocateCostEstimator = true;
}

GeneratorData::GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor,
		int nBifurcationTest, int vesselFunction, bool resetDLim){
	this->nLevelTest = nLevelTest;
	this->nTerminalTrial = nTerminalTrial;
	this->dLimReductionFactor = dLimReductionFactor;
	this->perfusionAreaFactor = perfusionAreaFactor;
	this->closeNeighborhoodFactor = closeNeighborhoodFactor;
	this->midPointDlimFactor = midPointDlimFactor;
	this->nBifurcationTest = nBifurcationTest;

	this->dLimCorrectionFactor = 1.0;
	this->vesselFunction = vesselFunction;
	this->resetsDLim = resetDLim;
	this->costEstimator = new VolumetricCostEstimator();
	this->didAllocateCostEstimator = true;
}

GeneratorData::GeneratorData(int nLevelTest, int nTerminalTrial, double dLimReductionFactor, double perfusionAreaFactor, double closeNeighborhoodFactor, double midPointDlimFactor,
		int nBifurcationTest, int vesselFunction, bool resetDLim, AbstractCostEstimator *costEstimator){
	this->nLevelTest = nLevelTest;
	this->nTerminalTrial = nTerminalTrial;
	this->dLimReductionFactor = dLimReductionFactor;
	this->perfusionAreaFactor = perfusionAreaFactor;
	this->closeNeighborhoodFactor = closeNeighborhoodFactor;
	this->midPointDlimFactor = midPointDlimFactor;
	this->nBifurcationTest = nBifurcationTest;

	this->dLimCorrectionFactor = 1.0;
	this->vesselFunction = vesselFunction;
	this->resetsDLim = resetDLim;
	this->costEstimator = costEstimator;
	this->didAllocateCostEstimator = false;
}

GeneratorData::~GeneratorData() {
	if(didAllocateCostEstimator) {
		delete this->costEstimator;
	}
}