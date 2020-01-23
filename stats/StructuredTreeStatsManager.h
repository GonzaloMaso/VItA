/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeStatsManager.h
 *
 *  Created on: Feb 9, 2018
 *      Author: gonzalo
 */

#ifndef STATS_STRUCTUREDTREESTATSMANAGER_H_
#define STATS_STRUCTUREDTREESTATSMANAGER_H_

#include "../structures/tree/AbstractStructuredCCOTree.h"
#include "VesselStructHandler.h"

/**
 * Facade class used to access to all statistical analysis services implemented.
 * The class aims to return vectors with the statistical results and no plots or
 * any other outputs. Outputs are expected to be generated in the CCORendering library
 * to decouple GUI from data and model components. The last is key to support the
 * CCO library in clusters or environments with no GUI capabilities.
 */
class StructuredTreeStatsManager {
	/** Tree that undergo through statistical analysis. */
	AbstractStructuredCCOTree * tree;
public:
	/**
	 * Creates an instance of the @class TreeStatsManager for the given @p tree.
	 * @param tree Tree for statistical analysis.
	 */
	StructuredTreeStatsManager(AbstractStructuredCCOTree *tree);

	/**
	 * Load the vector @p levels, @p means and @p stds with the number of bifurcation level
	 * and mean and standard deviation of a vessel field across that levels.
	 * @param levels Levels of bifurcation (0,..., N-1).
	 * @param means Mean for each level.
	 * @param stds Standard deviation for each level.
	 * @param att Field over which the mean and std is computed.
	 */
	void getMeanPerLevel(vector<double> *levels, vector<double> *means, vector<double> *stds, VesselStructHandler::ATTRIBUTE att);

};

#endif /* STATS_STRUCTUREDTREESTATSMANAGER_H_ */
