/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeStatsManager.h
 *
 *  Created on: Feb 9, 2018
 *      Author: gonzalo
 */

#ifndef STATS_OBJECTTREESTATSMANAGER_H_
#define STATS_OBJECTTREESTATSMANAGER_H_

#include "../structures/tree/AbstractObjectCCOTree.h"
#include "VesselObjectHandler.h"

/**
 * Facade class used to access to all statistical analysis services implemented.
 * The class aims to return vectors with the statistical results and no plots or
 * any other outputs. Outputs are expected to be generated in the CCORendering library
 * to decouple GUI from data and model components. The last is key to support the
 * CCO library in clusters or environments with no GUI capabilities.
 */
class ObjectTreeStatsManager {
	/** Tree that undergo through statistical analysis. */
	AbstractObjectCCOTree * tree;
public:
	/**
	 * Creates an instance of the @class TreeStatsManager for the given @p tree.
	 * @param tree Tree for statistical analysis.
	 */
	ObjectTreeStatsManager(AbstractObjectCCOTree *tree);

	/**
	 * Load the vector @p levels, @p means and @p stds with the number of bifurcation level
	 * and mean and standard deviation of a vessel field across that levels.
	 * @param levels Levels of bifurcation (0,..., N-1).
	 * @param means Mean for each level.
	 * @param stds Standard deviation for each level.
	 * @param att Field over which the mean and std is computed.
	 */
	void getMeanPerLevel(vector<double> *levels, vector<double> *means, vector<double> *stds, VesselObjectHandler::ATTRIBUTE att);

	/**
	 * Load @p values with all attributes indicated in @p att.
	 * @param levels Levels of bifurcation (0,..., N-1).
	 * @param values Violin data for each level.
	 * @param att Field for sampling.
	 */
	void getAttributesPerLevel(vector<double> *levels, vector<vector<double>> *values, vector<VesselObjectHandler::ATTRIBUTE> att);

	/**
	 * Load @p values with all attributes indicated in @p att in the subtree rooted at @p root.
	 * @param root Subtree to be scanned.
	 * @param levels Levels of bifurcation (0,..., N-1).
	 * @param values Violin data for each level.
	 * @param att Field for sampling.
	 */
	void getAttributesPerLevel(SingleVessel *root, vector<double> *levels, vector<vector<double>> *values, vector<VesselObjectHandler::ATTRIBUTE> att);

	/**
	 * Load load @p values with the attributes @p att. Only segments from stage @p stage and its children are processed.
	 * @param branches Vector with the identifier of each sample associated to its branch.
	 * @param levels Levels of bifurcation (0,..., N-1).
	 * @param values Violin data for each level.
	 * @param att Field for sampling.
	 */
	void getBranchesAttributesPerLevel(vector<double> *banches, vector<double> *levels, vector<vector<double>> *values, vector<VesselObjectHandler::ATTRIBUTE> att, int stage);

};

#endif /* STATS_OBJECTTREESTATSMANAGER_H_ */
