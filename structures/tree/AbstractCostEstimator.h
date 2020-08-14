/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractCostEstimator.h
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_ABSTRACTCOSTESTIMATOR_H_
#define TREE_ABSTRACTCOSTESTIMATOR_H_

#include "AbstractObjectCCOTree.h"
#include "../vascularElements/AbstractVascularElement.h"

/**
 * Estimates the functional cost of a tree.
 */
class AbstractCostEstimator {
public:
	/**
	 * Common constructor.
	 */
	AbstractCostEstimator();
	/**
	 * Common destructor.
	 */
	virtual ~AbstractCostEstimator();
	/**
	 * Clones the current estimator instance.
	 * @return Cloned instance.
	 */
	virtual AbstractCostEstimator *clone() = 0;

	/**
	 * Extracts information of the tree at the previous step.
	 * @param root Root of the tree at the previous step.
	 * @param parent	Vascular element where the new vessel will be connected.
	 * @param iNew	Distal position of the new vessel.
	 * @param iTest	Proximal position of the new vessel.
	 * @param dLim	Minimum radius distance from the new vessel to the tree.
	 */
	virtual void previousState(AbstractObjectCCOTree *tree, AbstractVascularElement *parent, point iNew, point iTest, double dLim) = 0;

	/**
	 * Computes the functional cost of the given tree.
	 * @param root	Root of the tree at the current step.
	 * @return Cost of the given tree.
	 */
	virtual double computeCost(AbstractObjectCCOTree *tree) = 0;
	/**
	 * Virtual method to be used by StagedFRROTreeGeneratorLogger
	 */	 
	virtual void logCostEstimator(FILE *fp) = 0;
};

#endif /* TREE_ABSTRACTCOSTESTIMATOR_H_ */
