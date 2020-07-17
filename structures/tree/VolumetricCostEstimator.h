/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VolumetricCostEstimation.h
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_VOLUMETRICCOSTESTIMATOR_H_
#define TREE_VOLUMETRICCOSTESTIMATOR_H_

#include "AbstractCostEstimator.h"

#include "../vascularElements/AbstractVascularElement.h"

/**
 * Cost estimator that computes the tree volume.
 */
class VolumetricCostEstimator: public AbstractCostEstimator {
	/**	Volume at the previous step. */
	double previousVolume;
public:
	/**
	 * Common constructor.
	 */
	VolumetricCostEstimator();
	/**
	 * Common destructor.
	 */
	virtual ~VolumetricCostEstimator();
	/**
	 * Clones the current estimator instance.
	 * @return Cloned instance.
	 */
	AbstractCostEstimator *clone();

	/**
	 * Extracts information of the tree at the previous step.
	 * @param root	Root of the tree at the previous step.
	 * @param parent	Vascular element where the new vessel will be connected.
	 * @param iNew	Distal position of the new vessel.
	 * @param iTest	Proximal position of the new vessel.
	 * @param dLim	Minimum radius distance from the new vessel to the tree.
	 */
	void previousState(AbstractObjectCCOTree *tree, AbstractVascularElement *parent, point iNew, point iTest, double dLim);

	/**
	 * Computes the functional cost of the given tree.
	 * @param tree	Tree at the current step.
	 * @return Cost of the given tree.
	 */
	double computeCost(AbstractObjectCCOTree *tree);

	void logCostEstimator(FILE *fp);

private:
	/**
	 * Computes the volume for the tree with root @p root.
	 * @param root	Root of the tree.
	 * @return	Volume of the tree.
	 */
	double computeTreeCost(AbstractVascularElement* root);
};

#endif /* TREE_VOLUMETRICCOSTESTIMATOR_H_ */
