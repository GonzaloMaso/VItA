/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VolumetricCostEstimation.cpp
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "VolumetricCostEstimator.h"
#include "../vascularElements/SingleVessel.h"

VolumetricCostEstimator::VolumetricCostEstimator() : AbstractCostEstimator(){
	previousVolume = 0.0;
}

VolumetricCostEstimator::~VolumetricCostEstimator(){
}

double VolumetricCostEstimator::computeCost(AbstractObjectCCOTree *tree){
	return computeTreeCost(tree->getRoot()) - previousVolume;
}

void VolumetricCostEstimator::previousState(AbstractObjectCCOTree *tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
	previousVolume = ((SingleVessel *) tree->getRoot())->treeVolume;
}

double VolumetricCostEstimator::computeTreeCost(AbstractVascularElement* root) {

	double currentCost = ((SingleVessel *)root)->getVolume();
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* VolumetricCostEstimator::clone(){
	return (new VolumetricCostEstimator());
}

void VolumetricCostEstimator::logCostEstimator(FILE *fp) {
	fprintf(fp, "This domain uses VolumetricCostEstimator.\n");
}