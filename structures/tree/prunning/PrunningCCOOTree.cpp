/*
 * PrunningCCOOTree.cpp
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#include "PrunningCCOOTree.h"

PrunningCCOOTree::PrunningCCOOTree(SingleVesselCCOOTree *tree){
	this->tree = tree;
}

PrunningCCOOTree::~PrunningCCOOTree(){
}

void PrunningCCOOTree::pruneTree(AbstractPrunningRule *rule, SingleVessel *root){

	if(root && rule->needsPrunning(root)){
		pruneBranch(root);
	}
	else{
		for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
			SingleVessel *childVessel = (SingleVessel *)(*it);
			pruneTree(rule, childVessel);
		}
	}
}

void PrunningCCOOTree::pruneBranch(SingleVessel *branch){

	this->tree->remove(branch);

	//	Update post-order nLevel, flux, initial resistances and intial betas.
	this->tree->updateTree((SingleVessel *) this->tree->root, this->tree);

	double maxVariation = INFINITY;
	while (maxVariation > this->tree->variationTolerance) {
		this->tree->updateTreeViscositiesBeta((SingleVessel *) this->tree->root, &maxVariation);
	}

}
