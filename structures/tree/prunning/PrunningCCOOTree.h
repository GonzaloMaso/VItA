/*
 * PrunningCCOOTree.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNNING_PRUNNINGCCOOTREE_H_
#define STRUCTURES_TREE_PRUNNING_PRUNNINGCCOOTREE_H_

#include "../SingleVesselCCOOTree.h"
#include "../vascularElements/AbstractVascularElement.h"
#include "AbstractPrunningRule.h"

/**
 * Prunes the branches of a tree that meet a certain prunning criterion described by AbstractPrunningRule objects.
 */
class PrunningCCOOTree {
	SingleVesselCCOOTree *tree;
public:
	/**
	 * Constructor.
	 * @param tree Tree to be prunned.
	 */
	PrunningCCOOTree(SingleVesselCCOOTree *tree);
	virtual ~PrunningCCOOTree();

	/**
	 * Prunes a tree by in order recursion.
	 * @param rule Rule used to determine if a branch must be prunned.
	 * @param root root of the branch to be analysed for prunning.
	 */
	void pruneTree(AbstractPrunningRule *rule, SingleVessel *root);

private:
	/**
	 * Prunes a given branch after being tag as "to be prunned" by pruneTree method.
	 * @param branch Branch to be prunned.
	 */
	void pruneBranch(SingleVessel *branch);
};

#endif /* STRUCTURES_TREE_PRUNNING_PRUNNINGCCOOTREE_H_ */
