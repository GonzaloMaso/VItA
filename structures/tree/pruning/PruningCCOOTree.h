/*
 * PrunningCCOOTree.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNNING_PRUNINGCCOOTREE_H_
#define STRUCTURES_TREE_PRUNNING_PRUNINGCCOOTREE_H_

#include "../SingleVesselCCOOTree.h"
#include "../../vascularElements/AbstractVascularElement.h"
#include "AbstractPruningRule.h"

/**
 * Prunes the branches of a tree that meet a certain prunning criterion described by AbstractPrunningRule objects.
 */
class PruningCCOOTree {
	SingleVesselCCOOTree *tree;
public:
	/**
	 * Constructor.
	 * @param tree Tree to be prunned.
	 */
	PruningCCOOTree(SingleVesselCCOOTree *tree);
	virtual ~PruningCCOOTree();

	/**
	 * Prunes a tree by in order recursion.
	 * @param rule Rule used to determine if a branch must be pruned.
	 * @param root root of the branch to be analysed for pruning.
	 */
	void pruneTree(AbstractPruningRule *rule, SingleVessel *root);

private:
	/**
	 * Prunes a given branch after being tag as "to be pruned" by pruneTree method.
	 * @param branch Branch to be pruned.
	 */
	void pruneBranch(SingleVessel *branch);
};

#endif /* STRUCTURES_TREE_PRUNNING_PRUNINGCCOOTREE_H_ */
