/*
 * AbstractPrunningRule.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNING_ABSTRACTPRUNINGRULE_H_
#define STRUCTURES_TREE_PRUNING_ABSTRACTPRUNINGRULE_H_

#include "../../vascularElements/SingleVessel.h"

/**
 * Abstract class used in for the definition of rules to prune branches of a tree.
 */
class AbstractPruningRule {
public:
	AbstractPruningRule();
	virtual ~AbstractPruningRule();

	/**
	 * If true the branch (root and all daughter vessels) are marked for puning.
	 * @param root Root vessel where the pruning is applied.
	 * @return If the branch has to be pruned.
	 */
	virtual bool needsPruning(SingleVessel *root) = 0;
};

#endif /* STRUCTURES_TREE_PRUNING_ABSTRACTPRUNNINGRULE_H_ */
