/*
 * AbstractPrunningRule.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNNING_ABSTRACTPRUNNINGRULE_H_
#define STRUCTURES_TREE_PRUNNING_ABSTRACTPRUNNINGRULE_H_

#include "../../vascularElements/SingleVessel.h"

/**
 * Abstract class used in for the definition of rules to prune branches of a tree.
 */
class AbstractPrunningRule {
public:
	AbstractPrunningRule();
	virtual ~AbstractPrunningRule();

	/**
	 * If true the branch (root and all daughter vessels) are marked for punning.
	 * @param root Root vessel where the prunning is applied.
	 * @return If the branch has to be prunned.
	 */
	virtual bool needsPrunning(SingleVessel *root) = 0;
};

#endif /* STRUCTURES_TREE_PRUNNING_ABSTRACTPRUNNINGRULE_H_ */
