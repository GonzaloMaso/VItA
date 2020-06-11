/*
 * OldLeavesPrunningRule.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNINGRULE_H_
#define STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNINGRULE_H_

#include "AbstractPruningRule.h"

/**
 * Rule that evaluates if the tree leaves belong to a certain old stage. The purpose of this rule is to remove needle-like
 * vessels from earlier stages that do not sprout in further stages, resulting in slender vessels.
 */
class OldLeavesPruningRule: public AbstractPruningRule {
	int oldLeaveStage;
public:
	OldLeavesPruningRule(int oldLeaveStage);
	virtual ~OldLeavesPruningRule();

	bool needsPruning(SingleVessel *root);
};

#endif /* STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNINGRULE_H_ */
