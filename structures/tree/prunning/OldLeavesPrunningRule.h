/*
 * OldLeavesPrunningRule.h
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNNINGRULE_H_
#define STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNNINGRULE_H_

#include "AbstractPrunningRule.h"

/**
 * Rule that evaluates if the tree leaves belong to a certain old stage. The purpose of this rule is to remove needle-like
 * vessels from earlier stages that do not sprout in further stages, resulting in slender vessels.
 */
class OldLeavesPrunningRule: public AbstractPrunningRule {
	int oldLeaveStage;
public:
	OldLeavesPrunningRule(int oldLeaveStage);
	virtual ~OldLeavesPrunningRule();

	bool needsPrunning(SingleVessel *root);
};

#endif /* STRUCTURES_TREE_PRUNNING_OLDLEAVESPRUNNINGRULE_H_ */
