/*
 * OldLeavesPrunningRule.cpp
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#include "OldLeavesPruningRule.h"

OldLeavesPruningRule::OldLeavesPruningRule(int oldLeaveStage){
	this->oldLeaveStage = oldLeaveStage;
}

OldLeavesPruningRule::~OldLeavesPruningRule(){
}

bool OldLeavesPruningRule::needsPruning(SingleVessel *root){

	return root->getChildren().empty() && root->stage == this->oldLeaveStage;
}
