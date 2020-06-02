/*
 * OldLeavesPrunningRule.cpp
 *
 *  Created on: 2/06/2020
 *      Author: Gonzalo D. Maso Talou
 */

#include "OldLeavesPrunningRule.h"

OldLeavesPrunningRule::OldLeavesPrunningRule(int oldLeaveStage){
	this->oldLeaveStage = oldLeaveStage
}

OldLeavesPrunningRule::~OldLeavesPrunningRule(){
}

bool OldLeavesPrunningRule::needsPrunning(SingleVessel *root){

	return root->getChildren().empty() && root->stage == this->oldLeaveStage;
}
