/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilterByBranchingMode.cpp
 *
 *  Created on: 22/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "VesselFilterByBranchingMode.h"

VesselFilterByBranchingMode::VesselFilterByBranchingMode(AbstractVascularElement::BRANCHING_MODE mode) : AbstractVesselFilter() {
	this->mode = mode;
}

VesselFilterByBranchingMode::~VesselFilterByBranchingMode(){
	// TODO Auto-generated destructor stub
}

vector<SingleVessel*> VesselFilterByBranchingMode::apply(vector<SingleVessel*> vessels){
	vector<SingleVessel*> filteredVessels;
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if( (*it)->branchingMode == mode){
			filteredVessels.push_back(*it);
		}
	}
	return filteredVessels;
}
