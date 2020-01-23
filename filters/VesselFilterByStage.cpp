/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilterByStage.cpp
 *
 *  Created on: 29/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "VesselFilterByStage.h"

VesselFilterByStage::VesselFilterByStage(int stage) : AbstractVesselFilter(){
	this->stage = stage;
}

VesselFilterByStage::~VesselFilterByStage(){
}

vector<SingleVessel*> VesselFilterByStage::apply(vector<SingleVessel*> vessels){
	vector<SingleVessel*> filteredVessels;
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if( (*it)->stage == stage){
			filteredVessels.push_back(*it);
		}
	}
	return filteredVessels;
}
