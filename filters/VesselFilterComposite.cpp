/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilter.cpp
 *
 *  Created on: 22/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "VesselFilterComposite.h"

#include "../structures/vascularElements/SingleVessel.h"

VesselFilterComposite::VesselFilterComposite(vector<AbstractVesselFilter *> filters) : AbstractVesselFilter(){
	this->filters = filters;
}

VesselFilterComposite::~VesselFilterComposite(){
}

vector<SingleVessel *> VesselFilterComposite::apply(vector<SingleVessel *> vessels){
	vector<SingleVessel *> filteredVessels = vessels;

	for (std::vector<AbstractVesselFilter *>::iterator it = filters.begin(); it != filters.end(); ++it) {
		filteredVessels = (*it)->apply(filteredVessels);
	}

	return filteredVessels;
}
