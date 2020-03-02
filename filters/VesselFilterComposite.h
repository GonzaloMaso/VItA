/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilter.h
 *
 *  Created on: 22/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef FILTERS_VESSELFILTERCOMPOSITE_H_
#define FILTERS_VESSELFILTERCOMPOSITE_H_

#include <vector>

#include "../structures/vascularElements/AbstractVascularElement.h"
#include "AbstractVesselFilter.h"

using namespace std;

/**
 * Filters a set of vessels by a set of filters defined in the constructor. Implemented with Composite pattern.
 */
class VesselFilterComposite : public AbstractVesselFilter{
	vector<AbstractVesselFilter *> filters;
public:
	VesselFilterComposite(vector<AbstractVesselFilter *> filters);
	virtual ~VesselFilterComposite();

	vector<SingleVessel *> apply(vector<SingleVessel *> vessels);
};

#endif /* FILTERS_VESSELFILTERCOMPOSITE_H_ */
