/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilterByBranchingMode.h
 *
 *  Created on: 22/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef FILTERS_VESSELFILTERBYBRANCHINGMODE_H_
#define FILTERS_VESSELFILTERBYBRANCHINGMODE_H_

#include "AbstractVesselFilter.h"
#include "../structures/vascularElements/AbstractVascularElement.h"

class VesselFilterByBranchingMode : public AbstractVesselFilter{
	AbstractVascularElement::BRANCHING_MODE mode;
public:
	VesselFilterByBranchingMode(AbstractVascularElement::BRANCHING_MODE mode);
	virtual ~VesselFilterByBranchingMode();

	vector<SingleVessel *> apply(vector<SingleVessel *> vessels);
};

#endif /* FILTERS_VESSELFILTERBYBRANCHINGMODE_H_ */
