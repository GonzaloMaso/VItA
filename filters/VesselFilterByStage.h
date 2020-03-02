/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilterByStage.h
 *
 *  Created on: 29/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef FILTERS_VESSELFILTERBYSTAGE_H_
#define FILTERS_VESSELFILTERBYSTAGE_H_

#include "AbstractVesselFilter.h"

class VesselFilterByStage: public AbstractVesselFilter {
	int stage;
public:
	VesselFilterByStage(int stage);
	virtual ~VesselFilterByStage();

	vector<SingleVessel *> apply(vector<SingleVessel *> vessels);
};

#endif /* FILTERS_VESSELFILTERBYSTAGE_H_ */
