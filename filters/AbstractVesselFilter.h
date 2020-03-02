/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractVesselFilter.h
 *
 *  Created on: 22/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef FILTERS_ABSTRACTVESSELFILTER_H_
#define FILTERS_ABSTRACTVESSELFILTER_H_

#include <vector>
#include "../structures/vascularElements/SingleVessel.h"

using namespace std;


class AbstractVesselFilter {
public:
	AbstractVesselFilter();
	virtual ~AbstractVesselFilter();

	virtual vector<SingleVessel *> apply(vector<SingleVessel *> vessels) = 0;
};

#endif /* FILTERS_ABSTRACTVESSELFILTER_H_ */
