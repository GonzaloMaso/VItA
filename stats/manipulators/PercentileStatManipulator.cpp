/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * PercentileDiameterStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "PercentileStatManipulator.h"

PercentileStatManipulator::PercentileStatManipulator(int numPercentile) : Abstract0DStatManipulator()
{
	this->numPercentile = numPercentile;
}

double PercentileStatManipulator::compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att){
	vector<double> values;
	for (unsigned long long i = 0; i < vessels.size(); ++i) {
		values[i] = handler->getVesselAttribute(vessels[i],att);
	}
	int position = vessels.size()/numPercentile;
	return values[position];
}
