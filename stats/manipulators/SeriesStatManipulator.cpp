/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "SeriesStatManipulator.h"

SeriesStatManipulator::SeriesStatManipulator() : Abstract1DStatManipulator(){
}

vector<double> SeriesStatManipulator::compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att){

	vector<double> series;
	for (unsigned i = 0; i < vessels.size(); ++i) {
		series.push_back(handler->getVesselAttribute(vessels[i],att));
	}
	return series;
}
