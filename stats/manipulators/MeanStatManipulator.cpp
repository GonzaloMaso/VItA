/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "MeanStatManipulator.h"

MeanStatManipulator::MeanStatManipulator() : Abstract0DStatManipulator(){
}

double MeanStatManipulator::compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att){
	double mean = 0.0;
	for (unsigned i = 0; i < vessels.size(); ++i) {
		mean += handler->getVesselAttribute(vessels[i],att);
	}
	mean /= vessels.size();
	return mean;
}
