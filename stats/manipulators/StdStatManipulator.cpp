/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StdDiameterStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "StdStatManipulator.h"

StdStatManipulator::StdStatManipulator() : Abstract0DStatManipulator(){
	this->meanManipulator = new MeanStatManipulator();
}

double StdStatManipulator::compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att)	{

	double mean = meanManipulator->compute(vessels,att);
	double std = 0.0;

	for (unsigned long long i = 0; i < vessels.size(); ++i) {
		double centroidDistance = handler->getVesselAttribute(vessels[i],att) - mean;
		std += centroidDistance*centroidDistance;
	}
	std = sqrt( std / (double) vessels.size() );
	return std;
}
