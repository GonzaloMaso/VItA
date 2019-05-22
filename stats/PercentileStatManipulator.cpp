/*
 * PercentileDiameterStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "PercentileStatManipulator.h"

PercentileStatManipulator::PercentileStatManipulator(int numPercentile)
{
	this->numPercentile = numPercentile;
}

double PercentileStatManipulator::compute(vector<vessel*> vessels, VesselHandler::ATTRIBUTE att){
	vector<double> values;
	for (unsigned long long i = 0; i < vessels.size(); ++i) {
		values[i] = handler->getVesselAttribute(vessels[i],att);
	}
	int position = vessels.size()/numPercentile;
	return values[position];
}
