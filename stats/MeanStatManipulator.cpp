/*
 * MeanStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "MeanStatManipulator.h"

MeanStatManipulator::MeanStatManipulator() : AbstractStatManipulator(){
}

double MeanStatManipulator::compute(vector<vessel*> vessels, VesselHandler::ATTRIBUTE att){
	double mean = 0.0;
	for (unsigned i = 0; i < vessels.size(); ++i) {
		mean += handler->getVesselAttribute(vessels[i],att);
	}
	mean /= vessels.size();
	return mean;
}
