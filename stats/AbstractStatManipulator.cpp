/*
 * AbstractStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "AbstractStatManipulator.h"

AbstractStatManipulator::AbstractStatManipulator()
{
	handler = new VesselHandler();
}

AbstractStatManipulator::~AbstractStatManipulator(){
	delete handler;
}
