/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "AbstractStructuredStatManipulator.h"

AbstractStructuredStatManipulator::AbstractStructuredStatManipulator()
{
	this->handler = new VesselStructHandler();
}

AbstractStructuredStatManipulator::~AbstractStructuredStatManipulator(){
	delete handler;
}
