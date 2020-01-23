/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "Abstract1DStatManipulator.h"

Abstract1DStatManipulator::Abstract1DStatManipulator()
{
	this->handler = new VesselObjectHandler();
}

Abstract1DStatManipulator::~Abstract1DStatManipulator(){
	delete handler;
}
