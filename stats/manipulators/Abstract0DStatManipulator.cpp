/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractStatManipulator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "Abstract0DStatManipulator.h"

Abstract0DStatManipulator::Abstract0DStatManipulator()
{
	this->handler = new VesselObjectHandler();
}

Abstract0DStatManipulator::~Abstract0DStatManipulator(){
	delete handler;
}
