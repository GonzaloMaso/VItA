/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractVascularElement.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "AbstractVascularElement.h"

AbstractVascularElement::AbstractVascularElement(){
	stage = 0;
	qReservedFraction = 0.0;
	parent = NULL;
	vesselFunction = VESSEL_FUNCTION::DISTRIBUTION;
	branchingMode = BRANCHING_MODE::DEFORMABLE_PARENT;
	terminalType = TERMINAL_TYPE::COMMON;
	children.clear();
}

AbstractVascularElement::~AbstractVascularElement(){
	this->children.clear();
	this->vessels.clear();
}

vector<SingleVessel*>& AbstractVascularElement::getVessels(){
	return vessels;
}

void AbstractVascularElement::addChild(AbstractVascularElement *newChild) {
	children.push_back(newChild);
}

void AbstractVascularElement::removeChildren() {
	children.clear();
}

void AbstractVascularElement::setBranchingMode(BRANCHING_MODE mode) {
	this->branchingMode = mode;
}
