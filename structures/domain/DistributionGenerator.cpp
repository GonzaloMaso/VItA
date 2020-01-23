/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * DistributionGenerator.cpp
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "DistributionGenerator.h"

DistributionGenerator::DistributionGenerator(){
}

DistributionGenerator::~DistributionGenerator(){
}

void DistributionGenerator::initialize(int seed, double* boundingBox){
	this->boundingBox = boundingBox;
	this->generator = mt19937(seed);
}
