/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * UniformDistributionGenerator.cpp
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "UniformDistributionGenerator.h"

UniformDistributionGenerator::UniformDistributionGenerator() : DistributionGenerator(){
}

UniformDistributionGenerator::~UniformDistributionGenerator(){
}

void UniformDistributionGenerator::initialize(int seed, double *boundingBox) {
	DistributionGenerator::initialize(seed, boundingBox);
	distX = uniform_real_distribution<double>(boundingBox[0], boundingBox[1]);
	distY = uniform_real_distribution<double>(boundingBox[2], boundingBox[3]);
	distZ = uniform_real_distribution<double>(boundingBox[4], boundingBox[5]);
}

vector<point> UniformDistributionGenerator::getNPoints(int n){
	vector<point> randomInnerPoints;
	for (int i = 0; i < n; ++i) {
		point p = { distX(generator), distY(generator), distZ(generator) };
		randomInnerPoints.push_back(p);
	}

	return randomInnerPoints;
}
