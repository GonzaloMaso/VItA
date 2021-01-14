/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CompositeDistributionGenerator.cpp
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "CompositeDistributionGenerator.h"
#include <algorithm>

CompositeDistributionGenerator::CompositeDistributionGenerator(vector<DistributionGenerator *> distributions) : DistributionGenerator(){
	this->distributions = distributions;
}

CompositeDistributionGenerator::~CompositeDistributionGenerator(){
	this->distributions.clear();
}

void CompositeDistributionGenerator::initialize(int seed, double *boundingBox) {
	DistributionGenerator::initialize(seed, boundingBox);
	for (std::vector<DistributionGenerator *>::iterator it = distributions.begin(); it != distributions.end(); ++it) {
		(*it)->initialize(seed, boundingBox);
	}
}

vector<point> CompositeDistributionGenerator::getNPoints(int n){
	vector<point> randomInnerPoints;
	int nGenerators = distributions.size();
	int nPerGenerator = n/nGenerators;

	for (std::vector<DistributionGenerator *>::iterator it = distributions.begin(); it != distributions.end(); ++it) {
		vector<point> generatorPoints = (*it)->getNPoints(nPerGenerator);
		randomInnerPoints.insert(randomInnerPoints.end(), generatorPoints.begin(), generatorPoints.end());
	}

	shuffle(randomInnerPoints.begin(), randomInnerPoints.end(), generator);

	return randomInnerPoints;
}
