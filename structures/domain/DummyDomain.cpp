/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * DummyDomain.cpp
 *
 *  Created on: Jun 8, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "DummyDomain.h"

/**
 * Creates a domain with the geometry of the file @param filename with the fixed
 * random points @param points.
 * @param filename VTP file containing the geometry for the domain.
 * @param points Points used for the tree generation.
 */
DummyDomain::DummyDomain(string filename, vector<point> points,GeneratorData *instanceData) : SimpleDomain(filename,instanceData)
{
	for(vector<point>::iterator it = points.begin(); it != points.end(); ++it)
		randomInnerPoints.push_back(*it);
}

/**
 * Returns the first random point and removes it from the domain.
 * @return
 */
point DummyDomain::getRandomPoint(){
	++pointCounter;
	if(randomInnerPoints.size()==0)
		return {0,0,0};
	else{
		point p = randomInnerPoints.front();
		randomInnerPoints.pop_front();
		return p;
	}
}

void DummyDomain::logDomainFiles(FILE *fp) {
	fprintf(fp, "DummyDomain\n");
    fprintf(fp, "filename = %s\n", this->getFilename().c_str());
}