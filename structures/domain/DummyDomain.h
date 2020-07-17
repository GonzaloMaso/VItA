/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * DummyDomain.h
 *
 *  Created on: Jun 8, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DUMMYDOMAIN_H_
#define DUMMYDOMAIN_H_

#include "../CCOCommonStructures.h"
#include <vector>
#include "SimpleDomain.h"

/**
 * Creates a Domain that does not generate random points, instead it loads them from the constructor call.
 */
class DummyDomain: public SimpleDomain {
public:
	/**
	 * Constructor that load the domain shape from a vtkPolydata representation within @p filename
	 * and a preloaded set of points contained in @p points.
	 * @param filename File of the domain representation.
	 * @param points Set of preloaded points in the domain.
	 */
	DummyDomain(string filename, vector<point> points, GeneratorData *instanceData);
	/**
	 * Return the last preloaded point.
	 * @return Last preloaded point.
	 */
	virtual point getRandomPoint();

	void logDomainFiles(FILE *fp);
};

#endif /* DUMMYDOMAIN_H_ */
