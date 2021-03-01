/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * IntersectionVascularizedDomain.h
 *
 *  Created on: Mar 20, 2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_INTERSECTIONVASCULARIZEDDOMAIN_H_
#define DOMAIN_INTERSECTIONVASCULARIZEDDOMAIN_H_

#include <iostream>
#include <stdio.h>
#include <deque>
#include <random>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkOBBTree.h>

#include "AbstractDomain.h"
#include "../CCOCommonStructures.h"

using namespace std;

/**
 * This domain vascularises just the intersection between all vascularised regions.
 */
class IntersectionVascularizedDomain : public AbstractDomain {
	vector<string> filenameVR;
	/** Vascularized subdomains represented by vtkPolydata. */
	vector<vtkSmartPointer<vtkPolyData> > vtkVascularizedRegions;
	/** Cell locators responsible to determine if a segment is inside of each domain. */
	vector<vtkSmartPointer<vtkOBBTree> > sublocators;

	/**	Bounding box of the domain */
	double *boundingBox;
	/** Characteristic length for the domain. */
	double characteristicLength;
	/**	Amount of points randomly generated when no more points are available. */
	int nDraw;
	/** Random instance seed. */
	int seed = -1;
	/**	Random generator. */
	mt19937 generator;

public:
	/**
	 * Constructs a domain shaped by @p filenameHull (transport domain) with vascularized and non-vascularized regions contained in
	 * @p filenameVascularRegions and @p filenameNonVascularRegions respectively.
	 * @param filenameHull Domain description.
	 * @param filenameVascularRegions Vascularized regions descriptions.
	 * @param filenameNonVascularRegions Non-vascularized regions descriptions.
	 */
	IntersectionVascularizedDomain(vector<string> filenameVascularRegions, GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by @p filenameHull (transport domain) with vascularized and non-vascularized regions contained in
	 * @p filenameVascularRegions and @p filenameNonVascularRegions respectively. Also, itgenerates batches of
	 * @p nDraw random points each time that the domain is with no more points.
	 * @param filenameVascularRegions Vascularized regions descriptions.
	 * @param nDraw Random points generated at each batch.
	 */
	IntersectionVascularizedDomain(vector<string> filenameVascularRegions, int nDraw,GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by @p filenameHull (transport domain) with vascularized and non-vascularized regions contained in
	 * @p filenameVascularRegions and @p filenameNonVascularRegions respectively. Also, itgenerates batches of
	 * @p nDraw random points each time that the domain is with no more points.
	 * @param filenameVascularRegions Vascularized regions descriptions.
	 * @param nDraw Random points generated at each batch.
	 * @param seed Seed for random points generation.
	 */
	IntersectionVascularizedDomain(vector<string> filenameVascularRegions, int nDraw, int seed, GeneratorData *instanceData);
	/**
	 * Desctructor
	 */
	~IntersectionVascularizedDomain();
	/**
	 * Returns if the segment defined by the vertexes @p xs and @p xf is inside the current transport domain without intersecting non-vascularized regions.
	 * @param xs Start point of the segment.
	 * @param xf End point of the segment.
	 * @return 1 if the segment defined by the vertexes @p xs and @p xf is inside the current domain otherwise 0.
	 */
	int isSegmentInside(point xs, point xf);
	/**
	 * Estimates a characteristic length for the current domain. This length is useful to estimate the perfusion volume of the domain.
	 * @return Chracteristic length.
	 */
	double getCharacteristicLength();
	/**
	 * Minimum distance between a new vessel terminal and the tree based on the current terminal perfusion. This is lower bound criteria for the vessel length and also aids
	 * toward a spatially homogeneous terminal distribution.
	 * @param nVessels Quantity of terminals in the current tree.
	 * @param factor Factor used to scale a perfusion measure in order to estimate the DLim value.
	 * @return DLim value.
	 */
	double getDLim(long long int nVessels, double factor);
	/**
	 * Returns the local neighbors to the @p p point locus. The amount of neighbors is estimated based on the current terminal perfusion.
	 * @param p Central point of the neighborhood.
	 * @param nVessels Amount of terminals in the tree.
	 * @return Array of neighbor vessels.
	 */
	double *getLocalNeighborhood(point p, long long int nVessels);
	/**
	 * Computes the size of the domain.
	 * @return Size of the domain.
	 */
	double getSize();
	/**
	 * Return a random point inside the domain.
	 * @return Point inside the domain.
	 */
	virtual point getRandomPoint();
	/**
	 * Returns the locator used to test inside segments.
	 * @return Locator used to test inside segments.
	 */
	vtkSmartPointer<vtkOBBTree>& getLocator();
	/**
	 * Returns the number of random points created in each batch of random generations.
	 * @return Number of random points created in each batch of random generations.
	 */
	int getDraw();
	/**
	 * Return a set of inner domain points.
	 * @return Set of inner domain points.
	 */
	deque<point>& getRandomInnerPoints();
	/**
	 * Returns the vtkPolydata with the domain representation.
	 * @return vtkPolydata with the domain representation.
	 */
	vtkSmartPointer<vtkPolyData>& getVtkGeometry();
	/**
	 * Returns random generator seed. 
	 * @return @p seed
	 */
	int getSeed();

	vector<string> getFilenameVR();

	void logDomainFiles(FILE *fp);

	/**
	 * Returns vtkSelectEnclosedPoints for the first geometry passed to the domain.
	 */
	vtkSmartPointer<vtkSelectEnclosedPoints> getEnclosedPoints() override;

protected:
	deque<point> randomInnerPoints;
	void generateRandomPoints();
	void removeRandomOuterPoints();

};

#endif /* DOMAIN_INTERSECTIONVASCULARIZEDDOMAIN_H_ */
