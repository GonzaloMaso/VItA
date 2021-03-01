/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * PartiallyVascularizedDomain.h
 *
 *  Created on: Mar 21, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_PARTIALLYVASCULARIZEDDOMAIN_H_
#define DOMAIN_PARTIALLYVASCULARIZEDDOMAIN_H_

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
 * This domain features three different type of volumes: (i) transport volume;
 * (ii) vascularized volumes; (iii) non-vascularized volumes. A transport volume
 * allows that existent vessels within it, create terminals to vascularized volumes.
 * Thus, its space may prone anastomose from existent arteries to vascularized regions,
 * but no terminals can be created in it. A vascularized volume allowed the anastomose
 * in its existent vessels and the creation of terminals within it. Lastly,
 * non-vascularized volumes do not allowed vessels to enter exit nor generate terminals.
 * It is mandatory that vascularized regions be within the transport volume.
 */
class PartiallyVascularizedDomain : public AbstractDomain {
	string filenameHull;
	vector<string> filenameVR;
	vector<string> filenameNVR;
	/** vtkPolydata description of the transport domain. */
	vtkSmartPointer<vtkPolyData> vtkTransportRegion;
	/** Vascularized subdomains represented by vtkPolydata. */
	vector<vtkSmartPointer<vtkPolyData> > vtkVascularizedRegions;
	/** Non-vascularized subdomains represented by vtkPolydata. */
	vector<vtkSmartPointer<vtkPolyData> > vtkHollowRegions;
	/** Cell locator responsible to determine if a segment is inside the domain. */
	vtkSmartPointer<vtkOBBTree> locator;
	/** Cell locator responsible to determine if a segment is inside a non-vascularized region. */
	vector<vtkSmartPointer<vtkOBBTree> > hollowLocators;

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
	PartiallyVascularizedDomain(string filenameHull, vector<string> filenameVascularRegions, vector<string> filenameNonVascularRegions, GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by @p filenameHull (transport domain) with vascularized and non-vascularized regions contained in
	 * @p filenameVascularRegions and @p filenameNonVascularRegions respectively. Also, itgenerates batches of
	 * @p nDraw random points each time that the domain is with no more points.
	 * @param filenameHull Domain description.
	 * @param filenameVascularRegions Vascularized regions descriptions.
	 * @param filenameNonVascularRegions Non-vascularized regions descriptions.
	 * @param nDraw Random points generated at each batch.
	 */
	PartiallyVascularizedDomain(string filenameHull, vector<string> filenameVascularRegions, vector<string> filenameNonVascularRegions,
			int nDraw,GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by @p filenameHull (transport domain) with vascularized and non-vascularized regions contained in
	 * @p filenameVascularRegions and @p filenameNonVascularRegions respectively. Also, itgenerates batches of
	 * @p nDraw random points each time that the domain is with no more points.
	 * @param filenameHull Domain description.
	 * @param filenameVascularRegions Vascularized regions descriptions.
	 * @param filenameNonVascularRegions Non-vascularized regions descriptions.
	 * @param nDraw Random points generated at each batch.
	 * @param seed Seed for random points generation.
	 */
	PartiallyVascularizedDomain(string filenameHull, vector<string> filenameVascularRegions, vector<string> filenameNonVascularRegions,
			int nDraw, int seed, GeneratorData *instanceData);
	/**
	 * Returns if the segment defined by the vertexes @p xs and @p xf is inside the current transport domain without intersecting non-vascularized regions.
	 * @param xs Start point of the segment.
	 * @param xf End point of the segment.
	 * @return 1 if the segment defined by the vertexes @p xs and @p xf is inside the current domain otherwise 0.
	 */
	/**
	 * Destructor
	 */
	~PartiallyVascularizedDomain();
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
	 * Saves the random generated points that have not been used yet.
	 * @param filename Output file writen in VTK format.
	 */
	void savePoints(string filename);
	/**
	 * Returns random generator seed. 
	 * @return @p seed
	 */
	int getSeed();

	string getFilenameHull();

	vector<string> getFilenameVR();

	vector<string> getFilenameNVR();

	void logDomainFiles(FILE *fp);

	vtkSmartPointer<vtkSelectEnclosedPoints> getEnclosedPoints() override;

protected:
	//	TODO This should be generalised in the AbstractDomain class
	deque<point> randomInnerPoints;
	void generateRandomPoints();
	void removeRandomOuterPoints();

};

#endif /* DOMAIN_PARTIALLYVASCULARIZEDDOMAIN_H_ */
