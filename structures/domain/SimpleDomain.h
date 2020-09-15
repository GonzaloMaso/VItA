/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * Mesh.h
 *
 *  Created on: May 24, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef SIMPLEDOMAIN_H_
#define SIMPLEDOMAIN_H_

#include <iostream>
#include <stdio.h>
#include <deque>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkOBBTree.h>

#include "AbstractDomain.h"
#include "../CCOCommonStructures.h"
#include "DistributionGenerator.h"

using namespace std;

/**
 * Class that model 3D domains with homogeneous perfusion characteristics. It deals with
 * geometrical issues such as identify inner/outer points and segments, estimation of area/volume computation,
 * generation of inner random points and computation of perfusion volume.
 */
class SimpleDomain: public AbstractDomain {
	/** Name of the file from which domain is constructed */
	string filename;
	/** vtkPolydata description of the domain. */
	vtkSmartPointer<vtkPolyData> vtkGeometry;
	/** Cell locator responsible to determine if a segment is inside the domain. */
	vtkSmartPointer<vtkOBBTree> locator;
	/** Characteristic length for the domain. */
	double characteristicLength;
	/**	Amount of points randomly generated when no more points are available. */
	int nDraw;
	/** Random generator seed.*/
	int seed;
	/** Point generator */
	DistributionGenerator *distribution;
	bool didAllocateDistribution;
public:
	/**
	 * Constructs a domain shaped by the vtkPolydata within @p filename.
	 * @param filename File of the domain representation.
	 */
	SimpleDomain(string filename, GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by the vtkPolydata within @p filename which generates
	 * batches of @p nDraw random points each time that the domain is with no more points.
	 * @param filename File of the domain representation.
	 * @param nDraw Random points generated at each batch.
	 */
	SimpleDomain(string filename, int nDraw, GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by the vtkPolydata within @p filename which generates
	 * batches of @p nDraw random points each time that the domain is with no more points.
	 * @param filename File of the domain representation.
	 * @param nDraw Random points generated at each batch.
	 * @param seed Seed that determine a specific random instance for the generation of points.
	 */
	SimpleDomain(string filename, int nDraw, int seed, GeneratorData *instanceData);
	/**
	 * Constructs a domain shaped by the vtkPolydata within @p filename which generates
	 * batches of @p nDraw random points each time that the domain is with no more points.
	 * @param filename File of the domain representation.
	 * @param nDraw Random points generated at each batch.
	 * @param seed Seed that determine a specific random instance for the generation of points.
	 * @param instanceData Data of the current generation process.
	 * @param distribution Generator of terminal points following a specific distribution.
	 */
	SimpleDomain(string filename, int nDraw, int seed, GeneratorData *instanceData, DistributionGenerator *distribution);
	/**
	 * Returns if the segment defined by the vertexes @p xs and @p xf is inside the current domain.
	 * @param xs Start point of the segment.
	 * @param xf End point of the segment.
	 * @return 1 if the segment defined by the vertexes @p xs and @p xf is inside the current domain otherwise 0.
	 */
	/**
	 * Destructor
	 */
	~SimpleDomain();
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
	 * Returns random seed generator.
	 * @return @p seed.
	 */
	int getSeed();

	string getFilename();

	void logDomainFiles(FILE *fp);

	vtkSmartPointer<vtkSelectEnclosedPoints> getEnclosedPoints() override;

protected:
	deque<point> randomInnerPoints;
	void generateRandomPoints();
	void removeRandomOuterPoints();
};

#endif /* SIMPLEDOMAIN_H_ */
