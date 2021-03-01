/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractDomain.h
 *
 *  Created on: 3 de dez de 2017
 *      Author: gonzalo
 */

#ifndef DOMAIN_ABSTRACTDOMAIN_H_
#define DOMAIN_ABSTRACTDOMAIN_H_

#include <deque>
#include <random>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include<vtkSelectEnclosedPoints.h>

#include "../vascularElements/AbstractVascularElement.h"
#include "../CCOCommonStructures.h"
#include "../../core/GeneratorData.h"
#include "IDomainObservable.h"

/**
 * Abstract domain class. Several domain methods guarantee
 * the eficiency during the CCO generation, e.g., getLocalNeighborhood and getDLim,
 * which models the perfusion area of certain terminal in order to establish the minimum
 * length criterion for the new vessel and the local set of vessels where it can be attached
 * to the tree.
 */

class AbstractDomain : public IDomainObservable {
protected:
	/**	Vessels with stages contained in such list are allowed to grow. If empty, all vessels are allowed to grow. */
	vector<int> growingStages;
	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;
	/** Quantity of points that have been consumed. */
	long long int pointCounter;
	/** Boolean value that models if the domain is convex or not. For convex domain some optimizations are applied. */
	bool isConvexDomain;
	/**	Vascular volume.	*/
	double volume;
	/**	Minimum bifurcation angle for vessels generated inside this domain */
	double minAngle;
	/** Boolean value that models if the domain present any restriction regarding the opening angle between the bifurcation plane and a new daughter vessel. */
	bool isBifPlaneContrained;
	/**	Minimum angle between the bifurcation plane and a new daughter vessel. */
	double minPlaneAngle;

public:
	/**
	 * Constructor.
	 * @param instanceData Wrapper with parameters associated to a tree generation process.
	 */
	AbstractDomain(GeneratorData *instanceData);
	/**
	 * Constructor.
	 * @param instanceData Wrapper with parameters associated to a tree generation process.
	 * @param growingStages Stages allowed to grow.
	 */
	AbstractDomain(GeneratorData *instanceData, vector<int> growingStages);
	/**
	 * Destructor.
	 */
	virtual ~AbstractDomain();

	/**
	 * Returns if the segment defined by the vertexes @p xs and @p xf is inside the current domain.
	 * @param xs Start point of the segment.
	 * @param xf End point of the segment.
	 * @return 1 if the segment defined by the vertexes @p xs and @p xf is inside the current domain otherwise 0.
	 */
	virtual int isSegmentInside(point xs, point xf) = 0;

	/**
	 * Returns if the vascular element is eligible to be a parent vessel.
	 * @param element Vascular element.
	 * @return 1 if it is eligible.
	 */
	virtual int isValidElement(AbstractVascularElement *element);
	/**
	 * Returns the maximum opening angle that the hosted tree can generate.
	 * @return Maximum bifurcation angle for vessels generated inside this domain.
	 */
	virtual double getMinBifurcationAngle();
	/**
	 * Sets the maximum opening angle that the hosted tree can generate.
	 * @param minAngle Maximum bifurcation angle for vessels generated inside this domain.
	 */
	virtual void setMinBifurcationAngle(double minAngle);
	/**
	 * Estimates a characteristic length for the current domain. This length is useful to estimate the perfusion volume of the domain.
	 * @return Chracteristic length.
	 */
	virtual double getCharacteristicLength() = 0;
	/**
	 * Minimum distance between a new vessel terminal and the tree based on the current terminal perfusion. This is lower bound criteria for the vessel length and also aids
	 * toward a spatially homogeneous terminal distribution.
	 * @param nVessels Quantity of terminals in the current tree.
	 * @param factor Factor used to scale a perfusion measure in order to estimate the DLim value.
	 * @return DLim value.
	 */
	virtual double getDLim(long long int nVessels, double factor) = 0;
	/**
	 * Returns the local neighbors to the @p p point locus. The amount of neighbors is estimated based on the current terminal perfusion.
	 * @param p Central point of the neighborhood.
	 * @param nVessels Amount of terminals in the tree.
	 * @return Array of neighbor vessels.
	 */
	virtual double *getLocalNeighborhood(point p, long long int nVessels) = 0;
	/**
	 * Computes the size of the domain.
	 * @return Size of the domain.
	 */
	virtual double getSize() = 0;
	/**
	 * Return a random point inside the domain.
	 * @return Point inside the domain.
	 */
	virtual point getRandomPoint() = 0;
	/**
	 * Return a set of inner domain points.
	 * @return Set of inner domain points.
	 */
	virtual deque<point>& getRandomInnerPoints() = 0;
	/**
	 * Returns the vtkPolydata with the domain representation.
	 * @return vtkPolydata with the domain representation.
	 */
	virtual vtkSmartPointer<vtkPolyData>& getVtkGeometry() = 0;
	/**
	 * Returns the quantity of points that have been consumed.
	 * @return Quantity of points consumed.
	 */
	long long int getPointCounter() const;
	/**
	 * Getter for @p isConvexDomain variable.
	 * @return Returns if the domain is setted as convex.
	 */
	bool isIsConvexDomain() const;
	/**
	 * Setter for @p isConvexDomain variable.
	 * @param isConvexDomain Sets the domain as convex or non-convex.
	 */
	void setIsConvexDomain(bool isConvexDomain);
	/**
	 * Default getter of @p instanceData.
	 * @return @p instanceData
	 */
	GeneratorData* getInstanceData();
	/**
	 * Default setter of @p instanceData.
	 * @param instanceData Instance for @p instanceData.
	 */
	void setInstanceData(GeneratorData* instanceData);
	/**
	 * Virtual method to be used by StagedFRROTreeGeneratorLogger.
	 */

	virtual vtkSmartPointer<vtkSelectEnclosedPoints> getEnclosedPoints() = 0;

	virtual void logDomainFiles(FILE *fp) = 0;
	const vector<int>& getGrowingStages() const;
	void setGrowingStages(const vector<int>& growingStages);
	bool isIsBifPlaneContrained() const;
	void setIsBifPlaneContrained(bool isBifPlaneContrained);
	double getMinPlaneAngle();
	void setMinPlaneAngle(double minPlaneAngle);
	virtual int getDraw() = 0;
	virtual int getSeed() = 0;
};

#endif /* DOMAIN_ABSTRACTDOMAIN_H_ */
