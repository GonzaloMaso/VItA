/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StagedDomain.h
 *
 *  Created on: Mar 14, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_STAGEDDOMAIN_H_
#define DOMAIN_STAGEDDOMAIN_H_

#include "AbstractDomain.h"

/**
 * Composite class that modify the domain depending in the amount of terminals that
 * have been included.
 */
class StagedDomain : public AbstractDomain {
	/** Domain at each stage. */
	vector< AbstractDomain *> domainStage;
	/** Terminals generated at each stage. */
	vector<long long int> terminalsPerStage;
	/** Terminals counter. */
	long long int currentTerminals;
	/**	Terminals at the end of the previous stage. */
	long long int terminalAtPrevStage;
	/**	Index of the current stage. */
	int currentStage;

	int initialStage;
public:
	/**
	 * Constructor.
	 */
	StagedDomain();
	/**
	 * Descturctor.
	 */
	~StagedDomain();
	/**
	 * Add an additional stage to the domain.
	 * @param terminals Amount of terminals for the added stage.
	 * @param domain Domain used for the added stage.
	 */
	void addStage(long long int terminals, AbstractDomain *domain);
	/**
	 * Increase the amount of @p currentTerminals by 1. It must be called after the addition
	 * of a new terminal to the tree.
	 */
	void update();
	/**
	 * Returns the total amount of terminals that will be generated along all stages.
	 * @return Total amount of terminals generated along all stages.
	 */
	long long int getTerminalsAfterGeneration();
	/**
	 * Returns if the segment defined by the vertexes @p xs and @p xf is inside the current domain.
	 * @param xs Start point of the segment.
	 * @param xf End point of the segment.
	 * @return 1 if the segment defined by the vertexes @p xs and @p xf is inside the current domain otherwise 0.
	 */
	int isSegmentInside(point xs, point xf);
	/**
	 * Returns if the vascular element is eligible to be a parent vessel.
	 * @param element Vascular element.
	 * @return 1 if it is eligible.
	 */
	int isValidElement(AbstractVascularElement *element);
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
	 * Returns the maximum angle that the hosted tree can generate.
	 * @return Maximum bifurcation angle for vessels generated inside this domain.
	 */
	double getMinBifurcationAngle();
	/**
	 * Returns the maximum opening angle that the hosted tree can generate.
	 * @return Maximum opening angle for vessels generated inside this domain.
	 */
	double getMinPlaneAngle();
	/**
	 * Computes the size of the domain.
	 * @return Size of the domain.
	 */
	double getSize();
	/**
	 * Return a random point inside the domain.
	 * @return Point inside the domain.
	 */
	point getRandomPoint();
	/**
	 * Return a set of inner domain points.
	 * @return Set of inner domain points.
	 */
	deque<point>& getRandomInnerPoints();
	/**
	 * Returns the quantity of points that have been consumed.
	 * @return Quantity of points consumed.
	 */
	long long int getPointCounter() const;
	/**
	 * Returns the vtkPolydata with the domain representation.
	 * @return vtkPolydata with the domain representation.
	 */
	vtkSmartPointer<vtkPolyData>& getVtkGeometry();
	/**
	 * Returns the current stage.
	 * @return Current stage number;
	 */
	int getCurrentStage() const;
	/**
	 * Sets the initial stage as @p currentStage.
	 * @param currentStage Current stage.
	 */
	void setInitialStage(int currentStage);
	/**
	 * Returns the domains in the StagedDomain.
	 * @return vector<AbstractDomain *> domains.
	 */
	vector< AbstractDomain *>* getDomains();
	/**
	 * Returns the number of terminals in each stage.
	 * @return  vector<long long int> n_term.
	 */
	vector<long long int>* getNTerminals();
	/**
	 * This function serves no purpose in StagedDomain.
	 */
	int getDraw();
	/**
	 * This function serves no purpose in StagedDomain.
	 */
	int getSeed();

	void logDomainFiles(FILE *fp);

	vtkSmartPointer<vtkSelectEnclosedPoints> getEnclosedPoints() override;
};

#endif /* DOMAIN_STAGEDDOMAIN_H_ */
