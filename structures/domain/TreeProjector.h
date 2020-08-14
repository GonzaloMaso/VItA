/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeProjector.h
 *
 *  Created on: 8/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_TREEPROJECTOR_H_
#define DOMAIN_TREEPROJECTOR_H_

#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellLocator.h>

#include "../vascularElements/SingleVessel.h"
#include "../CCOCommonStructures.h"

using namespace std;

/**
 * Projects a tree over or into a geometry. Mesh describing the geometry should have outgoing normals for the offset functionality.
 */
class TreeProjector {
	/** vtkPolydata description of the domain. */
	vtkSmartPointer<vtkPolyData> vtkGeometry;
	/** Cell locator responsible to determine if a segment is inside the domain. */
	vtkSmartPointer<vtkCellLocator> locator;
	/**	Penetration offset for projected points */
	double offset;
public:
	/**
	 * Creates a projector associated to the geometrical domain described @p filename.
	 * @param filename VTK file which contains the surface of the projection domain.
	 */
	TreeProjector(string filename);
	/**
	 * Creates a projector associated to the geometrical domain described @p filename. All projected points will be pushed inside
	 * @param filename VTK file which contains the surface of the projection domain.
	 * @param offset
	 */
	TreeProjector(string filename, double offset);
	virtual ~TreeProjector();

	/**
	 * Projects the terminals of the vessels list to the closest point within the domain.
	 * @param vessels
	 */
	void projectTerminals(vector<SingleVessel *> vessels);
	void projectVessel(vector<SingleVessel *> vessels);

};

#endif /* DOMAIN_TREEPROJECTOR_H_ */
