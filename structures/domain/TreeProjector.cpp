/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeProjector.cpp
 *
 *  Created on: 8/05/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "TreeProjector.h"
#include <vtkPolyDataReader.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>

TreeProjector::TreeProjector(string filename){
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
			vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Generate normals for the geometry
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(vtkGeometry);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	vtkGeometry = normalGenerator->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkCellLocator>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	this->offset = 1E-4;

}

TreeProjector::TreeProjector(string filename, double offset){
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
			vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Generate normals for the geometry
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(vtkGeometry);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	vtkGeometry = normalGenerator->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkCellLocator>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	this->offset = offset;
}

TreeProjector::~TreeProjector(){
	// TODO Auto-generated destructor stub
}

void TreeProjector::projectTerminals(vector<SingleVessel *> vessels){
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		point terminal = (*it)->xDist;
		point projection,normal;
		vtkIdType closeCellId;
		int subId;
		double distance;
		locator->FindClosestPoint(terminal.p,projection.p,closeCellId,subId,distance);

		//	Get normal of the projected element - // TESTME
		vtkSmartPointer<vtkDataArray> cellNormalsRetrieved = vtkGeometry->GetCellData()->GetNormals();
		cellNormalsRetrieved->GetTuple(closeCellId,normal.p);

		point displacement = projection-terminal;
		//	Displacement versor
		displacement = displacement / sqrt(displacement^displacement);
		//	Projection + offset - Inner product checks that the offset is applied towards the geometry interior.
		if((normal^displacement)<0)
			projection = projection + displacement * offset;
		else
			projection = projection - displacement * offset;
		(*it)->xDist = projection;	//	Need to update the VTK segment!
	}
}

void TreeProjector::projectVessel(vector<SingleVessel *> vessels){
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		point distal = (*it)->xDist;
		point proximal = (*it)->xProx;
		point projectionProx, projectionDist, normalProx, normalDist,normal;
		vtkIdType closeCellIdProx;
		vtkIdType closeCellIdDist;
		int subId;
		double distance;
		locator->FindClosestPoint(distal.p,projectionProx.p,closeCellIdProx,subId,distance);
		locator->FindClosestPoint(distal.p,projectionDist.p,closeCellIdDist,subId,distance);

		//	Get normal of the projected element
		vtkSmartPointer<vtkDoubleArray> cellNormalsRetrieved = vtkDoubleArray::SafeDownCast(vtkGeometry->GetCellData()->GetNormals());
		cellNormalsRetrieved->GetTuple(closeCellIdProx,normalProx.p);
		cellNormalsRetrieved->GetTuple(closeCellIdDist,normalDist.p);
		normal = (normalProx + normalDist)/2;

		point displacementDist = projectionDist-distal;
		point displacementProx = projectionProx-proximal;
		point displacement = (displacementDist + displacementProx ) / 2;
		//	Displacement versor
		displacement = displacement / sqrt(displacement^displacement);
		//	Projection + offset - Inner product checks that the offset is applied towards the geometry interior.
		if((normal^displacement)<0){
			projectionDist = projectionDist + displacement * offset;
			projectionProx = projectionProx + displacement * offset;
		}
		else{
			projectionDist = projectionDist - displacement * offset;
			projectionProx = projectionProx - displacement * offset;
		}
		(*it)->xDist = projectionDist;	//	Need to update the VTK segment!
		(*it)->xProx = projectionProx;  //	Need to update the VTK segment!
	}
}
