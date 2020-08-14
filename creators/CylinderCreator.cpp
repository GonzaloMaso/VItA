/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CylinderCreator.cpp
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#include "CylinderCreator.h"

#include <vtkCylinderSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>

CylinderCreator::CylinderCreator(vector<double> center, double radius, double height, int resolution) {
	this->center = center;
	this->radius = radius;
	this->height = height;
	this->resolution = resolution;
}

CylinderCreator::~CylinderCreator() {
}

void CylinderCreator::create(string filename) {
	vtkSmartPointer<vtkCylinderSource> cylinderSource =
			vtkSmartPointer<vtkCylinderSource>::New();
	cylinderSource->SetCenter(center[0], center[1], center[2]);
	cylinderSource->SetHeight(height);
	cylinderSource->SetRadius(radius);
	cylinderSource->SetResolution(resolution);
	cylinderSource->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(cylinderSource->GetOutput());
    triangleFilter->Update();

	//	Write in file
	vtkSmartPointer<vtkPolyDataWriter> writer1 =
			vtkSmartPointer<vtkPolyDataWriter>::New();
	writer1->SetFileName(filename.c_str());
	writer1->SetInputData(triangleFilter->GetOutput());
	writer1->SetFileTypeToASCII();
	writer1->Write();
}
