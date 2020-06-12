/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SphereCreator.cpp
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#include "SphereCreator.h"

#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>

SphereCreator::SphereCreator(vector<double> center, double radius, double phiResolution, double thetaResolution){
	this->center = center;
	this->radius = radius;
	this->phiResolution = phiResolution;
	this->thetaResolution = thetaResolution;
}

SphereCreator::~SphereCreator() {
}

void SphereCreator::create(string filename) {
	vtkSmartPointer<vtkSphereSource> sphereSource1 =
			vtkSmartPointer<vtkSphereSource>::New();
	sphereSource1->SetCenter(center[0], center[1], center[2]);
	sphereSource1->SetRadius(radius);
	sphereSource1->SetPhiResolution(phiResolution);
	sphereSource1->SetThetaResolution(thetaResolution);
	sphereSource1->Update();
	//	Write in file
	vtkSmartPointer<vtkPolyDataWriter> writer1 =
			vtkSmartPointer<vtkPolyDataWriter>::New();
	writer1->SetFileName(filename.c_str());
	writer1->SetInputData(sphereSource1->GetOutput());
	writer1->SetFileTypeToASCII();
	writer1->Write();
}

