/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * ParallelepipedCreator.cpp
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#include "ParallelepipedCreator.h"

#include <vtkCubeSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>

ParallelepipedCreator::ParallelepipedCreator(double *lb, double *ub) : AbstractCreator(){
	this->lb = lb;
	this->ub = ub;
}

ParallelepipedCreator::~ParallelepipedCreator() {
}

void ParallelepipedCreator::create(string filename) {
	double midPoint[3];
	midPoint[0] = (lb[0] + ub[0])/2.;
	midPoint[1] = (lb[1] + ub[1])/2.;
	midPoint[2] = (lb[2] + ub[2])/2.;

	vtkSmartPointer<vtkCubeSource> cubeSource =
			vtkSmartPointer<vtkCubeSource>::New();
	cubeSource->SetCenter(midPoint[0], midPoint[1], midPoint[2]);
	cubeSource->SetXLength(abs(lb[0] - ub[0]));
	cubeSource->SetYLength(abs(lb[1] - ub[1]));
	cubeSource->SetZLength(abs(lb[2] - ub[2]));
	cubeSource->Update();

	//	Write in file
	vtkSmartPointer<vtkPolyDataWriter> writer1 =
			vtkSmartPointer<vtkPolyDataWriter>::New();
	writer1->SetFileName(filename.c_str());
	writer1->SetInputData(cubeSource->GetOutput());
	writer1->SetFileTypeToASCII();
	writer1->Write();
}
