/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKConverter.cpp
 *
 *  Created on: Mar 9, 2018
 *      Author: gonzalo
 */

#include "VTKConverter.h"

#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkOBJReader.h>
#include <vtkPolyDataWriter.h>

void VTKConverter::STL2VTKConverter(string stlFile, string vtkFile, int writeASCIIFormat){

	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(stlFile.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(vtkFile.c_str());
	writer->SetInputData(reader->GetOutput());
	if(!writeASCIIFormat)
		writer->SetFileTypeToBinary();
	writer->Write();
}

void VTKConverter::OBJ2VTKConverter(string objFile, string vtkFile, int writeASCIIFormat){

	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader->SetFileName(objFile.c_str());
	reader->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(vtkFile.c_str());
	writer->SetInputData(reader->GetOutput());
	if(!writeASCIIFormat)
		writer->SetFileTypeToBinary();
	writer->Write();
}
