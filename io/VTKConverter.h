/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * STL2VTPConverter.h
 *
 *  Created on: Mar 9, 2018
 *      Author: gonzalo
 */

#ifndef VTKCONVERTER_H_
#define VTKCONVERTER_H_

#include <iostream>

using namespace std;

/**
 * Converts different mesh formats to VTK using VTK library.
 */
class VTKConverter {

public:
	/**
	 * Converts a mesh in STL format to VTK format.
	 * @param stlFile Input STL file.
	 * @param vtkFile Output VTK file.
	 * @param writeASCIIFormat True for ASCII VTK output, otherwise the output will be a binary VTK.
	 */
	static void STL2VTKConverter(string stlFile, string vtkFile, int writeASCIIFormat);

	/**
	 * Converts a mesh in OBJ format to VTK format.
	 * @param objFile Input OBJ file.
	 * @param vtkFile Output VTK file.
	 * @param writeASCIIFormat True for ASCII VTK output, otherwise the output will be a binary VTK.
	 */
	static void OBJ2VTKConverter(string objFile, string vtkFile, int writeASCIIFormat);
};

#endif /* VTKCONVERTER_H_ */
