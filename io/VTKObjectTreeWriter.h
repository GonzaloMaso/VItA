/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKTreeWriter.h
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VTKOBJECTTREEWRITER_H_
#define VTKOBJECTTREEWRITER_H_

#include <iostream>
#include "../structures/tree/AbstractObjectCCOTree.h"

using namespace std;

class VTKObjectTreeWriter {
public:
	VTKObjectTreeWriter();
	virtual void write(string filename, AbstractObjectCCOTree *tree) = 0;
	virtual ~VTKObjectTreeWriter();
};

#endif /* VTKOBJECTTREEWRITER_H_ */
