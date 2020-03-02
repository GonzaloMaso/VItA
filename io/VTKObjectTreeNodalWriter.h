/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeNodalWriter.h
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VTKOBJECTTREENODALWRITER_H_
#define VTKOBJECTTREENODALWRITER_H_

#include "VTKObjectTreeWriter.h"

class VTKObjectTreeNodalWriter: public VTKObjectTreeWriter {
public:
	VTKObjectTreeNodalWriter();
	virtual void write(string filename, AbstractObjectCCOTree *tree);
	virtual ~VTKObjectTreeNodalWriter();
};

#endif /* VTKOBJECTTREENODALWRITER_H_ */
