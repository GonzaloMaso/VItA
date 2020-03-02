/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeSplinesNodalWriter.h
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VTKOBJECTTREESPLINESNODALWRITER_H_
#define VTKOBJECTTREESPLINESNODALWRITER_H_

#include "VTKObjectTreeWriter.h"

class VTKObjectTreeSplinesNodalWriter: public VTKObjectTreeWriter {
	int resolution;
public:
	VTKObjectTreeSplinesNodalWriter();
	VTKObjectTreeSplinesNodalWriter(int resolution);
	virtual void write(string filename, AbstractObjectCCOTree *tree);
	virtual ~VTKObjectTreeSplinesNodalWriter();
private:
	void computeDerivatives(SingleVessel *vessel, point &proximal, point &distal);
};

#endif /* VTKOBJECTTREESPLINESNODALWRITER_H_ */
