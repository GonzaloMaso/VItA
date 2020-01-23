/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeElementalWriter.h
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VTKOBJECTTREEELEMENTALWRITER_H_
#define VTKOBJECTTREEELEMENTALWRITER_H_

#include "VTKObjectTreeWriter.h"

class VTKObjectTreeElementalWriter: public VTKObjectTreeWriter {
public:
	VTKObjectTreeElementalWriter();
	virtual void write(string filename, AbstractObjectCCOTree *tree);
	virtual ~VTKObjectTreeElementalWriter();
private:
	/**
	 * Creates an ordered vector with all VTK identifiers of the vessels in the current tree.
	 * @return Vector with all VTK identifiers of the vessels in the current tree.
	 */
	vector<SingleVessel *> createVTKIndex(AbstractObjectCCOTree* tree);

};

#endif /* VTKOBJECTTREEELEMENTALWRITER_H_ */
