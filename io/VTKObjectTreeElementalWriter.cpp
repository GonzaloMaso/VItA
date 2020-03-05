/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeElementalWriter.cpp
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "VTKObjectTreeElementalWriter.h"

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkLine.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkParametricSpline.h>
#include <vtkCardinalSpline.h>
#include <vtkSpline.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

#include <vector>

#include "../structures/vascularElements/SingleVessel.h"

VTKObjectTreeElementalWriter::VTKObjectTreeElementalWriter()
{
	// TODO Auto-generated constructor stub

}

/*
 * This will not work when vtkIndexes are not consecutive and complete. Maybe use a map ordered.
 */
vector<SingleVessel*> VTKObjectTreeElementalWriter::createVTKIndex(AbstractObjectCCOTree* tree) {
	unsigned long long int amount = 0;
	for (auto it = tree->getSegments().begin(); it != tree->getSegments().end(); ++it) {
		amount += (it->second)->getVessels().size();
	}
	vector<SingleVessel *> index(amount, NULL);
	vector<SingleVessel *> vessels = tree->getVessels();
	for (auto it = vessels.begin(); it != vessels.end(); ++it) {
		index[(*it)->vtkSegmentId] = (*it);
	}
	return index;
}

void VTKObjectTreeElementalWriter::write(string filename, AbstractObjectCCOTree* tree){
	vector<SingleVessel *> vtkIndex = createVTKIndex(tree);
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());

	vtkSmartPointer<vtkDoubleArray> cellDataLevel = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLevel->SetName("level");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLevel->InsertNextValue((*it)->nLevel);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataLevel);

	vtkSmartPointer<vtkDoubleArray> cellDataFlux = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataFlux->SetName("flow");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataFlux->InsertNextValue((*it)->flow);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataFlux);

	vtkSmartPointer<vtkDoubleArray> cellDataPressure = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataPressure->SetName("pressure");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataPressure->InsertNextValue((*it)->pressure);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataPressure);

	vtkSmartPointer<vtkDoubleArray> cellDataLResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLResistance->SetName("resistance");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLResistance->InsertNextValue((*it)->localResistance);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataLResistance);

	vtkSmartPointer<vtkDoubleArray> cellDataViscosity = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataViscosity->SetName("viscosity");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataViscosity->InsertNextValue((*it)->viscosity);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataViscosity);

	vtkSmartPointer<vtkDoubleArray> cellDataResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataResistance->SetName("eqResistance");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataResistance->InsertNextValue((*it)->resistance);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataResistance);

	vtkSmartPointer<vtkDoubleArray> cellDataID = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataID->SetName("ID_Vessel");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataID->InsertNextValue((*it)->ID);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataID);

	tree->computeTreeCost(tree->getRoot());
	vtkSmartPointer<vtkDoubleArray> cellDataRadius = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataRadius->SetName("radius");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataRadius->InsertNextValue((*it)->radius);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataRadius);

	vtkSmartPointer<vtkDoubleArray> cellDataLength = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLength->SetName("length");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLength->InsertNextValue((*it)->length);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataLength);

	vtkSmartPointer<vtkDoubleArray> cellDataBeta = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataBeta->SetName("beta");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataBeta->InsertNextValue((*it)->beta);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataBeta);

	vtkSmartPointer<vtkDoubleArray> cellDataGeomConst = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataGeomConst->SetName("geomConst");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataGeomConst->InsertNextValue((*it)->length - 2 * (*it)->radius);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataGeomConst);

	vtkSmartPointer<vtkDoubleArray> cellDataStage = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataStage->SetName("Stage");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataStage->InsertNextValue((*it)->stage);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataStage);

	vtkSmartPointer<vtkDoubleArray> cellDataBranch = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataBranch->SetName("branchingMode");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataBranch->InsertNextValue((int) (*it)->branchingMode);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataBranch);

	vtkSmartPointer<vtkDoubleArray> cellDataTermType = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataTermType->SetName("terminalType");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataTermType->InsertNextValue((int) (*it)->terminalType);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataTermType);

	vtkSmartPointer<vtkDoubleArray> cellDataVesselFunction = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataVesselFunction->SetName("vesselFunction");
	for (vector<SingleVessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataVesselFunction->InsertNextValue((int) (*it)->vesselFunction);
	}
	tree->getVtkTree()->GetCellData()->AddArray(cellDataVesselFunction);

	writer->SetInputData(tree->getVtkTree());
	writer->SetDataModeToBinary();
	writer->Write();
}

VTKObjectTreeElementalWriter::~VTKObjectTreeElementalWriter()
{
	// TODO Auto-generated destructor stub
}

