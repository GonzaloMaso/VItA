/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeNodalWriter.cpp
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "VTKObjectTreeNodalWriter.h"

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

VTKObjectTreeNodalWriter::VTKObjectTreeNodalWriter()
{
	// TODO Auto-generated constructor stub

}

void VTKObjectTreeNodalWriter::write(string filename, AbstractObjectCCOTree* tree) {
	tree->computeTreeCost(tree->getRoot());

	vtkSmartPointer<vtkPolyData> vtkNewTree = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkDoubleArray> nodeDataRadius = vtkSmartPointer<vtkDoubleArray>::New();
	nodeDataRadius->SetName("radius");
	vtkSmartPointer<vtkDoubleArray> nodeDataFlow = vtkSmartPointer<vtkDoubleArray>::New();
	nodeDataFlow->SetName("flow");
	vtkSmartPointer<vtkDoubleArray> nodeDataDistal = vtkSmartPointer<vtkDoubleArray>::New();
	nodeDataDistal->SetName("isDistal");
	vtkSmartPointer<vtkDoubleArray> nodeDataTerminal = vtkSmartPointer<vtkDoubleArray>::New();
	nodeDataTerminal->SetName("isTerminal");
	vtkSmartPointer<vtkDoubleArray> cellDataLevel = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLevel->SetName("level");
	vtkSmartPointer<vtkDoubleArray> cellDataFlow = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataFlow->SetName("flow");
	vtkSmartPointer<vtkDoubleArray> cellDataPressure = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataPressure->SetName("pressure");
	vtkSmartPointer<vtkDoubleArray> cellDataResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataResistance->SetName("resistance");
	vtkSmartPointer<vtkDoubleArray> cellDataViscosity = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataViscosity->SetName("viscosity");
	vtkSmartPointer<vtkDoubleArray> cellDataEqResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataEqResistance->SetName("eqResistance");
	vtkSmartPointer<vtkDoubleArray> cellDataID = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataID->SetName("ID_Vessel");
	vtkSmartPointer<vtkDoubleArray> cellDataRadius = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataRadius->SetName("radius");
	vtkSmartPointer<vtkDoubleArray> cellDataLength = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLength->SetName("length");
	vtkSmartPointer<vtkDoubleArray> cellDataBeta = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataBeta->SetName("beta");
	vtkSmartPointer<vtkDoubleArray> cellDataGeoConst = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataGeoConst->SetName("geomConst");
	vtkSmartPointer<vtkDoubleArray> cellDataStage = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataStage->SetName("stage");
	vtkSmartPointer<vtkDoubleArray> cellDataBranch = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataBranch->SetName("branchingType");
	vtkSmartPointer<vtkDoubleArray> cellDataTerminalType = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataTerminalType->SetName("terminalType");
	vtkSmartPointer<vtkDoubleArray> cellDataVesselFunction = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataVesselFunction->SetName("vesselFunction");

	for (auto it = tree->getSegments().begin(); it != tree->getSegments().end(); ++it) {
		for (vector<SingleVessel *>::iterator it2 = (it->second)->getVessels().begin(); it2 != (it->second)->getVessels().end(); ++it2) {
			SingleVessel *currentSegment = *it2;

			vtkIdType idProx = pts->InsertNextPoint(currentSegment->xProx.p);
			nodeDataRadius->InsertNextValue(currentSegment->radius);
			nodeDataFlow->InsertNextValue(currentSegment->flow);
			nodeDataDistal->InsertNextValue(0.0);
			nodeDataTerminal->InsertNextValue(0.0);
			vtkIdType idDist = pts->InsertNextPoint(currentSegment->xDist.p);
			nodeDataRadius->InsertNextValue(currentSegment->radius);
			nodeDataFlow->InsertNextValue(currentSegment->flow);
			nodeDataDistal->InsertNextValue(1.0);
			if(currentSegment->getChildren().size() == 0){
				nodeDataTerminal->InsertNextValue(1.0);
			}else{
				nodeDataTerminal->InsertNextValue(0.0);
			}

			vtkSmartPointer<vtkLine> vtkSegment = vtkSmartPointer<vtkLine>::New();
			vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
			vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
			lines->InsertNextCell(vtkSegment);

			cellDataLevel->InsertNextValue(currentSegment->nLevel);
			cellDataFlow->InsertNextValue(currentSegment->flow);
			cellDataPressure->InsertNextValue(currentSegment->pressure);
			cellDataResistance->InsertNextValue(currentSegment->localResistance);
			cellDataViscosity->InsertNextValue(currentSegment->viscosity);
			cellDataEqResistance->InsertNextValue(currentSegment->resistance);
			cellDataID->InsertNextValue(currentSegment->ID);
			cellDataRadius->InsertNextValue(currentSegment->radius);
			cellDataLength->InsertNextValue(currentSegment->length);
			cellDataBeta->InsertNextValue(currentSegment->beta);
			cellDataGeoConst->InsertNextValue(currentSegment->length - 2 * currentSegment->radius);
			cellDataStage->InsertNextValue(currentSegment->stage);
			cellDataBranch->InsertNextValue(currentSegment->branchingMode);
			cellDataTerminalType->InsertNextValue(currentSegment->terminalType);
			cellDataVesselFunction->InsertNextValue(currentSegment->vesselFunction);

		}

	}

	vtkNewTree->SetPoints(pts);
	vtkNewTree->SetLines(lines);

	vtkNewTree->GetPointData()->AddArray(nodeDataRadius);
	vtkNewTree->GetPointData()->AddArray(nodeDataFlow);
	vtkNewTree->GetPointData()->AddArray(nodeDataDistal);
	vtkNewTree->GetPointData()->AddArray(nodeDataTerminal);

	vtkNewTree->GetCellData()->AddArray(cellDataLevel);
	vtkNewTree->GetCellData()->AddArray(cellDataFlow);
	vtkNewTree->GetCellData()->AddArray(cellDataPressure);
	vtkNewTree->GetCellData()->AddArray(cellDataResistance);
	vtkNewTree->GetCellData()->AddArray(cellDataViscosity);
	vtkNewTree->GetCellData()->AddArray(cellDataEqResistance);
	vtkNewTree->GetCellData()->AddArray(cellDataID);
	vtkNewTree->GetCellData()->AddArray(cellDataRadius);
	vtkNewTree->GetCellData()->AddArray(cellDataLength);
	vtkNewTree->GetCellData()->AddArray(cellDataBeta);
	vtkNewTree->GetCellData()->AddArray(cellDataGeoConst);
	vtkNewTree->GetCellData()->AddArray(cellDataStage);
	vtkNewTree->GetCellData()->AddArray(cellDataBranch);
	vtkNewTree->GetCellData()->AddArray(cellDataTerminalType);
	vtkNewTree->GetCellData()->AddArray(cellDataVesselFunction);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(vtkNewTree);
	writer->SetDataModeToBinary();
	writer->Write();
}

VTKObjectTreeNodalWriter::~VTKObjectTreeNodalWriter()
{
	// TODO Auto-generated destructor stub
}

