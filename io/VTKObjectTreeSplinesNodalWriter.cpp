/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VTKObjectTreeSplinesNodalWriter.cpp
 *
 *  Created on: 9/08/2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "VTKObjectTreeSplinesNodalWriter.h"

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

VTKObjectTreeSplinesNodalWriter::VTKObjectTreeSplinesNodalWriter()
{
	this->resolution = 10;
}

VTKObjectTreeSplinesNodalWriter::VTKObjectTreeSplinesNodalWriter(int resolution)
		{
	this->resolution = resolution;
}

void VTKObjectTreeSplinesNodalWriter::write(string filename, AbstractObjectCCOTree* tree) {
	tree->computeTreeCost(tree->getRoot());

	vtkSmartPointer<vtkPolyData> vtkNewTree = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkDoubleArray> nodeDataRadius = vtkSmartPointer<vtkDoubleArray>::New();
	nodeDataRadius->SetName("radius");
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

			point proximalDirection, distalDirection;
			computeDerivatives(currentSegment, proximalDirection, distalDirection);

			/* Splines resampling*/
			vtkSmartPointer<vtkPoints> localPoints = vtkSmartPointer<vtkPoints>::New();
			vtkIdType idProx = localPoints->InsertNextPoint(currentSegment->xProx.p);
			vtkIdType idDist = localPoints->InsertNextPoint(currentSegment->xDist.p);

			vtkSmartPointer<vtkCardinalSpline> splineX =
					vtkSmartPointer<vtkCardinalSpline>::New();
			splineX->ClosedOff();
			splineX->SetLeftConstraint(1);
			splineX->SetLeftValue(proximalDirection.p[0]);
			splineX->SetRightConstraint(1);
			splineX->SetRightValue(distalDirection.p[0]);
			splineX->AddPoint(0, currentSegment->xProx.p[0]);
			splineX->AddPoint(1, currentSegment->xDist.p[0]);

			vtkSmartPointer<vtkCardinalSpline> splineY =
					vtkSmartPointer<vtkCardinalSpline>::New();
			splineY->ClosedOff();
			splineY->SetLeftConstraint(1);
			splineY->SetLeftValue(proximalDirection.p[1]);
			splineY->SetRightConstraint(1);
			splineY->SetRightValue(distalDirection.p[1]);
			splineY->AddPoint(0, currentSegment->xProx.p[1]);
			splineY->AddPoint(1, currentSegment->xDist.p[1]);

			vtkSmartPointer<vtkCardinalSpline> splineZ =
					vtkSmartPointer<vtkCardinalSpline>::New();
			splineZ->ClosedOff();
			splineZ->SetLeftConstraint(1);
			splineZ->SetLeftValue(proximalDirection.p[2]);
			splineZ->SetRightConstraint(1);
			splineZ->SetRightValue(distalDirection.p[2]);
			splineZ->AddPoint(0, currentSegment->xProx.p[2]);
			splineZ->AddPoint(1, currentSegment->xDist.p[2]);

			point previous = currentSegment->xProx;
			point next;
			double increment = 1. / (double) this->resolution;
			for (double t = increment; t <= 1 + 0.1 * increment; t += increment) {
				next.p[0] = splineX->Evaluate(t);
				next.p[1] = splineY->Evaluate(t);
				next.p[2] = splineZ->Evaluate(t);
				vtkIdType idProx = pts->InsertNextPoint(previous.p);
				nodeDataRadius->InsertNextValue(currentSegment->radius);
				vtkIdType idDist = pts->InsertNextPoint(next.p);
				nodeDataRadius->InsertNextValue(currentSegment->radius);

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

				previous = next;
			}
		}

	}

	vtkNewTree->SetPoints(pts);
	vtkNewTree->SetLines(lines);

	vtkNewTree->GetPointData()->AddArray(nodeDataRadius);

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

VTKObjectTreeSplinesNodalWriter::~VTKObjectTreeSplinesNodalWriter()
{
	// TODO Auto-generated destructor stub
}

void VTKObjectTreeSplinesNodalWriter::computeDerivatives(SingleVessel* vessel, point& proximal, point& distal) {
	SingleVessel *parent = (SingleVessel *) vessel->parent;

	//	Proximal derivative considering the ratio of the radii.
	if (parent) {
		point proximalParent = parent->xDist - parent->xProx;
		double rP = parent->radius;
		point proximalVessel = vessel->xDist - vessel->xProx;
		double rV = vessel->radius;
		proximal = proximalParent * (rP/(rP+rV)) + proximalVessel * (rV/(rP+rV));
	}
	else {
		proximal = vessel->xDist - vessel->xProx;
	}

	distal = vessel->xDist - vessel->xProx;
}
