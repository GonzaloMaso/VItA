/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractCCOTree.cpp
 *
 *  Created on: Feb 5, 2018
 *      Author: gonzalo
 */

#include "AbstractStructuredCCOTree.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

AbstractStructuredCCOTree::AbstractStructuredCCOTree(GeneratorData *instanceData) {

	this->instanceData = instanceData;

	this->pointCounter = 0;
	this->refPressure = 0.0;

	this->vtkTree = vtkSmartPointer<vtkPolyData>::New();
	this->vtkTreeLocator = vtkSmartPointer<vtkCellLocator>::New();
}

AbstractStructuredCCOTree::AbstractStructuredCCOTree(point xi, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu, double minAngle, double refPressure, GeneratorData *instanceData) {

	this->instanceData = instanceData;

	this->pointCounter = 0;

	this->xPerf = xi;
	this->qProx = qi;

	this->gam = gam;
	this->epsLim = epsLim;
	this->nu = nu;
	this->minAngle = minAngle;
	this->refPressure = refPressure;

	this->psiFactor = 0.0;
	this->dp = 0.0;
	this->nTerms = 0;

	this->root = NULL;

	this->vtkTree = vtkSmartPointer<vtkPolyData>::New();
	this->vtkTreeLocator = vtkSmartPointer<vtkCellLocator>::New();
}

AbstractStructuredCCOTree::~AbstractStructuredCCOTree() {

	for (unsigned long long int i = 0; i < segments.size(); ++i) {
		delete[] segments[i];
	}
	segments.clear();
}

AbstractConstraintFunction<double, int> *AbstractStructuredCCOTree::getEpsLim() {
	return epsLim;
}

AbstractConstraintFunction<double, int> *AbstractStructuredCCOTree::getGam() {
	return gam;
}

AbstractConstraintFunction<double, int> *AbstractStructuredCCOTree::getNu() {
	return nu;
}

double AbstractStructuredCCOTree::getQProx() {
	return qProx;
}

vessel *AbstractStructuredCCOTree::getRoot() {
	return root;
}

const vector<vessel*>& AbstractStructuredCCOTree::getSegments() const {
	return segments;
}

vtkSmartPointer<vtkPolyData> AbstractStructuredCCOTree::getVtkTree() {
	return vtkTree;
}

point AbstractStructuredCCOTree::getXProx() {
	return xPerf;
}

void AbstractStructuredCCOTree::save(string filename) {
	ofstream treeFile;

	treeFile.open(filename.c_str(), ios::out);
	treeFile.setf(ios::scientific, ios::floatfield);
	treeFile.precision(16);

	treeFile << "*Tree" << endl;
	saveTree(&treeFile);
	treeFile << endl << endl;

	treeFile << "*Vessels" << endl << segments.size() << endl;

	for (unsigned int i = 0; i < segments.size(); ++i) {
		vessel *currentVessel = segments[i];
		saveVessel(currentVessel, &treeFile);
		treeFile << endl;
	}
	treeFile << endl;

	treeFile << "*Connectivity" << endl;
	for (unsigned int i = 0; i < segments.size(); ++i) {
		vessel *currentVessel = segments[i];
		treeFile << currentVessel->vtkSegmentId << " ";
		for (unsigned int j = 0; j < currentVessel->anastomose.size(); ++j) {
			vessel *currentAnastomosis = currentVessel->anastomose[j];
			if (currentAnastomosis)
				treeFile << currentAnastomosis->vtkSegmentId << " ";
			else
				treeFile << -1 << " ";
		}
		treeFile << endl;
	}

	treeFile.flush();
	treeFile.close();
}

void AbstractStructuredCCOTree::saveVessel(vessel* v, ofstream *outFile) {
	*outFile << v->vtkSegmentId << " " << v->xProx.p[0] << " " << v->xProx.p[1] << " " << v->xProx.p[2] << " " << v->xDist.p[0] << " " << v->xDist.p[1] << " " << v->xDist.p[2] << " " << v->nLevel
			<< " " << v->radius << " " << v->beta << " " << v->length << " " << v->resistance << " " << v->flux << " " << v->pressure << " " << v->ID << " " << v->treeVolume;
}

void AbstractStructuredCCOTree::saveTree(ofstream* outFile) {

	*outFile << xPerf.p[0] << " " << xPerf.p[1] << " " << xPerf.p[2] << " " << qProx << " " << minAngle << " " << psiFactor << " " << dp << " " << nTerms << " " << refPressure << " "
			<< (long long) pointCounter << " ";

}

double AbstractStructuredCCOTree::getDp() const {
	return dp;
}

double AbstractStructuredCCOTree::getMinAngle() const {
	return minAngle;
}

void AbstractStructuredCCOTree::generateVTKstructures(vessel* root) {
	if (root) {
		cout << "Regenerating vessel #" << root->vtkSegmentId << endl;
		if (!root->anastomose[0]) {
			//	Update tree geometry
			vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
			vtkIdType idProx = pts->InsertNextPoint(root->xProx.p);
			vtkIdType idDist = pts->InsertNextPoint(root->xDist.p);
			vtkTree->SetPoints(pts);

			root->vtkSegment = vtkSmartPointer<vtkLine>::New();
			root->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
			root->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
			vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
			root->vtkSegmentId = lines->InsertNextCell(root->vtkSegment);
			cout << "New ID #" << root->vtkSegmentId << endl;
			vtkTree->SetLines(lines);

			//	Update tree locator
			vtkTreeLocator->SetDataSet(vtkTree);
			vtkTreeLocator->BuildLocator();
		} else {
			vessel *parent = root->anastomose[0];
			//	Update tree geometry
			vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(root->xDist.p);

			root->vtkSegment = vtkSmartPointer<vtkLine>::New();
			root->vtkSegment->GetPointIds()->SetId(0, parent->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
			root->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

			root->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(root->vtkSegment);
			cout << "New ID #" << root->vtkSegmentId << endl;

			vtkTree->BuildCells();
			vtkTree->Modified();

//			cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
//			cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

			vtkTreeLocator->Update();
		}
		if (root->anastomose.size() > 1) {
			generateVTKstructures(root->anastomose[1]);
			generateVTKstructures(root->anastomose[2]);
		}
	}
}

vector<vessel*> AbstractStructuredCCOTree::createVTKIndex() {
	vector<vessel *> index(segments.begin(), segments.end());
	for (vector<vessel *>::iterator it = segments.begin(); it != segments.end(); ++it) {
		index[(*it)->vtkSegmentId] = *it;
	}
	return index;
}

void AbstractStructuredCCOTree::setEpsLim(AbstractConstraintFunction<double, int> *epsLim) {
	this->epsLim = epsLim;
}

void AbstractStructuredCCOTree::setGam(AbstractConstraintFunction<double, int> *gam) {
	this->gam = gam;
}

void AbstractStructuredCCOTree::setMinAngle(double minAngle) {
	this->minAngle = minAngle;
}

void AbstractStructuredCCOTree::setNu(AbstractConstraintFunction<double, int> *nu) {
	this->nu = nu;
}

long long int AbstractStructuredCCOTree::getPointCounter() const {
	return pointCounter;
}

void AbstractStructuredCCOTree::setPointCounter(long long int pointCounter) {
	this->pointCounter = pointCounter;
}

long long int AbstractStructuredCCOTree::getNTerminals() {
	return countTerminals(this->root);
}

long long int AbstractStructuredCCOTree::countTerminals(vessel* root) {
	if (root->anastomose.size() > 1) {
		return countTerminals(root->anastomose[1]) + countTerminals(root->anastomose[2]);
	} else {
		return 1;
	}
}

void AbstractStructuredCCOTree::storeVTK(string filename) {
	vector<vessel *> vtkIndex = createVTKIndex();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());

	vtkSmartPointer<vtkDoubleArray> cellDataLevel = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLevel->SetName("level");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLevel->InsertNextValue((*it)->nLevel);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataLevel);

	vtkSmartPointer<vtkDoubleArray> cellDataFlux = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataFlux->SetName("flow");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataFlux->InsertNextValue((*it)->flux);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataFlux);

	vtkSmartPointer<vtkDoubleArray> cellDataPressure = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataPressure->SetName("pressure");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataPressure->InsertNextValue((*it)->pressure);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataPressure);

	vtkSmartPointer<vtkDoubleArray> cellDataLResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLResistance->SetName("resistance");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLResistance->InsertNextValue((*it)->localResistance);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataLResistance);

	vtkSmartPointer<vtkDoubleArray> cellDataViscosity = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataViscosity->SetName("viscosity");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataViscosity->InsertNextValue((*it)->viscosity);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataViscosity);

	vtkSmartPointer<vtkDoubleArray> cellDataResistance = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataResistance->SetName("eqResistance");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataResistance->InsertNextValue((*it)->resistance);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataResistance);

	vtkSmartPointer<vtkDoubleArray> cellDataID = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataID->SetName("ID_Vessel");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataID->InsertNextValue((*it)->ID);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataID);

	computeTreeCost(root, 1);
	vtkSmartPointer<vtkDoubleArray> cellDataRadius = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataRadius->SetName("radius");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataRadius->InsertNextValue((*it)->radius);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataRadius);

	vtkSmartPointer<vtkDoubleArray> cellDataLength = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataLength->SetName("length");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataLength->InsertNextValue((*it)->length);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataLength);

	vtkSmartPointer<vtkDoubleArray> cellDataBeta = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataBeta->SetName("beta");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataBeta->InsertNextValue((*it)->beta);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataBeta);

	vtkSmartPointer<vtkDoubleArray> cellDataGeomConst = vtkSmartPointer<vtkDoubleArray>::New();
	cellDataGeomConst->SetName("geomConst");
	for (vector<vessel *>::iterator it = vtkIndex.begin(); it != vtkIndex.end(); ++it) {
		cellDataGeomConst->InsertNextValue((*it)->length - 2 * (*it)->radius);
	}
	this->getVtkTree()->GetCellData()->AddArray(cellDataGeomConst);

	writer->SetInputData(this->getVtkTree());
	writer->SetDataModeToBinary();
	writer->Write();
}

void AbstractStructuredCCOTree::storeVTK(string filename, unsigned int mode){
	/*	Updates radius in the whole structure */
	computeTreeCost(root, 1);

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

	for(std::vector<vessel *>::iterator it = segments.begin(); it != segments.end(); ++it) {
		vessel *currentSegment = *it;

		vtkIdType idProx = pts->InsertNextPoint(currentSegment->xProx.p);
		nodeDataRadius->InsertNextValue(currentSegment->radius);
		vtkIdType idDist = pts->InsertNextPoint(currentSegment->xDist.p);
		nodeDataRadius->InsertNextValue(currentSegment->radius);

		vtkSmartPointer<vtkLine> vtkSegment = vtkSmartPointer<vtkLine>::New();
		vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		lines->InsertNextCell(vtkSegment);
		cellDataLevel->InsertNextValue(currentSegment->nLevel);
		cellDataFlow->InsertNextValue(currentSegment->flux);
		cellDataPressure->InsertNextValue(currentSegment->pressure);
		cellDataResistance->InsertNextValue(currentSegment->localResistance);
		cellDataViscosity->InsertNextValue(currentSegment->viscosity);
		cellDataEqResistance->InsertNextValue(currentSegment->resistance);
		cellDataID->InsertNextValue(currentSegment->ID);
		cellDataRadius->InsertNextValue(currentSegment->radius);
		cellDataLength->InsertNextValue(currentSegment->length);
		cellDataBeta->InsertNextValue(currentSegment->beta);
		cellDataGeoConst->InsertNextValue(currentSegment->length - 2 * currentSegment->radius);
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

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(vtkNewTree);
	writer->SetDataModeToBinary();
	writer->Write();

}


double AbstractStructuredCCOTree::computeTreeCost(vessel* root, double parentRadius) {

	root->radius = root->beta * parentRadius;
	//	Geometric constraint
	if (root->anastomose.size() > 1) {
		return M_PI * root->radius * root->radius * root->length + computeTreeCost(root->anastomose[1], root->radius) + computeTreeCost(root->anastomose[2], root->radius);
//		return root->length + computeTreeCost(root->anastomose[1], root->radius) + computeTreeCost(root->anastomose[2], root->radius);
	} else {
		if (root->radius * 4 < root->length) {
			return M_PI * root->radius * root->radius * root->length;
//			return root->length;
		} else {
//			cout << "Segment " << root->vtkSegmentId << " violate geometric constraint" << endl;
			return INFINITY;
		}
	}

}

void AbstractStructuredCCOTree::computePressure(vessel* root) {
	if (root->anastomose.size() > 1) {
		computePressure(root->anastomose[1]);
		computePressure(root->anastomose[2]);
		root->localResistance = 8 * root->viscosity / M_PI * root->length;
		root->pressure = root->flux * root->localResistance + root->anastomose[1]->pressure + root->anastomose[2]->pressure;
	} else {
		root->localResistance = 8 * root->viscosity / M_PI * root->length;
		root->pressure = root->flux * root->localResistance;
	}
}

void AbstractStructuredCCOTree::updateSegmentVtkLines() {
	vtkTree->GetLines()->InitTraversal();
	for (std::vector<vessel *>::iterator it = segments.begin(); it != segments.end(); ++it) {
		(*it)->vtkSegment = vtkSmartPointer<vtkLine>::New();
		vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
		vtkTree->GetLines()->GetNextCell(pts);
		(*it)->vtkSegment->GetPointIds()->SetId(0, pts->GetId(0));
		(*it)->vtkSegment->GetPointIds()->SetId(1, pts->GetId(1));
	}
}

void AbstractStructuredCCOTree::printVtkTree() {

	vtkTree->BuildCells();
	int numLines = vtkTree->GetNumberOfLines();
	int numPoints = vtkTree->GetNumberOfPoints();

	cout << "Tree with " << numLines << " segments and " << numPoints << " points." << endl;

	vtkSmartPointer<vtkCellArray> lines = vtkTree->GetLines();
	lines->InitTraversal();
	vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
	int i = 0;
	while (lines->GetNextCell(pts)) {
		cout << i++ << ": ";
		for (int j = 0; j < pts->GetNumberOfIds(); ++j) {
			cout << pts->GetId(j) << " ";
		}
		cout << endl;
	}

	vtkSmartPointer<vtkPoints> points = vtkTree->GetPoints();
	for (int i = 0; i < numPoints; ++i) {
		cout << i << ": (" << points->GetPoint(i)[0] << ", " << points->GetPoint(i)[1] << ", " << points->GetPoint(i)[2] << ")" << endl;
	}

}

GeneratorData* AbstractStructuredCCOTree::getInstanceData(){
	return instanceData;
}

void AbstractStructuredCCOTree::setInstanceData(GeneratorData* instanceData) {
	this->instanceData = instanceData;
}

vector<vector<double> > AbstractStructuredCCOTree::getVertices(){
	vector<vector<double> > points;
	vtkSmartPointer<vtkPoints> vtkPointsData = vtkTree->GetPoints(); // Points
	for (unsigned int i = 0; i < vtkPointsData->GetNumberOfPoints(); ++i) {
		vector<double> point;
		double *coordinates = vtkPointsData->GetPoint(i+1);
		point.push_back(coordinates[0]);
		point.push_back(coordinates[1]);
		point.push_back(coordinates[2]);
		points.push_back(point);
	}

	return points;

}

vector<vector<int> > AbstractStructuredCCOTree::getConnectivity(){
	vector<vector<int> > lines;
	for (std::vector<vessel *>::iterator it = segments.begin(); it != segments.end(); ++it) {
			vessel *currentSegment = *it;
			vector<int> line;
			vector<double> radius;

			vtkSmartPointer<vtkIdList> vtkLineData = currentSegment->vtkSegment->GetPointIds(); //	Lines
			line.push_back(vtkLineData->GetId(0));
			line.push_back(vtkLineData->GetId(1));
			lines.push_back(line);
		}
	return lines;
}
