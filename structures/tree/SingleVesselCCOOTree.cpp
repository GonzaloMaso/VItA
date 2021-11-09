/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SingleVesselCCOOTree.cpp
 *
 *  Created on: Mar 29, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "SingleVesselCCOOTree.h"

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkIdList.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLReader.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "AbstractCostEstimator.h"
#include "../CCOCommonStructures.h"
#include "../domain/AbstractDomain.h"
#include "../vascularElements/SingleVessel.h"
#include "../../constrains/AbstractConstraintFunction.h"
#include "../../core/GeneratorData.h"

SingleVesselCCOOTree::SingleVesselCCOOTree(point xi, double rootRadius, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
		AbstractConstraintFunction<double, int> *nu, double refPressure, double resistanceVariationTolerance, GeneratorData *instanceData) :
		AbstractObjectCCOTree(xi, qi, gam, epsLim, nu, refPressure, instanceData) {
	this->filenameCCO = "";
	this->rootRadius = rootRadius;
	this->variationTolerance = resistanceVariationTolerance;
	this->nCommonTerminals = 0;
}

SingleVesselCCOOTree::SingleVesselCCOOTree(string filenameCCO, GeneratorData *instanceData, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
		AbstractConstraintFunction<double, int> *nu) :
		AbstractObjectCCOTree(instanceData) {
	this->filenameCCO = filenameCCO;
	ifstream treeFile;

	treeFile.open(filenameCCO.c_str(), ios::in);
	string token;
	treeFile >> token;
	while (token.compare("*Tree") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}

	cout << "Reading tree data..." << endl;
	//	Read tree
	treeFile >> xPerf.p[0];
	treeFile >> xPerf.p[1];
	treeFile >> xPerf.p[2];
	treeFile >> qProx;
	treeFile >> psiFactor;
	treeFile >> dp;
	treeFile >> nTerms;
	treeFile >> refPressure;
	treeFile >> pointCounter;
	treeFile >> rootRadius;
	treeFile >> variationTolerance;

	treeFile >> token;
	while (token.compare("*Vessels") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}

	cout << "Reading vessel data..." << endl;
	//	Read each vessel
	int numVessels;
	treeFile >> numVessels;
	//	Patch ------------------
	this->qReservedFactor = 0.0;
	//	------------------------
	double accReservedFlowFraction = 0.0;
	double currentReservedFlow = 0.0;
	for (int i = 0; i < numVessels; ++i) {
		int branchingMode = 0;
		int vesselType = 0;
		SingleVessel *v = new SingleVessel();
		v->qReservedFraction = 0.0;
		treeFile >> v->vtkSegmentId;
		cout << v->vtkSegmentId << endl;
		treeFile >> v->xProx.p[0];
		treeFile >> v->xProx.p[1];
		treeFile >> v->xProx.p[2];
		treeFile >> v->xDist.p[0];
		treeFile >> v->xDist.p[1];
		treeFile >> v->xDist.p[2];
		treeFile >> token; 						//	Tappering
		treeFile >> token;						//	Proximal_Distal_Radius_switched
		treeFile >> token;						//	Regions
		treeFile >> currentReservedFlow;		//	Reserved fraction -> Ask Paulo if it is alright
		if(currentReservedFlow > 0.0){
			v->terminalType = AbstractVascularElement::TERMINAL_TYPE::RESERVED;
			v->qReservedFraction = currentReservedFlow;
			accReservedFlowFraction += currentReservedFlow;
		}
		treeFile >> branchingMode; 				//	Branching - 0:NO_BRANCHING, 1:RIGID_PARENT, 2:DEFORMABLE_PARENT, 3:DISTAL_BRANCHING, 4:ONLY_AT_PARENT_HOTSPOTS
		v->branchingMode = static_cast<AbstractVascularElement::BRANCHING_MODE>(branchingMode);
		treeFile >> v->radius;
		treeFile >> token;						//	Resistance 1
		treeFile >> token;						//	Resistance 2
		treeFile >> token;						//	Capacitance
		treeFile >> token;						//	Pressure
		treeFile >> vesselType;					//	Perforators - 0:DISTRIBUTION, 1:PERFORATOR, 2:TRANSPORT
		v->vesselFunction = static_cast<AbstractVascularElement::VESSEL_FUNCTION>(vesselType);
		treeFile >> token;						//	Heart
		treeFile >> token;						//	Valves_SResistors_codes
		treeFile >> v->stage;					//	Stage
		this->elements[v->vtkSegmentId] = v;
	}

	this->qReservedFactor = accReservedFlowFraction;

	//	Load connectivity among segments
	treeFile >> token;
	cout << token << endl;
	while (token.compare("*Connectivity") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}
	getline(treeFile, token);

	cout << "Reading connectivity data..." << endl;
	int rootIndex = 0;
	for (long long i = 0; i < numVessels; ++i) {
		getline(treeFile, token);
		stringstream ss;
		ss << token;
		int vtkId, parentId, childId;

		ss >> vtkId;
		cout << "Vessel ID = " << vtkId;

		//	Parent parsing
		ss >> parentId;
		cout << " - P = " << parentId;
		if (parentId == -1) {
			elements[vtkId]->parent = NULL;
			rootIndex = vtkId;
		} else {
			elements[vtkId]->parent = elements[parentId];
		}

		//	Children parsing
		cout << " - Children : ";
		while (ss >> childId) {			
			cout << childId << " " ;
			elements[vtkId]->addChild(elements[childId]);
		}
		cout << endl;
	}

	cout << "Assembling tree..." << endl;
	//	Tree values
	this->root = elements[rootIndex];
	this->nTerms = this->getNTerminals();
	this->nCommonTerminals = getNTerminals(AbstractVascularElement::TERMINAL_TYPE::COMMON);
	cout << "Terminals " << nTerms << " - Common terminals " << getNTerminals(AbstractVascularElement::TERMINAL_TYPE::COMMON) << " - Reserved terminals " << getNTerminals(AbstractVascularElement::TERMINAL_TYPE::RESERVED) << endl;
	cout << "qReserved fraction " << qReservedFactor << endl;

	this->nu = nu;
	this->gam = gam;
	this->epsLim = epsLim;

	cout << "Tree successfully loaded" << endl;
	cout << "Creating VTK structure..." << endl;

	createSegmentVtkLines(root);

	vtkTree->BuildCells();
	vtkTree->Modified();

	cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
	cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

	vtkTreeLocator->SetDataSet(vtkTree);
	vtkTreeLocator->BuildLocator();
	vtkTreeLocator->Update();

	treeFile.close();
}

SingleVesselCCOOTree::SingleVesselCCOOTree(string filenameCCO, GeneratorData* instanceData, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
		AbstractConstraintFunction<double, int> *nu, double refPressure, double viscosityTolerance) :
		AbstractObjectCCOTree(instanceData) {

	this->filenameCCO = filenameCCO;

	ifstream treeFile;

	treeFile.open(filenameCCO.c_str(), ios::in);

	string token;
	treeFile >> token;
	while (token.compare("*Vessels") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}

	//	Read each vessel
	int numVessels;
	treeFile >> numVessels;
	//	Patch ------------------
	this->qReservedFactor = 0.0;
	//	------------------------
	double accReservedFlowFraction = 0.0;
	double currentReservedFlow = 0.0;
	for (int i = 0; i < numVessels; ++i) {
		int branchingMode = 0;
		int vesselType = 0;
		SingleVessel *v = new SingleVessel();
		v->qReservedFraction = 0.0;
		treeFile >> v->vtkSegmentId;
		cout << v->vtkSegmentId << endl;
		treeFile >> v->xProx.p[0];
		treeFile >> v->xProx.p[1];
		treeFile >> v->xProx.p[2];
		treeFile >> v->xDist.p[0];
		treeFile >> v->xDist.p[1];
		treeFile >> v->xDist.p[2];
		treeFile >> token; 						//	Tappering
		treeFile >> token;						//	Proximal_Distal_Radius_switched
		treeFile >> token;						//	Regions
		treeFile >> currentReservedFlow;		//	Reserved fraction -> Ask Paulo if it is alright
		if(currentReservedFlow > 0.0){
			v->terminalType = AbstractVascularElement::TERMINAL_TYPE::RESERVED;
			v->qReservedFraction = currentReservedFlow;
			accReservedFlowFraction += currentReservedFlow;
		}
		treeFile >> branchingMode; 				//	Branching - 0:NO_BRANCHING, 1:RIGID_PARENT, 2:DEFORMABLE_PARENT, 3:DISTAL_BRANCHING, 4:ONLY_AT_PARENT_HOTSPOTS
		v->branchingMode = static_cast<AbstractVascularElement::BRANCHING_MODE>(branchingMode);
		treeFile >> v->radius;
		treeFile >> token;						//	Resistance 1
		treeFile >> token;						//	Resistance 2
		treeFile >> token;						//	Capacitance
		treeFile >> token;						//	Pressure
		treeFile >> vesselType;					//	Perforators - 0:DISTRIBUTION, 1:PERFORATOR, 2:TRANSPORT
		v->vesselFunction = static_cast<AbstractVascularElement::VESSEL_FUNCTION>(vesselType);
		treeFile >> token;						//	Heart
		treeFile >> token;						//	Valves_SResistors_codes

		v->stage = -1;
		this->elements[v->vtkSegmentId] = v;
	}

	this->qReservedFactor = accReservedFlowFraction;

	//	Load connectivity among segments
	treeFile >> token;
	cout << token << endl;
	while (token.compare("*Connectivity") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}
	getline(treeFile, token);

	long long rootId = 0;
	for (long long i = 0; i < numVessels; ++i) {
		getline(treeFile, token);
		stringstream ss;
		ss << token;
		int vtkId, parentId, childId;

		ss >> vtkId;
		cout << "Vessel ID = " << vtkId;

		//	Parent parsing
		ss >> parentId;
		cout << " - P = " << parentId;
		if (parentId == -1) {
			elements[vtkId]->parent = NULL;
			rootId = vtkId;
		} else {
			elements[vtkId]->parent = elements[parentId];	//	BUG! vtkId is not the indexation of the vector
		}

		//	Children parsing
		cout << " - Children : ";
		while (ss >> childId) {			
			cout << childId << " " ;
			elements[vtkId]->addChild(elements[childId]);
//			ss >> childId;
		}
		cout << endl;
	}

	//	Tree values
	this->root = elements[rootId];
	this->xPerf = ((SingleVessel *) root)->xProx;
	this->rootRadius = ((SingleVessel *) root)->radius;
	this->qProx = qi;
	this->psiFactor = 0;
	this->dp = 0.0;
	this->nTerms = this->getNTerminals();
	this->nCommonTerminals = getNTerminals(AbstractVascularElement::TERMINAL_TYPE::COMMON);
	cout << "Terminals " << nTerms << " - Common terminals " << getNTerminals(AbstractVascularElement::TERMINAL_TYPE::COMMON) << " - Reserved terminals " << getNTerminals(AbstractVascularElement::TERMINAL_TYPE::RESERVED) << endl;
	cout << "qReserved fraction " << qReservedFactor << endl;
	this->refPressure = refPressure;
	this->pointCounter = 0l;
	this->variationTolerance = viscosityTolerance;

	this->nu = nu;
	this->gam = gam;
	this->epsLim = epsLim;

	cout << "Tree successfully loaded" << endl;
	cout << "Creating VTK structure..." << endl;

	createSegmentVtkLines(root);

	vtkTree->BuildCells();
	vtkTree->Modified();

	cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
	cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

	vtkTreeLocator->SetDataSet(vtkTree);
	vtkTreeLocator->BuildLocator();
	vtkTreeLocator->Update();

	treeFile.close();
}

SingleVesselCCOOTree::SingleVesselCCOOTree(SingleVesselCCOOTree *baseTree) : 
	AbstractObjectCCOTree(baseTree->instanceData) {
	// AbstractTreeCCOOTree attributes
	this->xPerf = baseTree->xPerf;
	this->qProx = baseTree->qProx;
	this->qReservedFactor = this->qReservedFactor;
	this->gam = baseTree->gam;
	this->epsLim = baseTree->epsLim;
	this->nu = baseTree->nu;
	this->refPressure = baseTree->refPressure;

	// SingleVesselCCOOTree attributes
	this->filenameCCO = this->filenameCCO;
	this->rootRadius = baseTree->rootRadius;
	this->variationTolerance = baseTree->variationTolerance;
}

SingleVesselCCOOTree::~SingleVesselCCOOTree() {
}

double SingleVesselCCOOTree::getRootRadius() {
	return rootRadius;
}

void SingleVesselCCOOTree::getClosestTreePoint(point xNew, point *xBif, double *dist) {

	vtkIdType closeCellId;
	int subId;
	vtkTreeLocator->FindClosestPoint(xNew.p, xBif->p, closeCellId, subId, *dist);

	*dist = sqrt(*dist);
}

void SingleVesselCCOOTree::addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction) {

	nTerms++;
	nCommonTerminals++;

	//	Root
	if (!parent) {

		SingleVessel * newRoot = new SingleVessel();

		//	Nodal quantities
		newRoot->xDist = xDist;
		newRoot->xProx = this->xPerf;
		point dist = newRoot->xDist - newRoot->xProx;
		newRoot->nLevel = 0;
		newRoot->beta = rootRadius;
		newRoot->radius = rootRadius;
		newRoot->length = sqrt(dist ^ dist);
		newRoot->viscosity = nu->getValue(newRoot->nLevel);
		newRoot->resistance = (8 * newRoot->viscosity / M_PI) * newRoot->length;
		newRoot->flow = qProx;
		newRoot->treeVolume = M_PI * newRoot->length * rootRadius * rootRadius;
		newRoot->parent = NULL;
		newRoot->ID = nTerms;
		newRoot->stage = currentStage;
		newRoot->pressure = newRoot->resistance * newRoot->flow + refPressure;
		newRoot->vesselFunction = vesselFunction;

		//	Tree quantities
		psiFactor = pow(newRoot->beta, 4) / newRoot->flow;	//	Not used
		dp = newRoot->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(newRoot->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(newRoot->xDist.p);
		vtkTree->SetPoints(pts);

		newRoot->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newRoot->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		newRoot->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		newRoot->vtkSegmentId = lines->InsertNextCell(newRoot->vtkSegment);
		vtkTree->SetLines(lines);
		elements[newRoot->vtkSegmentId] = newRoot;

		root = newRoot;

		//	Update tree locator
		vtkTreeLocator->SetDataSet(vtkTree);
		vtkTreeLocator->BuildLocator();
	}
	//	Non-root case & distal branching
	else if(parent->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING){

		if(parent->getChildren().empty()){
			nTerms--;
			nCommonTerminals--;
		}

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		parent->addChild(iNew);

		//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
		updateTree(((SingleVessel *) root), this);

		//	Update resistance, pressure and betas
		double maxVariation = INFINITY;
		while (maxVariation > variationTolerance) {
			updateTreeViscositiesBeta(((SingleVessel *) root), &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		elements[iNew->vtkSegmentId] = iNew;

		vtkTree->BuildCells();
		vtkTree->Modified();

		//	Update tree locator
		vtkTreeLocator->Update();

	}
	//	Non-root case & not distal branching
	else{

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;
		point dCon = ((SingleVessel *) parent)->xDist - xProx;
		point dBif = xProx - ((SingleVessel *) parent)->xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		SingleVessel *iCon = new SingleVessel();
		iCon->xProx = xProx;
		iCon->xDist = ((SingleVessel *) parent)->xDist;
		iCon->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iCon->length = sqrt(dCon ^ dCon);
		iCon->viscosity = nu->getValue(iCon->nLevel);
		iCon->parent = parent;
		iCon->ID = ((SingleVessel *) parent)->ID;
		iCon->branchingMode = parent->branchingMode;
		iCon->stage = ((SingleVessel *) parent)->stage;
		iCon->vesselFunction = ((SingleVessel *) parent)->vesselFunction;

		vector<AbstractVascularElement *> prevChildrenParent = parent->getChildren();
		if (prevChildrenParent.empty()) {
			iCon->resistance = 8 * iCon->viscosity / M_PI * iCon->length;
		} else {
			for (vector<AbstractVascularElement *>::iterator it = prevChildrenParent.begin(); it != prevChildrenParent.end(); ++it) {
				iCon->addChild(*it);
				(*it)->parent = iCon;
			}
			parent->removeChildren();
		}
		parent->addChild(iNew);
		parent->addChild(iCon);

		((SingleVessel *) parent)->xDist = xProx;
		((SingleVessel *) parent)->length = sqrt(dBif ^ dBif);

		//	Update post-order nLevel and flow, and determine initial resistance and beta values.
		updateTree(((SingleVessel *) root), this);

		//	Update resistance, pressure and betas
		double maxVariation = INFINITY;
		while (maxVariation > variationTolerance) {
			updateTreeViscositiesBeta(((SingleVessel *) root), &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idProx = vtkTree->GetPoints()->InsertNextPoint(xProx.p);
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, idProx); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iCon->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iCon->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		iCon->vtkSegment->GetPointIds()->SetId(1, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second 1 is the index of xDist

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		iCon->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iCon->vtkSegment);

		elements[iNew->vtkSegmentId] = iNew;
		elements[iCon->vtkSegmentId] = iCon;

//		cout << "Parent VTK Cell ids : " << vtkTree->GetCell(parent->vtkSegmentId)->GetPointIds()->GetNumberOfIds() << endl;
//		cout << "Intented modified id " << parent->vtkSegment->GetPointId(1) << endl;
		vtkTree->ReplaceCellPoint(((SingleVessel *) parent)->vtkSegmentId, ((SingleVessel *) parent)->vtkSegment->GetPointId(1), idProx);
		((SingleVessel *) parent)->vtkSegment->GetPointIds()->SetId(1, idProx);

		vtkTree->BuildCells();
		vtkTree->Modified();

//		cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
//		cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

		//	Update tree locator
		vtkTreeLocator->Update();
	}

}

void SingleVesselCCOOTree::addVesselMergeFast(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction,
	unordered_map<string, SingleVessel *>* stringToPointer) {
	printf("SingleVesselCCOOTree::addVesselMergeFast\n");
	nTerms++;
	nCommonTerminals++;

	//	Root
	if (!parent) {
		printf("Root parent\n");
		SingleVessel * newRoot = new SingleVessel();

		//	Nodal quantities
		newRoot->xDist = xDist;
		newRoot->xProx = this->xPerf;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(newRoot->coordToString(), newRoot));
		if(!didInsert.second) {
			printf("Did not insert new root vessel!\n");
		}
		point dist = newRoot->xDist - newRoot->xProx;
		newRoot->nLevel = 0;
		newRoot->beta = rootRadius;
		newRoot->radius = rootRadius;
		newRoot->length = sqrt(dist ^ dist);
		newRoot->viscosity = nu->getValue(newRoot->nLevel);
		newRoot->resistance = (8 * newRoot->viscosity / M_PI) * newRoot->length;
		newRoot->flow = qProx;
		newRoot->treeVolume = M_PI * newRoot->length * rootRadius * rootRadius;
		newRoot->parent = NULL;
		newRoot->ID = nTerms;
		newRoot->stage = currentStage;
		newRoot->pressure = newRoot->resistance * newRoot->flow + refPressure;
		newRoot->vesselFunction = vesselFunction;

		//	Tree quantities
		psiFactor = pow(newRoot->beta, 4) / newRoot->flow;	//	Not used
		dp = newRoot->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(newRoot->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(newRoot->xDist.p);
		vtkTree->SetPoints(pts);

		newRoot->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newRoot->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		newRoot->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		newRoot->vtkSegmentId = lines->InsertNextCell(newRoot->vtkSegment);
		vtkTree->SetLines(lines);
		elements[newRoot->vtkSegmentId] = newRoot;

		root = newRoot;

		//	Update tree locator
		vtkTreeLocator->SetDataSet(vtkTree);
		vtkTreeLocator->BuildLocator();
	}
	//	Non-root case & distal branching
	else if(parent->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING){
		printf("Distal branching\n");
		if(parent->getChildren().empty()){
			nTerms--;
			nCommonTerminals--;
		}

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iNew->coordToString(), iNew));
		if(!didInsert.second) {
			printf("Did not insert new distal vessel!\n");
		}
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		parent->addChild(iNew);

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		elements[iNew->vtkSegmentId] = iNew;

		vtkTree->BuildCells();
		vtkTree->Modified();

		//	Update tree locator
		vtkTreeLocator->Update();

	}
	//	Non-root case & not distal branching
	else{
		printf("Versatile branching\n");
		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;
		point dCon = ((SingleVessel *) parent)->xDist - xProx;
		point dBif = xProx - ((SingleVessel *) parent)->xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iNew->coordToString(), iNew));
		if(!didInsert.second) {
			printf("Did not insert new vessel!\n");
		}
		else {
			printf("Added new vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				iNew->xProx.p[0], iNew->xProx.p[1], iNew->xProx.p[2],
				iNew->xDist.p[0], iNew->xDist.p[1], iNew->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		SingleVessel *iCon = new SingleVessel();
		iCon->xProx = xProx;
		iCon->xDist = ((SingleVessel *) parent)->xDist;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iCon->coordToString(), iCon));
		if(!didInsert.second) {
			printf("Did not insert new sibling vessel!\n");
		}
		else {
			printf("Added sibling vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				iCon->xProx.p[0], iCon->xProx.p[1], iCon->xProx.p[2],
				iCon->xDist.p[0], iCon->xDist.p[1], iCon->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		iCon->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iCon->length = sqrt(dCon ^ dCon);
		iCon->viscosity = nu->getValue(iCon->nLevel);
		iCon->parent = parent;
		iCon->ID = ((SingleVessel *) parent)->ID;
		iCon->branchingMode = parent->branchingMode;
		iCon->stage = ((SingleVessel *) parent)->stage;
		iCon->vesselFunction = ((SingleVessel *) parent)->vesselFunction;

		vector<AbstractVascularElement *> prevChildrenParent = parent->getChildren();
		if (prevChildrenParent.empty()) {
			iCon->resistance = 8 * iCon->viscosity / M_PI * iCon->length;
		} else {
			for (vector<AbstractVascularElement *>::iterator it = prevChildrenParent.begin(); it != prevChildrenParent.end(); ++it) {
				iCon->addChild(*it);
				(*it)->parent = iCon;
			}
			parent->removeChildren();
		}
		parent->addChild(iNew);
		parent->addChild(iCon);

		SingleVessel *parentSV = (SingleVessel *) parent;
		size_t didErase = stringToPointer->erase(parentSV->coordToString());
		if(!didErase) {
			printf("Failed to erase parent at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[1], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[1], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		else {
			printf("Erased parent vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[1], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[1], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		parentSV->xDist = xProx;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(parentSV->coordToString(), parentSV));
		if(!didInsert.second) {
			printf("Did not insert updated parent vessel!\n");
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		else {
			printf("Added new parent vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[1], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[1], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		((SingleVessel *) parent)->length = sqrt(dBif ^ dBif);

		//	Update tree geometry
		vtkIdType idProx = vtkTree->GetPoints()->InsertNextPoint(xProx.p);
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, idProx); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iCon->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iCon->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		iCon->vtkSegment->GetPointIds()->SetId(1, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second 1 is the index of xDist

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		iCon->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iCon->vtkSegment);

		elements[iNew->vtkSegmentId] = iNew;
		elements[iCon->vtkSegmentId] = iCon;

//		cout << "Parent VTK Cell ids : " << vtkTree->GetCell(parent->vtkSegmentId)->GetPointIds()->GetNumberOfIds() << endl;
//		cout << "Intented modified id " << parent->vtkSegment->GetPointId(1) << endl;
		vtkTree->ReplaceCellPoint(((SingleVessel *) parent)->vtkSegmentId, ((SingleVessel *) parent)->vtkSegment->GetPointId(1), idProx);
		((SingleVessel *) parent)->vtkSegment->GetPointIds()->SetId(1, idProx);

		vtkTree->BuildCells();
		vtkTree->Modified();

//		cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
//		cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

		//	Update tree locator
		vtkTreeLocator->Update();
	}

}

void SingleVesselCCOOTree::addVesselMerge(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction,
	unordered_map<string, SingleVessel *>* stringToPointer) {
	printf("SingleVesselCCOOTree::addVesselMerge\n");
	nTerms++;
	nCommonTerminals++;

	//	Root
	if (!parent) {
		printf("Root parent\n");
		SingleVessel * newRoot = new SingleVessel();

		//	Nodal quantities
		newRoot->xDist = xDist;
		newRoot->xProx = this->xPerf;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(newRoot->coordToString(), newRoot));
		if(!didInsert.second) {
			printf("Did not insert new root vessel!\n");
		}
		point dist = newRoot->xDist - newRoot->xProx;
		newRoot->nLevel = 0;
		newRoot->beta = rootRadius;
		newRoot->radius = rootRadius;
		newRoot->length = sqrt(dist ^ dist);
		newRoot->viscosity = nu->getValue(newRoot->nLevel);
		newRoot->resistance = (8 * newRoot->viscosity / M_PI) * newRoot->length;
		newRoot->flow = qProx;
		newRoot->treeVolume = M_PI * newRoot->length * rootRadius * rootRadius;
		newRoot->parent = NULL;
		newRoot->ID = nTerms;
		newRoot->stage = currentStage;
		newRoot->pressure = newRoot->resistance * newRoot->flow + refPressure;
		newRoot->vesselFunction = vesselFunction;

		//	Tree quantities
		psiFactor = pow(newRoot->beta, 4) / newRoot->flow;	//	Not used
		dp = newRoot->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(newRoot->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(newRoot->xDist.p);
		vtkTree->SetPoints(pts);

		newRoot->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newRoot->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		newRoot->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		newRoot->vtkSegmentId = lines->InsertNextCell(newRoot->vtkSegment);
		vtkTree->SetLines(lines);
		elements[newRoot->vtkSegmentId] = newRoot;

		root = newRoot;

		//	Update tree locator
		vtkTreeLocator->SetDataSet(vtkTree);
		vtkTreeLocator->BuildLocator();
	}
	//	Non-root case & distal branching
	else if(parent->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING){
		printf("Distal branching\n");
		if(parent->getChildren().empty()){
			nTerms--;
			nCommonTerminals--;
		}

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iNew->coordToString(), iNew));
		if(!didInsert.second) {
			printf("Did not insert new distal vessel!\n");
		}
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		parent->addChild(iNew);

		//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
		updateTree(((SingleVessel *) root), this);

		//	Update resistance, pressure and betas
		double maxVariation = INFINITY;
		while (maxVariation > variationTolerance) {
			updateTreeViscositiesBeta(((SingleVessel *) root), &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		elements[iNew->vtkSegmentId] = iNew;

		vtkTree->BuildCells();
		vtkTree->Modified();

		//	Update tree locator
		vtkTreeLocator->Update();

	}
	//	Non-root case & not distal branching
	else{
		printf("Versatile branching\n");
		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;
		point dCon = ((SingleVessel *) parent)->xDist - xProx;
		point dBif = xProx - ((SingleVessel *) parent)->xProx;

		SingleVessel *iNew = new SingleVessel();
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iNew->coordToString(), iNew));
		if(!didInsert.second) {
			printf("Did not insert new vessel!\n");
		}
		else {
			printf("Added new vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				iNew->xProx.p[0], iNew->xProx.p[1], iNew->xProx.p[2],
				iNew->xDist.p[0], iNew->xDist.p[1], iNew->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->parent = parent;
		iNew->ID = nTerms;
		iNew->stage = currentStage;
		iNew->vesselFunction = vesselFunction;

		SingleVessel *iCon = new SingleVessel();
		iCon->xProx = xProx;
		iCon->xDist = ((SingleVessel *) parent)->xDist;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(iCon->coordToString(), iCon));
		if(!didInsert.second) {
			printf("Did not insert new sibling vessel!\n");
		}
		else {
			printf("Added sibling vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				iCon->xProx.p[0], iCon->xProx.p[1], iCon->xProx.p[2],
				iCon->xDist.p[0], iCon->xDist.p[1], iCon->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		iCon->nLevel = ((SingleVessel *) parent)->nLevel + 1;
		iCon->length = sqrt(dCon ^ dCon);
		iCon->viscosity = nu->getValue(iCon->nLevel);
		iCon->parent = parent;
		iCon->ID = ((SingleVessel *) parent)->ID;
		iCon->branchingMode = parent->branchingMode;
		iCon->stage = ((SingleVessel *) parent)->stage;
		iCon->vesselFunction = ((SingleVessel *) parent)->vesselFunction;

		vector<AbstractVascularElement *> prevChildrenParent = parent->getChildren();
		if (prevChildrenParent.empty()) {
			iCon->resistance = 8 * iCon->viscosity / M_PI * iCon->length;
		} else {
			for (vector<AbstractVascularElement *>::iterator it = prevChildrenParent.begin(); it != prevChildrenParent.end(); ++it) {
				iCon->addChild(*it);
				(*it)->parent = iCon;
			}
			parent->removeChildren();
		}
		parent->addChild(iNew);
		parent->addChild(iCon);

		SingleVessel *parentSV = (SingleVessel *) parent;
		size_t didErase = stringToPointer->erase(parentSV->coordToString());
		if(!didErase) {
			printf("Failed to erase parent at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[2], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[2], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		else {
			printf("Erased parent vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[2], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[2], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		parentSV->xDist = xProx;
		didInsert = stringToPointer->insert(pair<string, SingleVessel *>(parentSV->coordToString(), parentSV));
		if(!didInsert.second) {
			printf("Did not insert updated parent vessel!\n");
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		else {
			printf("Added new parent vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
				parentSV->xProx.p[0], parentSV->xProx.p[2], parentSV->xProx.p[2],
				parentSV->xDist.p[0], parentSV->xDist.p[2], parentSV->xDist.p[2]);
			printf("stringToPointer.size = %lu\n", stringToPointer->size());
		}
		((SingleVessel *) parent)->length = sqrt(dBif ^ dBif);

		//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
		updateTree(((SingleVessel *) root), this);

		//	Update resistance, pressure and betas
		double maxVariation = INFINITY;
		while (maxVariation > variationTolerance) {
			updateTreeViscositiesBeta(((SingleVessel *) root), &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idProx = vtkTree->GetPoints()->InsertNextPoint(xProx.p);
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, idProx); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iCon->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iCon->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		iCon->vtkSegment->GetPointIds()->SetId(1, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second 1 is the index of xDist

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		iCon->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iCon->vtkSegment);

		elements[iNew->vtkSegmentId] = iNew;
		elements[iCon->vtkSegmentId] = iCon;

//		cout << "Parent VTK Cell ids : " << vtkTree->GetCell(parent->vtkSegmentId)->GetPointIds()->GetNumberOfIds() << endl;
//		cout << "Intented modified id " << parent->vtkSegment->GetPointId(1) << endl;
		vtkTree->ReplaceCellPoint(((SingleVessel *) parent)->vtkSegmentId, ((SingleVessel *) parent)->vtkSegment->GetPointId(1), idProx);
		((SingleVessel *) parent)->vtkSegment->GetPointIds()->SetId(1, idProx);

		vtkTree->BuildCells();
		vtkTree->Modified();

//		cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
//		cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

		//	Update tree locator
		vtkTreeLocator->Update();
	}

}

void SingleVesselCCOOTree::addValitatedVessel(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo) {
	
	(this->nTerms)++;
	(this->nCommonTerminals++);

	SingleVessel *parentInNewTree = copiedTo[(SingleVessel *) originalVessel->parent];

	//	Root
	if (!parentInNewTree) {

		//	Nodal quantities
		newVessel->xDist = originalVessel->xDist;
		newVessel->xProx = this->xPerf;
		point dist = newVessel->xDist - newVessel->xProx;
		newVessel->nLevel = 0;
		newVessel->beta = originalVessel->beta;
		newVessel->radius = originalVessel->radius;
		newVessel->length = originalVessel->length;
		newVessel->viscosity = originalVessel->viscosity;
		newVessel->resistance = originalVessel->resistance;
		newVessel->flow = originalVessel->flow;
		newVessel->treeVolume = originalVessel->treeVolume;
		newVessel->parent = nullptr;
		newVessel->ID = nTerms;
		newVessel->stage = originalVessel->stage;
		newVessel->pressure = originalVessel->pressure;
		newVessel->vesselFunction = originalVessel->vesselFunction;

		//	Tree quantities
		this->psiFactor = pow(newVessel->beta, 4) / newVessel->flow;	//	Not used
		this->dp = newVessel->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(newVessel->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(newVessel->xDist.p);
		this->vtkTree->SetPoints(pts);

		newVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newVessel->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		newVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		newVessel->vtkSegmentId = lines->InsertNextCell(newVessel->vtkSegment);
		this->vtkTree->SetLines(lines);
		this->elements[newVessel->vtkSegmentId] = newVessel;

		this->root = newVessel;

		//	Update tree locator
		this->vtkTreeLocator->SetDataSet(vtkTree);
		this->vtkTreeLocator->BuildLocator();
	}
	//	Non-root case 
	// Because the vessel is already validated, it will always be distal
	else {

		if(parentInNewTree->getChildren().empty()) {
			nTerms--;
			nCommonTerminals--;
		}

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = originalVessel->xDist - originalVessel->xProx;

		newVessel->xProx = originalVessel->xProx;
		newVessel->xDist = originalVessel->xDist;
		newVessel->parent = parentInNewTree;
		newVessel->nLevel = (parentInNewTree->nLevel) + 1;
		newVessel->length = originalVessel->length;
		newVessel->viscosity = originalVessel->viscosity;
		newVessel->resistance = originalVessel->resistance;
		newVessel->ID = this->nTerms;
		newVessel->stage = originalVessel->stage;
		newVessel->vesselFunction = originalVessel->vesselFunction;

		newVessel->parent->addChild(newVessel);

		//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
		updateTree(((SingleVessel *) this->root), this);

		//	Update resistance, pressure and betas
		double maxVariation = INFINITY;
		while (maxVariation > this->variationTolerance) {
			updateTreeViscositiesBeta(((SingleVessel *) this->root), &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(newVessel->xDist.p);

		newVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newVessel->vtkSegment->GetPointIds()->SetId(0, parentInNewTree->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		newVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		newVessel->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(newVessel->vtkSegment);
		this->elements[newVessel->vtkSegmentId] = newVessel;

		this->vtkTree->BuildCells();
		this->vtkTree->Modified();

		//	Update tree locator
		this->vtkTreeLocator->Update();

	}	
}

void SingleVesselCCOOTree::addValitatedVesselFast(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo) {
	
	(this->nTerms)++;
	(this->nCommonTerminals++);

	SingleVessel *parentInNewTree = copiedTo[(SingleVessel *) originalVessel->parent];

	//	Root
	if (!parentInNewTree) {

		//	Nodal quantities
		newVessel->xDist = originalVessel->xDist;
		newVessel->xProx = this->xPerf;
		point dist = newVessel->xDist - newVessel->xProx;
		newVessel->nLevel = 0;
		newVessel->beta = originalVessel->beta;
		newVessel->radius = originalVessel->radius;
		newVessel->length = originalVessel->length;
		newVessel->viscosity = originalVessel->viscosity;
		newVessel->resistance = originalVessel->resistance;
		newVessel->flow = originalVessel->flow;
		newVessel->treeVolume = originalVessel->treeVolume;
		newVessel->parent = nullptr;
		newVessel->ID = nTerms;
		newVessel->stage = originalVessel->stage;
		newVessel->pressure = originalVessel->pressure;
		newVessel->vesselFunction = originalVessel->vesselFunction;

		//	Tree quantities
		this->psiFactor = pow(newVessel->beta, 4) / newVessel->flow;	//	Not used
		this->dp = newVessel->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(newVessel->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(newVessel->xDist.p);
		this->vtkTree->SetPoints(pts);

		newVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newVessel->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		newVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		newVessel->vtkSegmentId = lines->InsertNextCell(newVessel->vtkSegment);
		this->vtkTree->SetLines(lines);
		this->elements[newVessel->vtkSegmentId] = newVessel;

		this->root = newVessel;

		//	Update tree locator
		this->vtkTreeLocator->SetDataSet(vtkTree);
		this->vtkTreeLocator->BuildLocator();
	}
	//	Non-root case 
	// Because the vessel is already validated, it will always be distal
	else {

		if(parentInNewTree->getChildren().empty()) {
			nTerms--;
			nCommonTerminals--;
		}

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = originalVessel->xDist - originalVessel->xProx;

		newVessel->xProx = originalVessel->xProx;
		newVessel->xDist = originalVessel->xDist;
		newVessel->parent = parentInNewTree;
		newVessel->nLevel = (parentInNewTree->nLevel) + 1;
		newVessel->length = originalVessel->length;
		newVessel->viscosity = originalVessel->viscosity;
		newVessel->resistance = originalVessel->resistance;
		newVessel->ID = this->nTerms;
		newVessel->stage = originalVessel->stage;
		newVessel->vesselFunction = originalVessel->vesselFunction;

		newVessel->parent->addChild(newVessel);

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(newVessel->xDist.p);

		newVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		newVessel->vtkSegment->GetPointIds()->SetId(0, parentInNewTree->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		newVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		newVessel->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(newVessel->vtkSegment);
		this->elements[newVessel->vtkSegmentId] = newVessel;

		this->vtkTree->BuildCells();
		this->vtkTree->Modified();

		//	Update tree locator
		this->vtkTreeLocator->Update();

	}	
}

//void SingleVesselCCOOTree::addVessel(point xDist, AbstractVascularElement *parent, AbstractVascularElement::BRANCHING_MODE mode,
//		AbstractVascularElement::VESSEL_FUNCTION vesselFunction) {
//	//	Root
//	if (!parent) {
//
//		SingleVessel * newRoot = new SingleVessel();
//
//		//	Nodal quantities
//		newRoot->xDist = xDist;
//		newRoot->xProx = this->xPerf;
//		point dist = newRoot->xDist - newRoot->xProx;
//		newRoot->nLevel = 0;
//		newRoot->beta = rootRadius;
//		newRoot->radius = rootRadius;
//		newRoot->length = sqrt(dist ^ dist);
//		newRoot->viscosity = nu->getValue(newRoot->nLevel);
//		newRoot->resistance = (8 * newRoot->viscosity / M_PI) * newRoot->length;
//		newRoot->flow = qProx;
//		newRoot->treeVolume = M_PI * newRoot->length * rootRadius * rootRadius;
//		newRoot->parent = NULL;
//		newRoot->ID = nTerms;
//		newRoot->pressure = newRoot->resistance * newRoot->flow + refPressure;
//		newRoot->stage = currentStage;
//		newRoot->branchingMode = mode;
//		newRoot->vesselFunction = vesselFunction;
//
//		//	Tree quantities
//		psiFactor = pow(newRoot->beta, 4) / newRoot->flow;
//		dp = newRoot->resistance / psiFactor;
//
//		//	Update tree geometry
//		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
//		vtkIdType idProx = pts->InsertNextPoint(newRoot->xProx.p);
//		vtkIdType idDist = pts->InsertNextPoint(newRoot->xDist.p);
//		vtkTree->SetPoints(pts);
//
//		newRoot->vtkSegment = vtkSmartPointer<vtkLine>::New();
//		newRoot->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
//		newRoot->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
//		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//		newRoot->vtkSegmentId = lines->InsertNextCell(newRoot->vtkSegment);
//		vtkTree->SetLines(lines);
//		elements[newRoot->vtkSegmentId] = newRoot;
//
//		root = newRoot;
//
//		//	Update tree locator
//		vtkTreeLocator->SetDataSet(vtkTree);
//		vtkTreeLocator->BuildLocator();
//	}
//	//	Non-root case
//	else {
//
//		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
//		point dNew = xDist - ((SingleVessel *) parent)->xDist;
//
//		SingleVessel *iNew = new SingleVessel();
//		iNew->xProx = ((SingleVessel *) parent)->xDist;
//		iNew->xDist = xDist;
//		iNew->nLevel = ((SingleVessel *) parent)->nLevel + 1;
//		iNew->length = sqrt(dNew ^ dNew);
//		iNew->viscosity = nu->getValue(iNew->nLevel);
//		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
//		iNew->parent = parent;
//		iNew->ID = nTerms;
//		iNew->branchingMode = mode;
//		iNew->stage = currentStage;
//
//		iNew->radius = ((SingleVessel *) parent)->radius;
//		iNew->beta = 1;
//		iNew->vesselFunction = vesselFunction;
//
//		parent->addChild(iNew);
//
//		nTerms = getNTerminals();
//		nCommonTerminals = getNTerminals(AbstractVascularElement::TERMINAL_TYPE::COMMON);
//
//		//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
//		updateTree(((SingleVessel *) root), this);
//
//		//	Update resistance, pressure and betas
//		double maxVariation = INFINITY;
//		while (maxVariation > variationTolerance) {
//			updateTreeViscositiesBeta(((SingleVessel *) root), &maxVariation);
//		}
//
//		//	Update tree geometry
//		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);
//
//		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
//		iNew->vtkSegment->GetPointIds()->SetId(0, ((SingleVessel *) parent)->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
//		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point
//
//		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
//
//		elements[iNew->vtkSegmentId] = iNew;
//
//		vtkTree->BuildCells();
//		vtkTree->Modified();
//
//		//	Update tree locator
//		vtkTreeLocator->Update();
//	}
//
//}

vector<AbstractVascularElement*> SingleVesselCCOOTree::getCloseSegments(point xNew, AbstractDomain *domain, int* nFound) {

	vtkSmartPointer<vtkIdList> idSegments = vtkSmartPointer<vtkIdList>::New();
	double *localBox = domain->getLocalNeighborhood(xNew, nCommonTerminals);

	vtkTreeLocator->FindCellsWithinBounds(localBox, idSegments);

	vector<AbstractVascularElement*> closerSegments;
	int nElements = (int) idSegments->GetNumberOfIds();
	for (int i = 0; i < nElements; ++i) {
		int elemIndex = idSegments->GetId(i);
		AbstractVascularElement *candidate = elements[elemIndex];
		if(domain->isValidElement(candidate))
			if(candidate->branchingMode != AbstractVascularElement::NO_BRANCHING)
				closerSegments.push_back(candidate);
	}
	*nFound = closerSegments.size();
	delete[] localBox;
	return closerSegments;
}

int SingleVesselCCOOTree::testVessel(point xNew, AbstractVascularElement *parent, AbstractDomain *domain, vector<AbstractVascularElement *> neighbors, double dLim, point* xBif, double* cost) {

	vector<point> bifPoints;
	parent->getBranchingPoints(&bifPoints, xNew);
	SingleVessel *pVessel = (SingleVessel *) parent;

	vector<double> costs(bifPoints.size(), INFINITY);
	for (unsigned int i = 0; i < bifPoints.size(); ++i) {
		point bif = bifPoints[i];
		//	TODO Implement the BIG if as a filter design pattern for testing vessels. IMPORTANT! Benchmark that implementation against the hardcoded version to evaluate the performance since
		//	its a highly covered piece of the code. Advantages: can dynamically modify the checks at different stages to enhance computation.
		// Branching is distal or angles are valid
		if (pVessel->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING || (areValidAngles(bif, xNew, pVessel, domain->getMinBifurcationAngle())
				&&	isValidOpeningAngle(bif, xNew, pVessel, domain->getMinPlaneAngle()))
			) {
			/* x_n, bif is inside the domain ANDAND
			((Vessel is perforator OR x_p,x_b is inside) AND
			x_b, x_p is inside)
			In other words
			v_new is inside the domain AND
			(parent vessel is distal OR
			((v_p is inside the domain OR parente vessel is perforator) AND
			v_s is inside the domain))
			*/
			if (domain->isSegmentInside(xNew, bif) && (pVessel->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING ||
					((pVessel->vesselFunction == AbstractVascularElement::VESSEL_FUNCTION::PERFORATOR ||  domain->isSegmentInside(pVessel->xProx, bif)) && domain->isSegmentInside(pVessel->xDist, bif)) ) ) {
				/* v_new, v_s and v_p do not intersect neighbouring vessel */
				if (!isIntersectingVessels(xNew, bif, pVessel, neighbors) &&
						!isIntersectingVessels(pVessel->xProx, bif, pVessel, neighbors) &&
						!isIntersectingVessels(pVessel->xDist, bif, pVessel, neighbors)) {
					// Is distal
					if(pVessel->branchingMode == AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING){
						costs[i] = evaluate(xNew, pVessel, dLim);
					}
					// Is rigid/deformable/no_branching
					else{
						costs[i] = evaluate(xNew, bif, pVessel, dLim);
//						cout << "Cost for xNew " << xNew << " and " << parent->vtkSegmentId << " with bifurcation at " << coordinates[majorIndex + j-1] << " is " << costs[majorIndex + j-1] << endl;
					}
				} else {
					costs[i] = INFINITY;
					// cout << "Intersection detected." << endl;
				}
			} else {
				costs[i] = INFINITY;
				// cout << "Cost for bifurcation outside the domain." << endl;
			}
		} else {
			costs[i] = INFINITY;
			// cout << "Small angle detected." << endl;
		}
//#pragma omp critical
//			cout << "Cost of bifurcation at coordinates " << coordinates[majorIndex + j] << " is " << costs[majorIndex + j] << endl;
	}

	*cost = INFINITY;
	*xBif = {INFINITY,INFINITY,INFINITY};
	for (unsigned int i = 0; i < bifPoints.size(); ++i) {
		if (costs[i] < *cost) {
			*cost = costs[i];
			*xBif = bifPoints[i];
		}
	}

	return *cost != INFINITY;
}

double SingleVesselCCOOTree::evaluate(point xNew, point xTest, SingleVessel *parent, double dLim) {

	SingleVesselCCOOTree *clonedTree = cloneUpTo(instanceData->nLevelTest, parent);
//	SingleVesselCCOOTree *clonedTree = this->clone();

	AbstractCostEstimator *localEstimator = instanceData->costEstimator->clone();
	localEstimator->previousState(clonedTree, parent, xNew, xTest, dLim);

	clonedTree->nTerms++;
	clonedTree->nCommonTerminals++;

	//	Fast-forward until parent in the cloned tree
	auto it = clonedTree->elements.begin();
	for (; ((SingleVessel *) (it->second))->vtkSegmentId != ((SingleVessel *) parent)->vtkSegmentId; ++it)
		;
	SingleVessel *clonedParent = (SingleVessel *) (it->second);

	//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
	point dNew = xNew - xTest;
	point dCon = clonedParent->xDist - xTest;
	point dBif = xTest - clonedParent->xProx;

	SingleVessel *iNew = new SingleVessel();
	iNew->nLevel = clonedParent->nLevel + 1;
	iNew->length = sqrt(dNew ^ dNew);
	iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
	iNew->parent = clonedParent;

	SingleVessel *iCon = new SingleVessel();
	iCon->nLevel = clonedParent->nLevel + 1;
	iCon->length = sqrt(dCon ^ dCon);
	iCon->parent = clonedParent;

	vector<AbstractVascularElement *> prevChildrenParent = clonedParent->getChildren();
	if (prevChildrenParent.empty()) {
		iCon->resistance = 8 * nu->getValue(iCon->nLevel) / M_PI * iCon->length;
	} else {
		for (vector<AbstractVascularElement *>::iterator it = prevChildrenParent.begin(); it != prevChildrenParent.end(); ++it) {
			iCon->addChild(*it);
			(*it)->parent = iCon;
		}
		clonedParent->removeChildren();
	}
	clonedParent->addChild(iNew);
	clonedParent->addChild(iCon);

	//	Not needed because the updates use the tree structure to visit and update (not the element structure)
//	clonedTree->elements.push_back(iNew);
//	clonedTree->elements.push_back(iCon);

	clonedParent->length = sqrt(dBif ^ dBif);

	//	Update post-order nLevel, flux, initial resistances and intial betas.
	updateTree((SingleVessel *) clonedTree->root, clonedTree);

	double maxVariation = INFINITY;
	while (maxVariation > variationTolerance) {
		updateTreeViscositiesBeta((SingleVessel *) clonedTree->root, &maxVariation);
	}

	//	Check the symmetry constraint only for the newest vessel.
	if (!isSymmetricallyValid(iCon->beta, iNew->beta, iCon->nLevel)) {
		delete localEstimator;
		delete clonedTree;
		delete iNew;
		delete iCon;
		return INFINITY;
	}

	//	Compute cost and checks the geometric constraint only at the terminals - if the last is violated, cost is INFINITY
	double diffCost = localEstimator->computeCost(clonedTree);

	delete localEstimator;
	delete clonedTree;

	// As iCon and iNew are not added to clonedTree->elements we have to manually delete it.
	delete iNew;
	delete iCon;

	return diffCost;

}

double SingleVesselCCOOTree::evaluate(point xNew, SingleVessel *parent, double dLim) {

	SingleVesselCCOOTree *clonedTree = cloneUpTo(instanceData->nLevelTest, parent);
//	SingleVesselCCOOTree *clonedTree = this->clone();

	AbstractCostEstimator *localEstimator = instanceData->costEstimator->clone();
	localEstimator->previousState(clonedTree, parent, xNew, parent->xDist, dLim);

	if(parent->getChildren().size()>0){
		clonedTree->nTerms++;
		clonedTree->nCommonTerminals++;
	}

	//	Fast-forward until parent in the cloned tree
	auto it = clonedTree->elements.begin();
	for (; ((SingleVessel *) (it->second))->vtkSegmentId != ((SingleVessel *) parent)->vtkSegmentId; ++it)
		;
	SingleVessel *clonedParent = (SingleVessel *) (it->second);

	//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
	point dNew = xNew - clonedParent->xDist;
	point dBif = clonedParent->xDist - clonedParent->xProx;

	SingleVessel *iNew = new SingleVessel();
	iNew->nLevel = clonedParent->nLevel + 1;
	iNew->length = sqrt(dNew ^ dNew);
	iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
	iNew->parent = clonedParent;

	vector<AbstractVascularElement *> prevChildrenParent = clonedParent->getChildren();
	clonedParent->addChild(iNew);

	//	Same as in the other case, element structure is not needed for the following updates in the cloned tree.
//	clonedTree->elements.push_back(iNew);

	clonedParent->length = sqrt(dBif ^ dBif);

	//	Update post-order nLevel, flux, initial resistances and intial betas.
	updateTree((SingleVessel *) clonedTree->root, clonedTree);

	double maxVariation = INFINITY;
	while (maxVariation > variationTolerance) {
		updateTreeViscositiesBeta((SingleVessel *) clonedTree->root, &maxVariation);
	}

	//	FIXME Define symmetry law for N-ary bifurcations (Most different betas?)
	//	Check the symmetry constraint only for the newest vessel.
	if (!isSymmetricallyValid( ((SingleVessel *)clonedParent->getChildren()[0])->beta, iNew->beta, iNew->nLevel)) {
		delete localEstimator;
		delete clonedTree;
		delete iNew;
		return INFINITY;
	}

	//	Compute cost and checks the geometric constraint only at the terminals - if the last is unsatisfied, cost is INFINITY
	double diffCost = localEstimator->computeCost(clonedTree);

	delete localEstimator;
	delete clonedTree;

	// As iNew is not added to clonedTree->elements we have to manually delete it
	delete iNew;

	return diffCost;

}

int SingleVesselCCOOTree::isSymmetricallyValid(double beta1, double beta2, int nLevel) {
	double epsRad;
	if (beta1 > beta2)
		epsRad = beta2 / beta1;
	else
		epsRad = beta1 / beta2;

	return epsRad >= epsLim->getValue(nLevel);
}

void SingleVesselCCOOTree::print() {
	cout << "Printing FixedRadiusRootCCOTree" << endl;
	cout << "Root at " << xPerf << " with a radius of " << rootRadius << " mm and a flux of " << qProx << " cm^3/s" << endl;
	cout << "Pressure=" << dp << " Pa and psi factor=" << psiFactor << " cm x s" << endl;
	cout << "End terminals=" << nTerms << endl;
	cout << endl << "Topology" << endl;

	// TODO each vascular element must implement its own object printing
//	for (vector<vessel *>::iterator it = elements.begin(); it != elements.end(); ++it) {
//		cout << (*it)->print() << endl;
//	}

}

int SingleVesselCCOOTree::isIntersectingVessels(point p1, point p2, SingleVessel* parent, vector<AbstractVascularElement *> neighbors) {

	for (vector<AbstractVascularElement *>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
		double uv[2];
		SingleVessel* currentNeighbor = (SingleVessel*) (*it);
		int isIntersecting = 0;
		if (currentNeighbor != parent) {
			// TO DO check if this was changed in VTK 8.2 or 9.0
			isIntersecting = vtkLine::Intersection3D(currentNeighbor->xProx.p, currentNeighbor->xDist.p, p1.p, p2.p, uv[0], uv[1]);
			if (isIntersecting && (uv[0] > 0 && uv[0] < 1) && (uv[1] > 0 && uv[1] < 1)) {
				return true;
			}
		}
	}
	return false;
}

SingleVesselCCOOTree* SingleVesselCCOOTree::clone() {
	SingleVesselCCOOTree *copy = new SingleVesselCCOOTree(this->xPerf, this->rootRadius, this->qProx, this->getGam(), this->getEpsLim(), this->getNu(), this->refPressure,
			this->variationTolerance, this->instanceData);
	copy->nCommonTerminals = this->nCommonTerminals;
	copy->currentStage = this->currentStage;
	copy->qReservedFactor = this->qReservedFactor;
	copy->psiFactor = this->psiFactor;
	copy->dp = this->dp;
	copy->nTerms = this->nTerms;

	copy->root = this->cloneTree((SingleVessel *) root, &(copy->elements));

	return copy;
}

SingleVessel* SingleVesselCCOOTree::cloneTree(SingleVessel* root, unordered_map<long long, AbstractVascularElement *> *segments) {

	SingleVessel *copy = new SingleVessel();

	copy->parent = NULL;
	copy->vtkSegmentId = root->vtkSegmentId;
	copy->xProx = root->xProx;
	copy->xDist = root->xDist;
	copy->nLevel = root->nLevel;
	copy->radius = root->radius;
	copy->beta = root->beta;
	copy->length = root->length;
	copy->resistance = root->resistance;
	copy->flow = root->flow;
	copy->viscosity = root->viscosity;
	copy->treeVolume = root->treeVolume;

	(*segments)[copy->vtkSegmentId] = copy;

	vector<AbstractVascularElement *> rootChildren = root->getChildren();
	for (unsigned int i = 0; i < rootChildren.size(); ++i) {
		copy->children.push_back(cloneTree((SingleVessel*) rootChildren[i], segments));
		copy->children[i]->parent = copy;
	}

	return copy;
}

void SingleVesselCCOOTree::updateTree(SingleVessel* root, SingleVesselCCOOTree* tree) {
	if (root->getChildren().empty()) {
		root->flow = root->getTerminalFlow(tree->qProx, tree->qProx * tree->qReservedFactor, tree->nCommonTerminals); //tree->qProx / tree->nTerms;
//		cout << tree->qProx << " " << tree->qReservedFactor << " " << tree->nCommonTerminals << " " << root->flow << endl;
	} else {
		vector<AbstractVascularElement *> rootChildren = root->getChildren();
		double totalFlow = 0.0;
		double invTotalResistance = 0.0;
		for (vector<AbstractVascularElement *>::iterator it = rootChildren.begin(); it != rootChildren.end(); ++it) {
			SingleVessel *currentVessel = (SingleVessel *) (*it);
			currentVessel->nLevel = root->nLevel + 1;
			updateTree(currentVessel, tree);
			totalFlow += currentVessel->flow;
			invTotalResistance += 1 / currentVessel->resistance;
		}
		root->flow = totalFlow;

		double invResistanceContributions = 0.0;
		if (rootChildren.size() == 1) {
			SingleVessel *currentVessel = (SingleVessel *) rootChildren[0];
			currentVessel->beta = 1.0;
			invResistanceContributions = 1 / currentVessel->resistance;
		} else {
			for (vector<AbstractVascularElement *>::iterator it = rootChildren.begin(); it != rootChildren.end(); ++it) {
				SingleVessel *currentVessel = (SingleVessel *) (*it);
				double siblingsFlow = totalFlow - currentVessel->flow;
				double siblingsResistance = 1 / (invTotalResistance - 1 / currentVessel->resistance);
				double betaRatio = sqrt(sqrt((siblingsFlow * siblingsResistance) / (currentVessel->flow * currentVessel->resistance)));

				currentVessel->beta = pow(1 + pow(betaRatio, gam->getValue(currentVessel->nLevel)), -1.0 / gam->getValue(currentVessel->nLevel));

				double betaSqr = currentVessel->beta * currentVessel->beta;
				//	Check! This expression
				invResistanceContributions += betaSqr * betaSqr / currentVessel->resistance;
			}
		}
		root->localResistance = 8 * nu->getValue(root->nLevel) / M_PI * root->length;
		//	Check! Is not 1/ (1/localResistance + invResistanceContribution)?
		root->resistance = root->localResistance + 1 / invResistanceContributions;
	}
}

int SingleVesselCCOOTree::areValidAngles(point xBif, point xNew, SingleVessel* parent, double minAngle) {

	point iNew = xNew - xBif;
	point iCon = parent->xDist - xBif;
	point iBif = parent->xProx - xBif;

	double maxAngle = M_PI_2 - minAngle;

	double arg1 = (iNew ^ iCon) / sqrt((iNew ^ iNew) * (iCon ^ iCon));
	arg1 = min(1.0, max(-1.0, arg1));
	double arg2 = (iNew ^ iBif) / sqrt((iNew ^ iNew) * (iBif ^ iBif));
	arg2 = min(1.0, max(-1.0, arg2));

	double angle1 = abs(acos(arg1) - M_PI_2);
	double angle2 = abs(acos(arg2) - M_PI_2);

	if (angle1 > maxAngle || angle2 > maxAngle)
		return 0;

	return 1;
}

void SingleVesselCCOOTree::updateTreeViscositiesBeta(SingleVessel* root, double* maxBetaVariation) {
	if (root->parent) {
		root->radius = root->beta * ((SingleVessel*) root->parent)->radius;
	} else {
		root->radius = root->beta;
	}

	vector<AbstractVascularElement *> rootChildren = root->getChildren();
	if (rootChildren.empty()) {
		root->viscosity = getNuFL(root->radius);
		root->resistance = 8 * root->viscosity / M_PI * root->length;
		root->pressure = root->resistance * root->flow + refPressure;
		root->treeVolume = root->radius * root->radius * M_PI * root->length;
		*maxBetaVariation = 0.0;
	} else {

		double totalChildrenFlow = 0.0;
		double totalChildrenVolume = 0.0;
		double invTotalResistance = 0.0;
		*maxBetaVariation = 0.0;
		for (vector<AbstractVascularElement *>::iterator it = rootChildren.begin(); it != rootChildren.end(); ++it) {
			SingleVessel *currentVessel = (SingleVessel *) (*it);
			double betaVariation;
			updateTreeViscositiesBeta(currentVessel, &betaVariation);

			if (betaVariation > *maxBetaVariation)
				*maxBetaVariation = betaVariation;

			totalChildrenFlow += currentVessel->flow;
			totalChildrenVolume += currentVessel->treeVolume;
			invTotalResistance += 1 / currentVessel->resistance;
		}

		double invResistanceContributions = 0.0;
		if (rootChildren.size() == 1) {
			SingleVessel *currentVessel = (SingleVessel *) rootChildren[0];

			double previousBeta = currentVessel->beta;
			currentVessel->beta = 1.0;
			double betaVariation = abs(currentVessel->beta - previousBeta);
			if (betaVariation > *maxBetaVariation)
				*maxBetaVariation = betaVariation;

			invResistanceContributions = 1 / currentVessel->resistance;
		} else {
			for (vector<AbstractVascularElement *>::iterator it = rootChildren.begin(); it != rootChildren.end(); ++it) {
				SingleVessel *currentVessel = (SingleVessel *) (*it);
				double siblingsResistance = 1 / (invTotalResistance - 1 / currentVessel->resistance);
				double betaRatio = sqrt(sqrt(((totalChildrenFlow - currentVessel->flow) * siblingsResistance) / (currentVessel->flow * currentVessel->resistance)));
				double previousBeta = currentVessel->beta;
				currentVessel->beta = pow(1 + pow(betaRatio, gam->getValue(currentVessel->nLevel)), -1 / gam->getValue(currentVessel->nLevel));

				double betaVariation = abs(currentVessel->beta - previousBeta);
				if (betaVariation > *maxBetaVariation)
					*maxBetaVariation = betaVariation;

				double betaSqr = currentVessel->beta * currentVessel->beta;
				invResistanceContributions += betaSqr * betaSqr / currentVessel->resistance;
			}
		}

		root->viscosity = getNuFL(root->radius);
		root->localResistance = 8 * root->viscosity / M_PI * root->length;
		root->resistance = root->localResistance + 1 / invResistanceContributions;
		root->treeVolume = root->radius * root->radius * M_PI * root->length + totalChildrenVolume;
	}

	double radiusSqr = root->radius * root->radius;
	root->pressure = root->resistance * root->flow / (radiusSqr * radiusSqr) + refPressure;

}

SingleVesselCCOOTree* SingleVesselCCOOTree::cloneUpTo(int levels, SingleVessel* parent) {

	SingleVessel *subtreeRoot = parent;

	//	Identify N levels root
	for (int i = 0; i < levels && subtreeRoot->parent; ++i) {
		subtreeRoot = (SingleVessel*) subtreeRoot->parent;
	}

	SingleVesselCCOOTree *copy = new SingleVesselCCOOTree(this->xPerf, this->rootRadius, this->qProx, this->getGam(), this->getEpsLim(), this->getNu(), this->refPressure,
			this->variationTolerance, this->instanceData);
	copy->nCommonTerminals = this->nCommonTerminals;
	copy->currentStage = this->currentStage;
	copy->qReservedFactor = this->qReservedFactor;
	copy->psiFactor = this->psiFactor;
	copy->dp = this->dp;
	copy->nTerms = this->nTerms;

	copy->root = this->cloneTree(subtreeRoot, &(copy->elements));
	((SingleVessel *) copy->root)->beta = subtreeRoot->radius;
	((SingleVessel *) copy->root)->radius = subtreeRoot->radius;

	return copy;

}

/**
 * Function for radius in milimeters
 * @param radius Vessel radius in millimeters
 * @return Viscosity in centipoise
 */
inline double SingleVesselCCOOTree::getNuFL(double radius) {
	//	viscosity is in cP units; diameter is in microns.

	double d = radius * 2000;
	if(isInCm)
		d *= 10;
//	cout << "current diameter in microns " << d << endl;
	double nuMixture = 6 * exp(-0.085 * d) - 2.44 * exp(-0.06 * pow(d, 0.645)) + 3.2;
	double nuPlasma = 1.1245;
	double relDSqr = d / (d - 1.1);
	relDSqr *= relDSqr;

	double viscosity = (nuPlasma * (1 + (nuMixture - 1) * relDSqr) * relDSqr);// / 100;

//	cout << "The current viscosity is " << viscosity << endl;
	return viscosity;
}

void SingleVesselCCOOTree::saveTree(ofstream *outFile) {
	this->AbstractObjectCCOTree::saveTree(outFile);
	*outFile << rootRadius << " " << variationTolerance << " ";
}

string SingleVesselCCOOTree::getTreeName() {
	return "SingleVesselCCOOTree";
}

void SingleVesselCCOOTree::createSegmentVtkLines(AbstractVascularElement *vessel) {

	SingleVessel *currentVessel = (SingleVessel *) vessel;
	point xp = currentVessel->xProx;
	point xd = currentVessel->xDist;
	point dNew = xp - xd;

	//	BUG!! Parent may not have vtkSegment since it was not visited yet.
	SingleVessel *currentParent = (SingleVessel *) currentVessel->parent;
	if (currentParent) {
		currentVessel->nLevel = currentParent->nLevel + 1;
		currentVessel->beta = currentVessel->radius / currentParent->radius;
		currentVessel->length = sqrt(dNew ^ dNew);
		currentVessel->viscosity = nu->getValue(currentVessel->nLevel);
		currentVessel->resistance = 8 * nu->getValue(currentVessel->nLevel) / M_PI * currentVessel->length;
		currentVessel->pressure = 0.0;

		//	Update tree geometry
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xd.p);

		currentVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		currentVessel->vtkSegment->GetPointIds()->SetId(0, currentParent->vtkSegment->GetPointId(1)); // the second index is the global index of the mesh point
		currentVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		currentVessel->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(currentVessel->vtkSegment);

		vtkTree->BuildCells();
		vtkTree->Modified();
	} else {
		currentVessel->nLevel = 0;
		currentVessel->beta = rootRadius;
		currentVessel->radius = rootRadius;
		currentVessel->length = sqrt(dNew ^ dNew);
		currentVessel->viscosity = nu->getValue(currentVessel->nLevel);
		currentVessel->resistance = (8 * currentVessel->viscosity / M_PI) * currentVessel->length;
		currentVessel->flow = qProx;
		currentVessel->treeVolume = M_PI * currentVessel->length * rootRadius * rootRadius;
		currentVessel->pressure = 0.0;

		//	Tree quantities
		psiFactor = pow(currentVessel->beta, 4) / currentVessel->flow;
		dp = currentVessel->resistance / psiFactor;

		//	Update tree geometry
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		vtkIdType idProx = pts->InsertNextPoint(currentVessel->xProx.p);
		vtkIdType idDist = pts->InsertNextPoint(currentVessel->xDist.p);
		vtkTree->SetPoints(pts);

		currentVessel->vtkSegment = vtkSmartPointer<vtkLine>::New();
		currentVessel->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		currentVessel->vtkSegment->GetPointIds()->SetId(1, idDist); // the second 1 is the index of xDist
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		currentVessel->vtkSegmentId = lines->InsertNextCell(currentVessel->vtkSegment);
		vtkTree->SetLines(lines);

		vtkTree->BuildCells();
		vtkTree->Modified();
	}

	vector<AbstractVascularElement *> children = currentVessel->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		createSegmentVtkLines(*it);
	}

}

//	NEVER TESTED
void SingleVesselCCOOTree::removeWitheredBranches(int stage) {
	for (auto it = elements.begin(); it != elements.end(); ++it) {
		SingleVessel *currentVessel = (SingleVessel *) (it->second);
		if (currentVessel->stage == stage) {
			if (isWithered(currentVessel)) {
				this->remove(currentVessel);
				elements.erase(it);
			}
		}
	}
}

void SingleVesselCCOOTree::remove(SingleVessel* vessel) {

	vector<AbstractVascularElement *> children = vessel->getChildren();
	printf("children.size() = %lu\n", children.size());
	for (vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		remove( (SingleVessel*) *it);
	}
	// Vessel is root
	// if (!vessel->parent) {
	// 	return;
	// }
	vector<AbstractVascularElement *> parentChildren = vessel->parent->getChildren();
	printf("parentChildren.size() = %lu\n", parentChildren.size());
	vector<AbstractVascularElement *>::iterator it = parentChildren.begin();
	while (it != parentChildren.end()) {
		if (*it == vessel) {
			// We need to update the iterator this way, else we might run into trouble
			it = parentChildren.erase(it);
		}
		else {
			++it;
		}
	}

	printf("vtkCellType = %d\n", vessel->vtkSegment->GetCellType());
	printf("GetNumberOfPoints = %lld\n", vessel->vtkSegment->GetNumberOfPoints());
	vtkTree->DeletePoint(vessel->vtkSegment->GetPointId(1));
	vtkTree->DeleteCell(vessel->vtkSegmentId);	
	
	delete vessel;
}

//	NEVER TESTED
bool SingleVesselCCOOTree::isWithered(SingleVessel* vessel) {
	int currentStage = vessel->stage;
	vector<AbstractVascularElement *> children = vessel->getChildren();
	if (!children.empty()) {
		for (vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
			SingleVessel *child = (SingleVessel *) *it;
			if(child->stage > currentStage || !isWithered(child))
				return false;
		}
	}
	return true;
}

int SingleVesselCCOOTree::isValidOpeningAngle(point xBif, point xNew, SingleVessel* parent, double minPlaneAngle){

	point iNew = xNew - xBif;
	point iCon = parent->xDist - xBif;
	point iBif = parent->xProx - xBif;

	point planeNormal;

	planeNormal.p[0] = iBif.p[1] * iCon.p[2] - iBif.p[2] * iCon.p[1];
	planeNormal.p[1] = iBif.p[2] * iCon.p[0] - iBif.p[0] * iCon.p[2];
	planeNormal.p[2] = iBif.p[0] * iCon.p[1] - iBif.p[1] * iCon.p[0];

	//	acos return result in the interval of [0,pi]
	//	We're using the abs of the inner product since we are not interested in the orientation of the normal vector.
	double openingAngle = M_PI/2 - acos( abs(iNew ^ planeNormal) / sqrt((iNew^iNew) * (planeNormal^planeNormal)) );

	return minPlaneAngle <= openingAngle;

}

double SingleVesselCCOOTree::getVariationTolerance()
{
	return this->variationTolerance;
}

string SingleVesselCCOOTree::getFilenameCCO() {
	return this->filenameCCO;
}
