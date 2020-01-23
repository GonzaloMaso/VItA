/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CCOTree.cpp
 *
 *  Created on: May 24, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "FRRVaViOptCCOSTree.h"

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

#include "../../core/GeneratorData.h"
#include "../domain/AbstractDomain.h"

FRRVaViOptCCOSTree::FRRVaViOptCCOSTree(point xi, double rootRadius, double qi, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *epsLim, AbstractConstraintFunction<double,int> *nu, double minAngle, double refPressure,
		double resistanceVariationTolerance, GeneratorData *instanceData) :
		AbstractStructuredCCOTree(xi, qi, gam, epsLim, nu, minAngle, refPressure, instanceData) {
	this->rootRadius = rootRadius;
	this->variationTolerance = resistanceVariationTolerance;
}

FRRVaViOptCCOSTree::FRRVaViOptCCOSTree(string filenameCCO, string filenameVTK, GeneratorData *instanceData) :
		AbstractStructuredCCOTree(instanceData) {
	ifstream treeFile;

	treeFile.open(filenameCCO.c_str(), ios::in);
	string token;
	treeFile >> token;
	while (token.compare("*Tree") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}

	//	Read tree
	treeFile >> xPerf.p[0];
	treeFile >> xPerf.p[1];
	treeFile >> xPerf.p[2];
	treeFile >> qProx;
	treeFile >> minAngle;
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

	//	Read each vessel
	int numVessels;
	treeFile >> numVessels;
	for (int i = 0; i < numVessels; ++i) {
		vessel *v = new vessel[1];
		treeFile >> v->vtkSegmentId;
		treeFile >> v->xProx.p[0];
		treeFile >> v->xProx.p[1];
		treeFile >> v->xProx.p[2];
		treeFile >> v->xDist.p[0];
		treeFile >> v->xDist.p[1];
		treeFile >> v->xDist.p[2];
		treeFile >> v->nLevel;
		treeFile >> v->radius;
		treeFile >> v->beta;
		treeFile >> v->length;
		treeFile >> v->resistance;
		treeFile >> v->flux;
		treeFile >> v->pressure;
		treeFile >> v->ID;
		treeFile >> v->treeVolume;
		this->segments.push_back(v);
	}

	this->root = segments[0];

	//	Load connectivity among segments
	treeFile >> token;
	while (token.compare("*Connectivity") != 0 && !treeFile.eof()) {
		treeFile >> token;
	}
	getline(treeFile, token);

	for (long long i = 0; i < numVessels; ++i) {
		getline(treeFile, token);
		stringstream ss;
		ss << token;
		int vtkId;
		ss >> vtkId;
		ss >> vtkId;
		while (!ss.eof()) {
			if (vtkId == -1)
				segments[i]->anastomose.push_back(NULL);
			else
				segments[i]->anastomose.push_back(segments[vtkId]);
			ss >> vtkId;
		}
	}
	cout << "Tree successfully loaded" << endl;
	cout << "Loading VTK structure..." << endl;

//	Load VTK tree from .vtk file
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filenameVTK.c_str());
	reader->Update();
	this->vtkTree = reader->GetOutput();
	cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
	cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

	updateSegmentVtkLines();

	vtkTree->BuildCells();
	vtkTree->Modified();

//	printVtkTree();

	vtkTreeLocator->SetDataSet(vtkTree);
	vtkTreeLocator->BuildLocator();
	vtkTreeLocator->Update();

	treeFile.close();

}

FRRVaViOptCCOSTree::~FRRVaViOptCCOSTree() {
}

double FRRVaViOptCCOSTree::getRootRadius() {
	return rootRadius;
}

void FRRVaViOptCCOSTree::getClosestTreePoint(point xNew, point *xBif, double *dist) {

	vtkIdType closeCellId;
	int subId;
	vtkTreeLocator->FindClosestPoint(xNew.p, xBif->p, closeCellId, subId, *dist);

	*dist = sqrt(*dist);
}

void FRRVaViOptCCOSTree::addVessel(point xProx, point xDist, vessel* parent) {
	nTerms++;
	//	Root
	if (!parent) {
		root = new vessel[1];

		//	Nodal quantities
		root->xDist = xDist;
		root->xProx = this->xPerf;
		point dist = root->xDist - root->xProx;
		root->nLevel = 0;
		root->beta = rootRadius;
		root->radius = rootRadius;
		root->length = sqrt(dist ^ dist);
		root->viscosity = nu->getValue(root->nLevel);
		root->resistance = (8 * root->viscosity / M_PI) * root->length;
		root->flux = qProx;
		root->treeVolume = M_PI * root->length * rootRadius * rootRadius;
		root->anastomose.push_back(NULL);
		root->ID = nTerms;
		//	TODO compute pressure
		root->pressure = 0.0;

		//	Tree quantities
		psiFactor = pow(root->beta, 4) / root->flux;
		dp = root->resistance / psiFactor;

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
		vtkTree->SetLines(lines);
		segments.push_back(root);

		//	Update tree locator
		vtkTreeLocator->SetDataSet(vtkTree);
		vtkTreeLocator->BuildLocator();
	}
	//	Non-root case
	else {

		//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
		point dNew = xDist - xProx;
		point dCon = parent->xDist - xProx;
		point dBif = xProx - parent->xProx;

		vessel *iNew = new vessel[1];
		iNew->xProx = xProx;
		iNew->xDist = xDist;
		iNew->nLevel = parent->nLevel + 1;
		iNew->length = sqrt(dNew ^ dNew);
		iNew->viscosity = nu->getValue(iNew->nLevel);
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->anastomose.push_back(parent);
		iNew->ID = nTerms;
		iNew->pressure = 0.0;

		vessel *iCon = new vessel[1];
		iCon->xProx = xProx;
		iCon->xDist = parent->xDist;
		iCon->nLevel = parent->nLevel + 1;
		iCon->length = sqrt(dCon ^ dCon);
		iCon->viscosity = nu->getValue(iCon->nLevel);
		iCon->anastomose.push_back(parent);
		iCon->ID = parent->ID;
		iCon->pressure = 0.0;

		if (parent->anastomose.size() > 1) {
			iCon->anastomose.push_back(parent->anastomose[1]);
			iCon->anastomose.push_back(parent->anastomose[2]);

			iCon->anastomose[1]->anastomose[0] = iCon;
			iCon->anastomose[2]->anastomose[0] = iCon;

			parent->anastomose[1] = iCon;
			parent->anastomose[2] = iNew;
		} else {
			iCon->resistance = 8 * iCon->viscosity / M_PI * iCon->length;
			parent->anastomose.push_back(iCon);
			parent->anastomose.push_back(iNew);
		}

		parent->xDist = xProx;
		parent->length = sqrt(dBif ^ dBif);

		segments.push_back(iNew);
		segments.push_back(iCon);

		//	Update post-order nLevel, flux and determine initial resistance and beta values.
		updateTree(this->root, this);

		//	Update resistance and betas
		double maxVariation = INFINITY;
		while (maxVariation > variationTolerance) {
			updateTreeViscositiesBeta(this->root, 1.0, &maxVariation);
		}

		//	Update tree geometry
		vtkIdType idProx = vtkTree->GetPoints()->InsertNextPoint(xProx.p);
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, idProx); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iCon->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iCon->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		iCon->vtkSegment->GetPointIds()->SetId(1, parent->vtkSegment->GetPointId(1)); // the second 1 is the index of xDist

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iNew->vtkSegment);
		iCon->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(iCon->vtkSegment);

//		cout << "Parent VTK Cell ids : " << vtkTree->GetCell(parent->vtkSegmentId)->GetPointIds()->GetNumberOfIds() << endl;
//		cout << "Intented modified id " << parent->vtkSegment->GetPointId(1) << endl;
		vtkTree->ReplaceCellPoint(parent->vtkSegmentId, parent->vtkSegment->GetPointId(1), idProx);
		parent->vtkSegment->GetPointIds()->SetId(1, idProx);

		vtkTree->BuildCells();
		vtkTree->Modified();

//		cout << "Points = " << vtkTree->GetNumberOfPoints() << endl;
//		cout << "Vessels = " << vtkTree->GetNumberOfLines() << endl;

		//	Update tree locator
		vtkTreeLocator->Update();
	}

}

vector<vessel*> FRRVaViOptCCOSTree::getCloseSegments(point xNew, AbstractDomain *domain, int* nFound) {

	vtkSmartPointer<vtkIdList> idSegments = vtkSmartPointer<vtkIdList>::New();
	double *localBox = domain->getLocalNeighborhood(xNew, nTerms);

	vtkTreeLocator->FindCellsWithinBounds(localBox, idSegments);

	vector<vessel*> closerSegments;
	int nElements = (int) idSegments->GetNumberOfIds();
	for (int i = 0; i < nElements; ++i) {
		closerSegments.push_back(segments[idSegments->GetId(i)]);
	}
	*nFound = closerSegments.size();
	delete[] localBox;
	return closerSegments;

}

int FRRVaViOptCCOSTree::testVessel(point xNew, vessel* parent, AbstractDomain *domain, vector<vessel *> neighbors, double dLim, point* xBif, double* cost) {
	point xP = parent->xProx;
	point xD = parent->xDist;

	int bifPartition = instanceData->nBifurcationTest;

	double ds = 1 / (double) (bifPartition - 1);
	int trials = (bifPartition) * (bifPartition + 1) / 2;
	double *costs = new double[trials];
	point *coordinates = new point[trials];

	//	All parametric coordinates are swept
	for (int i = 0; i < bifPartition; ++i) {
		for (int j = 0; j < bifPartition - i; ++j) {
			double eps = i * ds;
			double nu = j * ds;
			int majorIndex = i * bifPartition - i * (i - 1) / 2;

			coordinates[majorIndex + j] = xP * (1 - eps - nu) + xD * eps + xNew * nu;
			//	Invalid bifurcation positions
			if ((i == 0 && j == 0) || (i == bifPartition - 1) || (j == bifPartition - 1)) {
				costs[majorIndex + j] = INFINITY;
			} else {

				if (areValidAngles(coordinates[majorIndex + j], xNew, parent)) {
					if (domain->isSegmentInside(xNew, coordinates[majorIndex + j]) && domain->isSegmentInside(parent->xProx, coordinates[majorIndex + j])
							&& domain->isSegmentInside(parent->xDist, coordinates[majorIndex + j])) {
						if (!isIntersectingVessels(xNew, coordinates[majorIndex + j], parent, neighbors) && !isIntersectingVessels(parent->xProx, coordinates[majorIndex + j], parent, neighbors)
								&& !isIntersectingVessels(parent->xDist, coordinates[majorIndex + j], parent, neighbors)) {
							costs[majorIndex + j] = evaluate(xNew, coordinates[majorIndex + j], parent);
//						cout << "Cost for xNew " << xNew << " and " << parent->vtkSegmentId << " with bifurcation at " << coordinates[majorIndex + j-1] << " is " << costs[majorIndex + j-1] << endl;
						} else {
							costs[majorIndex + j] = INFINITY;
//					cout << "Intersection detected." << endl;
						}
					} else {
						costs[majorIndex + j] = INFINITY;
//				cout << "Cost for bifurcation at " << coordinates[majorIndex + j-1] << " connection outside the domain." << endl;
					}
				} else {
					costs[majorIndex + j] = INFINITY;
//					cout << "Small angle detected." << endl;
				}
			}
//#pragma omp critical
//			cout << "Cost of bifurcation at coordinates " << coordinates[majorIndex + j] << " is " << costs[majorIndex + j] << endl;
		}
	}

	*cost = INFINITY;
	*xBif = coordinates[0];
	for (int i = 0; i < trials; ++i) {
		if (costs[i] < *cost) {
			*cost = costs[i];
			*xBif = coordinates[i];
		}
	}

	delete[] costs;
	delete[] coordinates;

	return *cost != INFINITY;
}

double FRRVaViOptCCOSTree::evaluate(point xNew, point xTest, vessel* parent) {

	FRRVaViOptCCOSTree *clonedTree = cloneUpTo(instanceData->nLevelTest, parent);

	double previousCost = clonedTree->root->treeVolume;
	clonedTree->nTerms++;

	//	Fast-forward until parent in the cloned tree
	vector<vessel *>::iterator it = clonedTree->segments.begin();
	for (; (*it)->vtkSegmentId != parent->vtkSegmentId; ++it)
		;
	vessel *clonedParent = (*it);

	//	Add segment iNew, iCon and iBif in the cloned tree updating nLevel and lengths
	point dNew = xNew - xTest;
	point dCon = clonedParent->xDist - xTest;
	point dBif = xTest - clonedParent->xProx;

	vessel *iNew = new vessel[1];
	iNew->nLevel = clonedParent->nLevel + 1;
	iNew->length = sqrt(dNew ^ dNew);
	iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
	iNew->anastomose.push_back(clonedParent);

	vessel *iCon = new vessel[1];
	iCon->nLevel = clonedParent->nLevel + 1;
	iCon->length = sqrt(dCon ^ dCon);
	iCon->anastomose.push_back(clonedParent);
	if (clonedParent->anastomose.size() > 1) {
		iCon->anastomose.push_back(clonedParent->anastomose[1]);
		iCon->anastomose.push_back(clonedParent->anastomose[2]);

		iCon->anastomose[1]->anastomose[0] = iCon;
		iCon->anastomose[2]->anastomose[0] = iCon;

		clonedParent->anastomose[1] = iCon;
		clonedParent->anastomose[2] = iNew;
	} else {
		iCon->resistance = 8 * nu->getValue(iCon->nLevel) / M_PI * iCon->length;
		clonedParent->anastomose.push_back(iCon);
		clonedParent->anastomose.push_back(iNew);
	}
	clonedTree->segments.push_back(iNew);
	clonedTree->segments.push_back(iCon);

	clonedParent->length = sqrt(dBif ^ dBif);

	//	Update post-order nLevel, flux, initial resistances and intial betas.
	updateTree(clonedTree->root, clonedTree);

//	for (vector<vessel *>::iterator it = clonedTree->segments.begin(); it != clonedTree->segments.end(); ++it)
//		cout << (*it)->beta << endl;

	double maxVariation = INFINITY;
	while (maxVariation > variationTolerance) {
		updateTreeViscositiesBeta(clonedTree->root, 1.0, &maxVariation);
	}

	//	Check the symmetry constraint only for the newest vessel.
	if (!isSymmetricallyValid(iCon->beta, iNew->beta, iCon->nLevel)) {
		delete clonedTree;
		return INFINITY;
	}

	//	Compute cost and checks the geometric constraint only at the terminals - if the last is violated, cost is INFINITY
	double diffCost = computeTreeCost(clonedTree->root, 1.0) - previousCost;

	delete clonedTree;

	return diffCost;

}

int FRRVaViOptCCOSTree::isSymmetricallyValid(double beta1, double beta2, int nLevel) {
	double epsRad;
	if (beta1 > beta2)
		epsRad = beta2 / beta1;
	else
		epsRad = beta1 / beta2;

	return epsRad >= epsLim->getValue(nLevel);
}

void FRRVaViOptCCOSTree::print() {
	cout << "Printing FixedRadiusRootCCOTree" << endl;
	cout << "Root at " << xPerf << " with a radius of " << rootRadius << " mm and a flux of " << qProx << " cm^3/s" << endl;
	cout << "Pressure=" << dp << " Pa and psi factor=" << psiFactor << " cm x s" << endl;
	cout << "End terminals=" << nTerms << endl;
	cout << endl << "Topology" << endl;

	for (vector<vessel *>::iterator it = segments.begin(); it != segments.end(); ++it) {
		cout << (*it)->ID << ": " << (*it)->xProx << " -- " << (*it)->xDist << endl;
	}

}

int FRRVaViOptCCOSTree::isIntersectingVessels(point p1, point p2, vessel* parent, vector<vessel *> neighbors) {

	for (vector<vessel *>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
		double uv[2];
		int isIntersecting = 0;
		if (*it != parent) {
			isIntersecting = vtkLine::Intersection3D((*it)->xProx.p, (*it)->xDist.p, p1.p, p2.p, uv[0], uv[1]);
			if (isIntersecting && (uv[0] > 0 && uv[0] < 1) && (uv[1] > 0 && uv[1] < 1)) {
				return true;
			}
		}
	}
	return false;
}

FRRVaViOptCCOSTree* FRRVaViOptCCOSTree::clone() {
	FRRVaViOptCCOSTree *copy = new FRRVaViOptCCOSTree(this->xPerf, this->rootRadius, this->qProx, this->getGam(), this->getEpsLim(), this->getNu(), this->minAngle, this->refPressure,
			this->variationTolerance, this->instanceData);
	copy->psiFactor = this->psiFactor;
	copy->dp = this->dp;
	copy->nTerms = this->nTerms;

	copy->root = this->cloneTree(this->root, &(copy->segments));

	return copy;
}

vessel* FRRVaViOptCCOSTree::cloneTree(vessel* root, vector<vessel *> *segments) {

	vessel *copy = new vessel[1];

	copy->anastomose.push_back(NULL);
	copy->vtkSegmentId = root->vtkSegmentId;
	copy->xProx = root->xProx;
	copy->xDist = root->xDist;
	copy->nLevel = root->nLevel;
	copy->radius = root->radius;
	copy->beta = root->beta;
	copy->length = root->length;
	copy->resistance = root->resistance;
	copy->flux = root->flux;
	copy->viscosity = root->viscosity;
	copy->treeVolume = root->treeVolume;

	segments->push_back(copy);

	if (root->anastomose.size() > 1) {
		copy->anastomose.push_back(cloneTree(root->anastomose[1], segments));
		copy->anastomose.push_back(cloneTree(root->anastomose[2], segments));
		(copy->anastomose[1])->anastomose[0] = copy;
		(copy->anastomose[2])->anastomose[0] = copy;
	}

	return copy;
}

void FRRVaViOptCCOSTree::updateTree(vessel* root, FRRVaViOptCCOSTree* tree) {
	if (root->anastomose.size() > 1) {
		//	Is branch
		//	Update from subtrees
		vessel *vRight = root->anastomose[1];
		vessel *vLeft = root->anastomose[2];

		vRight->nLevel = root->nLevel + 1;
		vLeft->nLevel = root->nLevel + 1;

		updateTree(vRight, tree);
		updateTree(vLeft, tree);

		root->flux = vRight->flux + vLeft->flux;

		double betaRatio = sqrt(sqrt((vRight->flux * vRight->resistance) / (vLeft->flux * vLeft->resistance)));
		vRight->beta = pow(1 + pow(betaRatio, -gam->getValue(vRight->nLevel)), -1 / gam->getValue(vRight->nLevel));
		vLeft->beta = pow(1 + pow(betaRatio, gam->getValue(vLeft->nLevel)), -1 / gam->getValue(vLeft->nLevel));

		root->resistance = 8 * nu->getValue(root->nLevel) / M_PI * root->length + 1 / (pow(vRight->beta, 4) / vRight->resistance + pow(vLeft->beta, 4) / vLeft->resistance);

	} else {
		//	Is leaf
		root->flux = tree->qProx / tree->nTerms;
	}
}

int FRRVaViOptCCOSTree::areValidAngles(point xBif, point xNew, vessel* parent) {

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

int FRRVaViOptCCOSTree::isOverlapped(point xBif, point xNew, vessel* parent, double dLim) {

	point midPoint = xBif + (xNew - xBif) * 0.5;

	vtkSmartPointer<vtkIdList> idSegments = vtkSmartPointer<vtkIdList>::New();
	double *localBox = new double[6];
	double dLimScaled = dLim * instanceData->midPointDlimFactor;
	localBox[0] = midPoint.p[0] - dLimScaled;
	localBox[1] = midPoint.p[0] + dLimScaled;
	localBox[2] = midPoint.p[1] - dLimScaled;
	localBox[3] = midPoint.p[1] + dLimScaled;
	localBox[4] = midPoint.p[2] - dLimScaled;
	localBox[5] = midPoint.p[2] + dLimScaled;

	vtkTreeLocator->FindCellsWithinBounds(localBox, idSegments);
	delete[] localBox;

	vector<vessel*> closerSegments;
	int nElements = (int) idSegments->GetNumberOfIds();
	for (int i = 0; i < nElements; ++i) {
		vessel *current = segments[idSegments->GetId(i)];
		double *nearPoint;
		double t;
		if (current->vtkSegmentId != parent->vtkSegmentId)
			if (vtkLine::DistanceToLine(midPoint.p, current->xProx.p, current->xDist.p, t, nearPoint) < dLimScaled)
				return 1;
	}

	return 0;

}

void FRRVaViOptCCOSTree::updateTreeViscositiesResistance(vessel* root, double parentRadius, double* maxResistanceVariation) {
	double previousResistance = root->resistance;
	root->radius = root->beta * parentRadius;
	if (root->anastomose.size() > 1) {
		double maxResistanceVariation1, maxResistanceVariation2;
		vessel *vRight = root->anastomose[1];
		vessel *vLeft = root->anastomose[2];
		updateTreeViscositiesResistance(vRight, root->radius, &maxResistanceVariation1);
		updateTreeViscositiesResistance(vLeft, root->radius, &maxResistanceVariation2);

		// 	Estimate new betas and resistances
		double betaRatio = sqrt(sqrt((vRight->flux * vRight->resistance) / (vLeft->flux * vLeft->resistance)));
		vRight->beta = pow(1 + pow(betaRatio, -gam->getValue(vRight->nLevel)), -1 / gam->getValue(vRight->nLevel));
		vLeft->beta = pow(1 + pow(betaRatio, gam->getValue(vLeft->nLevel)), -1 / gam->getValue(vLeft->nLevel));

		root->resistance = 8 * getNuFL(root->radius) / M_PI * root->length + 1 / (pow(vRight->beta, 4) / vRight->resistance + pow(vLeft->beta, 4) / vLeft->resistance);

		root->treeVolume = root->radius * root->radius * M_PI * root->length + vLeft->treeVolume + vRight->treeVolume;

		*maxResistanceVariation = max(abs(root->resistance - previousResistance), max(maxResistanceVariation1, maxResistanceVariation2));
	} else {
		root->resistance = 8 * getNuFL(root->radius) / M_PI * root->length;
		root->treeVolume = root->radius * root->radius * M_PI * root->length;
		*maxResistanceVariation = abs(root->resistance - previousResistance);
	}
}

void FRRVaViOptCCOSTree::updateTreeViscositiesBeta(vessel* root, double parentRadius, double* maxBetaVariation) {
	root->radius = root->beta * parentRadius;
	if (root->anastomose.size() > 1) {
		double maxBetaVariation1, maxBetaVariation2;
		vessel *vRight = root->anastomose[1];
		vessel *vLeft = root->anastomose[2];
		updateTreeViscositiesBeta(vRight, root->radius, &maxBetaVariation1);
		updateTreeViscositiesBeta(vLeft, root->radius, &maxBetaVariation2);

		double previousBeta1 = vRight->beta;
		double previousBeta2 = vLeft->beta;

		// 	Estimate new betas and resistances
		double betaRatio = sqrt(sqrt((vRight->flux * vRight->resistance) / (vLeft->flux * vLeft->resistance)));
		vRight->beta = pow(1 + pow(betaRatio, -gam->getValue(vRight->nLevel)), -1 / gam->getValue(vRight->nLevel));
		vLeft->beta = pow(1 + pow(betaRatio, gam->getValue(vLeft->nLevel)), -1 / gam->getValue(vLeft->nLevel));

		double rBetaSqrt = vRight->beta * vRight->beta;
		double lBetaSqrt = vLeft->beta * vLeft->beta;
		root->viscosity = getNuFL(root->radius);
		root->resistance = 8 * root->viscosity / M_PI * root->length + 1 / (rBetaSqrt * rBetaSqrt / vRight->resistance + lBetaSqrt * lBetaSqrt / vLeft->resistance);

		root->treeVolume = root->radius * root->radius * M_PI * root->length + vLeft->treeVolume + vRight->treeVolume;

		double partialMaxVariation = max(abs(vRight->beta - previousBeta1), abs(vLeft->beta - previousBeta2));
		*maxBetaVariation = max(partialMaxVariation, max(maxBetaVariation1, maxBetaVariation2));
	} else {
		root->viscosity = getNuFL(root->radius);
		root->resistance = 8 * root->viscosity / M_PI * root->length;
		root->treeVolume = root->radius * root->radius * M_PI * root->length;
		*maxBetaVariation = 0.;
	}
}

FRRVaViOptCCOSTree* FRRVaViOptCCOSTree::cloneUpTo(int levels, vessel* parent) {

	vessel *subtreeRoot = parent;

	//	Identify N levels root
	for (int i = 0; i < levels && subtreeRoot->anastomose[0]; ++i) {
		subtreeRoot = subtreeRoot->anastomose[0];
	}

	FRRVaViOptCCOSTree *copy = new FRRVaViOptCCOSTree(this->xPerf, this->rootRadius, this->qProx, this->getGam(), this->getEpsLim(), this->getNu(), this->minAngle, this->refPressure,
			this->variationTolerance, this->instanceData);
	copy->psiFactor = this->psiFactor;
	copy->dp = this->dp;
	copy->nTerms = this->nTerms;

	copy->root = this->cloneTree(subtreeRoot, &(copy->segments));
	copy->root->beta = subtreeRoot->radius;
	copy->root->radius = subtreeRoot->radius;

	return copy;

}

inline double FRRVaViOptCCOSTree::getNuFL(double radius) {
	//	viscosity is in cP units; diameter is in microns.
	double d = radius * 20000;
//	cout << "current diameter in microns " << d << endl;
	double nuMixture = 6 * exp(-0.085 * d) - 2.44 * exp(-0.06 * pow(d, 0.645)) + 3.2;
	double nuPlasma = 1.1245;
	double relDSqr = d / (d - 1.1);
	relDSqr *= relDSqr;

	double viscosity = (nuPlasma * (1 + (nuMixture - 1) * relDSqr) * relDSqr) / 100;

//	cout << "The current viscosity is " << viscosity << endl;
	return viscosity;
}

void FRRVaViOptCCOSTree::saveTree(ofstream *outFile) {
	this->AbstractStructuredCCOTree::saveTree(outFile);
	*outFile << rootRadius << " " << variationTolerance << " ";
}

string FRRVaViOptCCOSTree::getTreeName() {
	return "FiRaRoVaViOptCCOTree";
}
