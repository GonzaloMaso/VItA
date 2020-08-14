/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CCOTree.cpp
 *
 *  Created on: May 24, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "FRRCCOSTree.h"

#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <math.h>

FRRCCOSTree::FRRCCOSTree(point xi,
		double rootRadius, double qi, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *epsLim,
		AbstractConstraintFunction<double,int> *nu, double minAngle, double refPressure, GeneratorData *instanceData) : AbstractStructuredCCOTree(xi, qi, gam, epsLim, nu, minAngle, refPressure, instanceData) {
	this->rootRadius = rootRadius;
}

double FRRCCOSTree::getRootRadius() {
	return rootRadius;
}

void FRRCCOSTree::getClosestTreePoint(point xNew, point *xBif,
		double *dist) {

	vtkIdType closeCellId;
	int subId;
	vtkTreeLocator->FindClosestPoint(xNew.p, xBif->p, closeCellId, subId,
			*dist);

	*dist = sqrt(*dist);
}

void FRRCCOSTree::addVessel(point xProx, point xDist,
		vessel* parent) {
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
		root->length = sqrt(dist ^ dist);
		root->resistance = (8 * nu->getValue(root->nLevel) / M_PI) * root->length;
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
		vtkSmartPointer<vtkCellArray> lines =
				vtkSmartPointer<vtkCellArray>::New();
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
		iNew->resistance = 8 * nu->getValue(iNew->nLevel) / M_PI * iNew->length;
		iNew->anastomose.push_back(parent);
		iNew->ID = nTerms;
		//	TODO compute pressure
		iNew->pressure = 0.0;

		vessel *iCon = new vessel[1];
		iCon->xProx = xProx;
		iCon->xDist = parent->xDist;
		iCon->nLevel = parent->nLevel + 1;
		iCon->length = sqrt(dCon ^ dCon);
		iCon->anastomose.push_back(parent);
		iCon->ID = parent->ID;
		//	TODO compute pressure
		iCon->pressure = 0.0;
		if (parent->anastomose.size() > 1) {
			iCon->anastomose.push_back(parent->anastomose[1]);
			iCon->anastomose.push_back(parent->anastomose[2]);

			parent->anastomose[1] = iCon;
			parent->anastomose[2] = iNew;
		} else {
			iCon->resistance = 8 * nu->getValue(iCon->nLevel) / M_PI * iCon->length;
			parent->anastomose.push_back(iCon);
			parent->anastomose.push_back(iNew);
		}

		parent->xDist = xProx;
		parent->length = sqrt(dBif ^ dBif);

		segments.push_back(iNew);
		segments.push_back(iCon);

		//	Update post-order nLevel, flux, resistance and beta symmetry constraints.
		updateTree(this->root, this);

		//	Update tree geometry
		vtkIdType idProx = vtkTree->GetPoints()->InsertNextPoint(xProx.p);
		vtkIdType idDist = vtkTree->GetPoints()->InsertNextPoint(xDist.p);

		iNew->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iNew->vtkSegment->GetPointIds()->SetId(0, idProx); // the second index is the global index of the mesh point
		iNew->vtkSegment->GetPointIds()->SetId(1, idDist); // the second index is the global index of the mesh point

		iCon->vtkSegment = vtkSmartPointer<vtkLine>::New();
		iCon->vtkSegment->GetPointIds()->SetId(0, idProx); // the second 0 is the index of xProx
		iCon->vtkSegment->GetPointIds()->SetId(1,
				parent->vtkSegment->GetPointId(1)); // the second 1 is the index of xDist

		iNew->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(
				iNew->vtkSegment);
		iCon->vtkSegmentId = vtkTree->GetLines()->InsertNextCell(
				iCon->vtkSegment);

//		cout << "Parent VTK Cell ids : " << vtkTree->GetCell(parent->vtkSegmentId)->GetPointIds()->GetNumberOfIds() << endl;
//		cout << "Intented modified id " << parent->vtkSegment->GetPointId(1) << endl;
		vtkTree->ReplaceCellPoint(parent->vtkSegmentId,
				parent->vtkSegment->GetPointId(1), idProx);
		parent->vtkSegment->GetPointIds()->SetId(1, idProx);
		vtkTree->BuildCells();
		vtkTree->Modified();

		//	Update tree locator
		vtkTreeLocator->Update();
	}

}

vector<vessel*> FRRCCOSTree::getCloseSegments(point xNew, AbstractDomain *domain,
		int* nFound) {

	vtkSmartPointer<vtkIdList> idSegments = vtkSmartPointer<vtkIdList>::New();
	double *localBox = domain->getLocalNeighborhood(xNew, nTerms);

	vtkTreeLocator->FindCellsWithinBounds(localBox, idSegments);

	vector<vessel*> closerSegments;
	int nElements = (int) idSegments->GetNumberOfIds();
	for (int i = 0; i < nElements; ++i) {
		closerSegments.push_back(segments[idSegments->GetId(i)]);
	}
	*nFound = closerSegments.size();
	return closerSegments;

}

/**
 *	Tests N_BIF_TRIES * (N_BIF_TRIES+1)/2 - 3 possible positions for the bifurcation site in the triangle
 *	delineated by the triplet (xNew, parent->xProx, parent->xDist). It returns the bifurcation site xBif
 *	which corresponds to the minimum cost in terms of volumetric area of the tree. If all bifurcations are
 *	invalid, cost will be DBL_MAX.
 */
int FRRCCOSTree::testVessel(point xNew, vessel* parent,
		AbstractDomain *domain, vector<vessel *> neighbors, double dLim,
		point* xBif, double* cost) {
	point xP = parent->xProx;
	point xD = parent->xDist;

	int bifPartition = instanceData->nBifurcationTest;

	double ds = 1 / (double) (bifPartition - 1);
	int trials = (bifPartition) * (bifPartition + 1) / 2;
	double *costs = new double[trials];
	point *coordinates = new point[trials];

//	cout << "Neighbor " << parent->vtkSegmentId << endl;
	//	All parametric coordinates are swept
	for (int i = 0; i < bifPartition; ++i) {
		for (int j = 0; j < bifPartition - i; ++j) {
			double eps = i * ds;
			double nu = j * ds;
			int majorIndex = i * bifPartition - i * (i - 1) / 2;

			coordinates[majorIndex + j] = xP * (1 - eps - nu) + xD * eps
					+ xNew * nu;
			//	Invalid bifurcation positions
			if ((i == 0 && j == 0) || (i == bifPartition - 1) || (j == bifPartition - 1)) {
				costs[majorIndex + j] = INFINITY;
			} else {

				//	The three created segmented must be within the domain: TODO Optimize by checking three segments at once
//			if (domain->isSegmentInside(xNew, coordinates[majorIndex + j ])
//					&& domain->isSegmentInside(parent->xProx,
//							coordinates[majorIndex + j ])
//					&& domain->isSegmentInside(parent->xDist,
//							coordinates[majorIndex + j ])) {
				//	If the proposed bifurcation point produces interpenetration or narrow bifurcation, the value set to INFINITY.
				if (!isIntersectingVessels(xNew, coordinates[majorIndex + j], parent, neighbors)
						&& !isIntersectingVessels(parent->xProx, coordinates[majorIndex + j], parent, neighbors)
						&& !isIntersectingVessels(parent->xDist, coordinates[majorIndex + j], parent, neighbors)) {
					if (areValidAngles(coordinates[majorIndex + j], xNew, parent)) {
//						if(!isOverlapped(coordinates[majorIndex + j], xNew, parent, dLim)){
							costs[majorIndex + j] = evaluate(xNew, coordinates[majorIndex + j], parent);
//					cout << "Cost for xNew " << xNew << " and " << parent << " with bifurcation at " << coordinates[majorIndex + j-1] << " is " << costs[majorIndex + j-1] << endl<< endl;
//						} else {
//							costs[majorIndex + j] = INFINITY;
//	//					cout << "Overlapping detected." << endl;
//						}
					} else {
						costs[majorIndex + j] = INFINITY;
//					cout << "Small angle detected." << endl;
					}
				} else {
					costs[majorIndex + j] = INFINITY;
//					cout << "Intersection detected." << endl;
				}
//			} else {
//				costs[majorIndex + j ] = INFINITY;
//				cout << "Cost for bifurcation at " << coordinates[majorIndex + j-1] << " connection outside the domain." << endl;
//			}
			}
		}
	}

	*cost = INFINITY;
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

double FRRCCOSTree::evaluate(point xNew, point xTest,
		vessel* parent) {

	FRRCCOSTree *clonedTree = this->clone();

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

		clonedParent->anastomose[1] = iCon;
		clonedParent->anastomose[2] = iNew;
	} else {
		iCon->resistance = 8 * nu->getValue(iCon->nLevel) / M_PI * iCon->length;
		clonedParent->anastomose.push_back(iCon);
		clonedParent->anastomose.push_back(iNew);
	}

	clonedParent->length = sqrt(dBif ^ dBif);

	//	Update post-order nLevel, flux, resistance and beta symmetry constraints.
	updateTree(clonedTree->root, clonedTree);

	//	Check the symmetry constraint only for the newest vessel.
	if (!isSymmetricallyValid(iCon->beta, iNew->beta, iCon->nLevel)) {
		delete clonedTree;
		return INFINITY;
	}

	//	Compute cost and checks the geometric constraint only at the terminals - if the last is violated, cost is INFINITY
	double cost = computeTreeCost(clonedTree->root, 1.0);

	delete clonedTree;

	return cost;

}

int FRRCCOSTree::isSymmetricallyValid(double beta1,
		double beta2, int nLevel) {
	double epsRad;
	if (beta1 > beta2)
		epsRad = beta2 / beta1;
	else
		epsRad = beta1 / beta2;

	return epsRad >= epsLim->getValue(nLevel);
}

void FRRCCOSTree::print() {
	cout << "Printing FixedRadiusRootCCOTree" << endl;
	cout << "Root at " << xPerf << " with a radius of " << rootRadius
			<< " mm and a flux of " << qProx << " cm^3/s" << endl;
	cout << "Pressure=" << dp << " Pa and psi factor=" << psiFactor << " cm x s"
			<< endl;
	cout << "End terminals=" << nTerms << endl;
	cout << endl << "Topology" << endl;

	for (vector<vessel *>::iterator it = segments.begin(); it != segments.end();
			++it) {
		cout << *it;
	}

}

int FRRCCOSTree::isIntersectingVessels(point p1, point p2,
		vessel* parent, vector<vessel *> neighbors) {

	for (vector<vessel *>::iterator it = neighbors.begin();
			it != neighbors.end(); ++it) {
		double u, v;
		int isIntersecting = 0;
		if (*it != parent) {
			isIntersecting = vtkLine::Intersection3D((*it)->xProx.p,
					(*it)->xDist.p, p1.p, p2.p, u, v);
			if (isIntersecting && (u > 0 && u < 1)) {
//				cout << "Intersects " << p1 << " - " << p2 << " with neighbor " << (void *) *it << " at " << u << " " << v << endl;
				return true;
			}
		}
	}
	return false;
}

FRRCCOSTree * FRRCCOSTree::clone() {
	FRRCCOSTree *copy = new FRRCCOSTree(
			this->xPerf, this->rootRadius, this->qProx, this->getGam(),
			this->getEpsLim(), this->getNu(), this->minAngle, this->refPressure, this->instanceData);
	copy->psiFactor = this->psiFactor;
	copy->dp = this->dp;
	copy->nTerms = this->nTerms;

	copy->root = this->cloneTree(this->root, &(copy->segments));

	return copy;
}

vessel* FRRCCOSTree::cloneTree(vessel* root,
		vector<vessel *> *segments) {

	vessel *copy = new vessel[1];

	copy->anastomose.push_back(NULL);
	copy->vtkSegmentId = root->vtkSegmentId;
	copy->xProx = root->xProx;
	copy->xDist = root->xDist;
	copy->nLevel = root->nLevel;
	copy->beta = root->beta;
	copy->length = root->length;
	copy->resistance = root->resistance;
	copy->flux = root->flux;

	segments->push_back(copy);

	if (root->anastomose.size() > 1) {
		copy->anastomose.push_back(cloneTree(root->anastomose[1], segments));
		copy->anastomose.push_back(cloneTree(root->anastomose[2], segments));
		(copy->anastomose[1])->anastomose[0] = copy;
		(copy->anastomose[2])->anastomose[0] = copy;
	}

	return copy;
}

void FRRCCOSTree::updateTree(vessel* root,
		FRRCCOSTree* tree) {
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

		double betaRatio = sqrt(sqrt(
				(vRight->flux * vRight->resistance)
						/ (vLeft->flux * vLeft->resistance)));
		vRight->beta = pow(1 + pow(betaRatio, -gam->getValue(vRight->nLevel)),
				-1 / gam->getValue(vRight->nLevel));
		vLeft->beta = pow(1 + pow(betaRatio, gam->getValue(vLeft->nLevel)),
				-1 / gam->getValue(vLeft->nLevel));

		root->resistance = 8 * nu->getValue(root->nLevel) / M_PI * root->length
				+ 1
						/ (pow(vRight->beta, 4) / vRight->resistance
								+ pow(vLeft->beta, 4) / vLeft->resistance);

	}
	else {
		//	Is leaf
		root->flux = tree->qProx / tree->nTerms;
	}
}

FRRCCOSTree::~FRRCCOSTree() {
}


int FRRCCOSTree::areValidAngles(point xBif, point xNew,
		vessel* parent) {
	point iNew = xNew - xBif;
	point iCon = parent->xDist - xBif;
	point iBif = parent->xProx - xBif;

	double maxAngle = M_PI_2 - minAngle;

	double arg1 = (iNew ^ iCon) / sqrt( (iNew ^ iNew) * (iCon ^ iCon) );
	arg1 = min(1.0,max(-1.0, arg1));
	double arg2 = (iNew ^ iBif) / sqrt( (iNew ^ iNew) * (iBif ^ iBif) );
	arg2 = min(1.0,max(-1.0, arg2));

	double angle1 = abs( acos( arg1 ) - M_PI_2 );
	double angle2 = abs( acos( arg2 ) - M_PI_2 );
	if (angle1 > maxAngle || angle2 > maxAngle)
		return 0;

	return 1;
}

int FRRCCOSTree::isOverlapped(point xBif, point xNew, vessel* parent, double dLim)	{

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
		if(current->vtkSegmentId != parent->vtkSegmentId)
			if(vtkLine::DistanceToLine(midPoint.p,current->xProx.p,current->xDist.p, t,nearPoint) < dLimScaled)
				return 1;
	}

	return 0;

}

void FRRCCOSTree::saveTree(ofstream *outFile) {
	this->AbstractStructuredCCOTree::saveTree(outFile);
	*outFile << rootRadius << " ";
}

string FRRCCOSTree::getTreeName() {
	return "FixedRadiusRootCCOTreeLegacy";
}
