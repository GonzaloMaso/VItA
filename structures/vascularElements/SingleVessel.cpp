/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SingleVessel.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "SingleVessel.h"

int SingleVessel::bifurcationTests = 6;

SingleVessel::SingleVessel() :
		AbstractVascularElement() {
	vessels.push_back(this);
	branchingMode = BRANCHING_MODE::DEFORMABLE_PARENT;
}

SingleVessel::~SingleVessel() {
}

AbstractVascularElement* SingleVessel::getParent() {
	return parent;
}

SingleVessel* SingleVessel::getParentVesselTo(point xp) {
	if (parent) {
		vector<SingleVessel*> *parents = parent->getVesselsConnectedTo(xp);
		SingleVessel *parent = (*parents)[0];
		delete parents;
		return parent;
	} else
		return NULL;
}

vector<AbstractVascularElement*>& SingleVessel::getChildren() {
	return children;
}

vector<SingleVessel*> *SingleVessel::getChildrenVesselTo(point xd) {
	vector<SingleVessel*> *vesselChildren = new vector<SingleVessel*>();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		vector<SingleVessel*> *partialChildren = (*it)->getVesselsConnectedTo(xd);
		vesselChildren->insert(vesselChildren->end(), partialChildren->begin(), partialChildren->end());
		delete partialChildren;
	}
	return vesselChildren;
}

vector<SingleVessel*> *SingleVessel::getVesselsConnectedTo(point p) {
	vector<SingleVessel*> *vessels = new vector<SingleVessel*>();
	if (p == xDist || p == xProx)
		vessels->push_back(this);
	return vessels;
}

long long int SingleVessel::getTerminals() {
	if (children.empty())
		return 1;
	else
		return 0;
}

double SingleVessel::getVolume() {
	radius = beta;
	if (parent) {
		radius *= parent->getDistalRadius();
	}
	return M_PI * radius * radius * length;
}

void SingleVessel::updatePressure() {
	localResistance = 8 * viscosity / M_PI * length;
	double distalPressure = 0;
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		distalPressure += (*it)->getProximalPressure();
	}
	pressure = flow * localResistance + distalPressure;
}

double SingleVessel::getDistalRadius() {
	return radius;
}

double SingleVessel::getProximalPressure() {
	return pressure;
}

void SingleVessel::saveVesselData(ofstream* treeFile) {
	*treeFile << vtkSegmentId << " " << xProx.p[0] << " " << xProx.p[1] << " " << xProx.p[2] << " " << xDist.p[0] << " " << xDist.p[1] << " " << xDist.p[2] << " "
			<< 0.0 << " " << 0.0 << " " << 0.0 << " " << qReservedFraction << " " << branchingMode << " " << radius << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "
			<< vesselFunction << " " << 0.0 << " " << 0.0 << " " << stage;
}

void SingleVessel::saveVesselConnectivity(ofstream* treeFile) {
	*treeFile << vtkSegmentId << " ";

	if (parent)
		*treeFile << ((SingleVessel *) parent)->vtkSegmentId;
	else
		*treeFile << -1;

	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		*treeFile << " " << ((SingleVessel *) (*it))->vtkSegmentId;
	}

}

void SingleVessel::getBranchingPoints(vector<point>* branchingPoints, point xNew) {

	int bifPartition = SingleVessel::bifurcationTests;
	double ds = 1 / (double) (bifPartition - 1);

	switch (branchingMode) {
	case RIGID_PARENT:
		for (int i = 1; i < bifPartition - 1; ++i) {
			double eps = i * ds;
			branchingPoints->push_back(xProx * (1 - eps) + xDist * eps);
		}
		break;
	case DISTAL_BRANCHING:
		branchingPoints->push_back(xDist);
		break;
	case NO_BRANCHING:
		branchingPoints->clear();
		break;
	default:

		for (int i = 0; i < bifPartition; ++i) {
			for (int j = 0; j < bifPartition - i; ++j) {
				double eps = i * ds;
				double nu = j * ds;

				//	Invalid bifurcation positions
				if (!(i == 0 && j == 0) && !(i == bifPartition - 1) && !(j == bifPartition - 1)) {
					branchingPoints->push_back(xProx * (1 - eps - nu) + xDist * eps + xNew * nu);
				}
			}
		}
		break;
	}
}

double SingleVessel::getTerminalFlow(double qProx, double qReserved, int commonTerminals){

	switch (terminalType) {
		case COMMON:
			return ( (qProx - qReserved) / commonTerminals);
			break;
		case RESERVED:
			return (qProx * qReservedFraction);
			break;
		default:
			return ( (qProx- qReserved) / commonTerminals);
			break;
	}
}

long long int SingleVessel::getTerminals(TERMINAL_TYPE type){
	if (children.empty() && terminalType == type)
		return 1;
	else
		return 0;

}

string SingleVessel::coordToString() {
    double coordArray[6] = {this->xProx.p[0], this->xProx.p[1], this->xProx.p[2],
		this->xDist.p[0], this->xDist.p[1], this->xDist.p[2]};
    char *coordCString = (char *) malloc(6 * sizeof(double));
    memcpy(coordCString, &coordArray, 6 * sizeof(double));
    string coordString(coordCString, (6 * sizeof(double) / sizeof(char)));
    free(coordCString);
    return coordString;
}