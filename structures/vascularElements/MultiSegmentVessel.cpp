/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MultiSegmentVessel.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "MultiSegmentVessel.h"

MultiSegmentVessel::MultiSegmentVessel(vector<SingleVessel*> innerSegments){
	parent = NULL;
	branchingMode = BRANCHING_MODE::DEFORMABLE_PARENT;
	length = 0;
	localResistance = INFINITY;
	treeVolume = 0;
	for(std::vector<SingleVessel *>::iterator it = innerSegments.begin(); it != innerSegments.end(); ++it) {
		vessels.push_back(*it);
		length += (*it)->length;
		localResistance += 1/(*it)->localResistance;
		treeVolume += (*it)->treeVolume;
	}
	localResistance = 1 / localResistance;
}

MultiSegmentVessel::~MultiSegmentVessel(){
}

AbstractVascularElement* MultiSegmentVessel::getParent(){
	return parent;
}

SingleVessel* MultiSegmentVessel::getParentVesselTo(point xp){
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if((*it)->xDist == xp)
			return *it;
	}
	return NULL;
}

vector<AbstractVascularElement*>& MultiSegmentVessel::getChildren(){
	return children;
}

vector<SingleVessel*>* MultiSegmentVessel::getChildrenVesselTo(point xd){
	vector<SingleVessel*> *childrenVessels = new vector<SingleVessel *>();
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if((*it)->xProx == xd)
			childrenVessels->push_back(*it);
	}

	return childrenVessels;
}

vector<SingleVessel*>* MultiSegmentVessel::getVesselsConnectedTo(point p){
	vector<SingleVessel*> *childrenVessels = new vector<SingleVessel *>();
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if((*it)->xProx == p || (*it)->xDist == p)
			childrenVessels->push_back(*it);
	}

	return childrenVessels;
}

long long int MultiSegmentVessel::getTerminals(){
	if ((vessels.back())->getChildren().size() == 0)
		return 1;
	else
		return 0;
}

double MultiSegmentVessel::getVolume(){
	double cost = 0.0;
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		SingleVessel *currentVessel = *it;
		cost += currentVessel->radius * currentVessel->radius * M_PI * currentVessel->length;
	}

	return cost;
}

void MultiSegmentVessel::updatePressure(){
	for(std::vector<SingleVessel *>::iterator it = vessels.end(); it != vessels.begin(); --it) {
		(*it)->updatePressure();
	}
}

double MultiSegmentVessel::getDistalRadius(){
	return (vessels.back())->radius;
}

double MultiSegmentVessel::getProximalPressure(){
	return (vessels.front())->pressure;
}

void MultiSegmentVessel::saveVesselData(ofstream* treeFile){
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		(*it)->saveVesselData(treeFile);
	}
}

void MultiSegmentVessel::saveVesselConnectivity(ofstream* treeFile){
	for(std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		(*it)->saveVesselConnectivity(treeFile);
	}
}
