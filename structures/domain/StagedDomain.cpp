/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StagedDomain.cpp
 *
 *  Created on: Mar 14, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "StagedDomain.h"

StagedDomain::StagedDomain() : AbstractDomain(NULL){
	currentTerminals = 0;
	terminalAtPrevStage = 0;
	currentStage = 0;
	initialStage = 0;
}

StagedDomain::~StagedDomain() {
	this->domainStage.clear();
}

void StagedDomain::addStage(long long int terminals, AbstractDomain* domain){
	domainStage.push_back(domain);
	if(terminalsPerStage.size()>0){
		terminalsPerStage.push_back(terminalsPerStage.back()+terminals);
	}
	else{
		terminalsPerStage.push_back(terminals);
		instanceData = domain->getInstanceData();
	}
}

void StagedDomain::update(){
	++currentTerminals;
	if(currentTerminals > terminalsPerStage[currentStage-initialStage]){
		terminalAtPrevStage = currentTerminals - 1;
		++currentStage;
		if(!domainStage[currentStage-initialStage]->getInstanceData()->resetsDLim)
			domainStage[currentStage-initialStage]->getInstanceData()->dLimCorrectionFactor = instanceData->dLimCorrectionFactor;
		instanceData = domainStage[currentStage-initialStage]->getInstanceData();
		notifyObservers();
	}
}

long long int StagedDomain::getTerminalsAfterGeneration(){
	return terminalsPerStage.back();
}

int StagedDomain::isSegmentInside(point xs, point xf){
	return domainStage[currentStage-initialStage]->isSegmentInside(xs,xf);
}

double StagedDomain::getCharacteristicLength(){
	return domainStage[currentStage-initialStage]->getCharacteristicLength();
}

double StagedDomain::getDLim(long long int nVessels, double factor)	{
	if(instanceData->resetsDLim)
		return domainStage[currentStage-initialStage]->getDLim(nVessels - terminalAtPrevStage,factor);
	else
		return domainStage[currentStage-initialStage]->getDLim(nVessels,factor);
}

double* StagedDomain::getLocalNeighborhood(point p, long long int nVessels)	{
	return domainStage[currentStage-initialStage]->getLocalNeighborhood(p,nVessels - terminalAtPrevStage);
}

double StagedDomain::getSize(){
	return domainStage[currentStage-initialStage]->getSize();
}

point StagedDomain::getRandomPoint(){
	return domainStage[currentStage-initialStage]->getRandomPoint();
}

deque<point>& StagedDomain::getRandomInnerPoints(){
	return domainStage[currentStage-initialStage]->getRandomInnerPoints();
}

vtkSmartPointer<vtkPolyData>& StagedDomain::getVtkGeometry(){
	return domainStage[currentStage-initialStage]->getVtkGeometry();
}

int StagedDomain::getCurrentStage() const
{
	return currentStage;
}

int StagedDomain::isValidElement(AbstractVascularElement* element){
	return domainStage[currentStage-initialStage]->isValidElement(element);
}

double StagedDomain::getMinBifurcationAngle(){
	return domainStage[currentStage-initialStage]->getMinBifurcationAngle();
}

void StagedDomain::setInitialStage(int currentStage){
	try{
		if(this->currentStage!=0)
			throw 1;
	} catch(int idxExeception){
		cout << "Initial stage can only be done once and before execution" << endl;
	}
	this->initialStage = currentStage;
	this->currentStage = this->initialStage;
}

double StagedDomain::getMinPlaneAngle(){
	return domainStage[currentStage-initialStage]->getMinPlaneAngle();
}

long long int StagedDomain::getPointCounter() const{
	return domainStage[currentStage-initialStage]->getPointCounter();
}

vector<AbstractDomain *>* StagedDomain::getDomains()
{
	return &(this->domainStage);
}

vector<long long int>* StagedDomain::getNTerminals()
{
	return &(this->terminalsPerStage);
}

int StagedDomain::getDraw()
{
	return -1;
}

int StagedDomain::getSeed()
{
	return -1;
}

void StagedDomain::logDomainFiles(FILE *fp) {
	fprintf(fp, "StagedDomain\n");
}

vtkSmartPointer<vtkSelectEnclosedPoints> StagedDomain::getEnclosedPoints() {
	return this->domainStage[(this->currentStage)-(this->initialStage)]->getEnclosedPoints();
}