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
	if(currentTerminals > terminalsPerStage[currentStage]){
		terminalAtPrevStage = currentTerminals - 1;
		++currentStage;
		if(!domainStage[currentStage]->getInstanceData()->resetsDLim)
			domainStage[currentStage]->getInstanceData()->dLimCorrectionFactor = instanceData->dLimCorrectionFactor;
		instanceData = domainStage[currentStage]->getInstanceData();
		notifyObservers();
	}
}

long long int StagedDomain::getTerminalsAfterGeneration(){
	return terminalsPerStage.back();
}

int StagedDomain::isSegmentInside(point xs, point xf){
	return domainStage[currentStage]->isSegmentInside(xs,xf);
}

double StagedDomain::getCharacteristicLength(){
	return domainStage[currentStage]->getCharacteristicLength();
}

double StagedDomain::getDLim(long long int nVessels, double factor)	{
	if(instanceData->resetsDLim)
		return domainStage[currentStage]->getDLim(nVessels - terminalAtPrevStage,factor);
	else
		return domainStage[currentStage]->getDLim(nVessels,factor);
}

double* StagedDomain::getLocalNeighborhood(point p, long long int nVessels)	{
	return domainStage[currentStage]->getLocalNeighborhood(p,nVessels - terminalAtPrevStage);
}

double StagedDomain::getSize(){
	return domainStage[currentStage]->getSize();
}

point StagedDomain::getRandomPoint(){
	++pointCounter;
	return domainStage[currentStage]->getRandomPoint();
}

deque<point>& StagedDomain::getRandomInnerPoints(){
	return domainStage[currentStage]->getRandomInnerPoints();
}

vtkSmartPointer<vtkPolyData>& StagedDomain::getVtkGeometry(){
	return domainStage[currentStage]->getVtkGeometry();
}

int StagedDomain::getCurrentStage() const
{
	return currentStage;
}