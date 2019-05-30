/*
 * AbstractDomain.cpp
 *
 *  Created on: 3 de dez de 2017
 *      Author: gonzalo
 */

#include "AbstractDomain.h"

#include <algorithm>

AbstractDomain::AbstractDomain(GeneratorData *instanceData) {
	this->instanceData = instanceData;
	pointCounter = 0;
	isConvexDomain = false;
	volume = 0.0;
	minAngle = M_PI * 3. / 18.;
}

AbstractDomain::AbstractDomain(GeneratorData* instanceData, vector<int> growingStages){
	this->instanceData = instanceData;
	pointCounter = 0;
	isConvexDomain = false;
	volume = 0.0;
	minAngle = M_PI * 3. / 18.;
	this->growingStages = growingStages;
}

AbstractDomain::~AbstractDomain() {
}

long long int AbstractDomain::getPointCounter() const {
	return pointCounter;
}

bool AbstractDomain::isIsConvexDomain() const {
	return isConvexDomain;
}

void AbstractDomain::setIsConvexDomain(bool isConvexDomain)	{
	this->isConvexDomain = isConvexDomain;
}

GeneratorData* AbstractDomain::getInstanceData(){
	return instanceData;
}

void AbstractDomain::setInstanceData(GeneratorData* instanceData) {
	this->instanceData = instanceData;
}

int AbstractDomain::isValidElement(AbstractVascularElement* element){
	return (growingStages.empty() || find(growingStages.begin(), growingStages.end(), element->stage) != growingStages.end());
}

const vector<int>& AbstractDomain::getGrowingStages() const
{
	return growingStages;
}

void AbstractDomain::setGrowingStages(const vector<int>& growingStages){
	this->growingStages = growingStages;
}

double AbstractDomain::getMinBifurcationAngle(){
	return minAngle;
}

void AbstractDomain::setMinBifurcationAngle(double minAngle){
	this->minAngle = minAngle;
}
