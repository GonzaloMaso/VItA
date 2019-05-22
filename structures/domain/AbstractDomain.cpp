/*
 * AbstractDomain.cpp
 *
 *  Created on: 3 de dez de 2017
 *      Author: gonzalo
 */

#include "AbstractDomain.h"

AbstractDomain::AbstractDomain(GeneratorData *instanceData) {
	this->instanceData = instanceData;
	pointCounter = 0;
	isConvexDomain = false;
	volume = 0.0;
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
