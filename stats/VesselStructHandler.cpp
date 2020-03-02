/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselHandler.cpp
 *
 *  Created on: Feb 14, 2018
 *      Author: gonzalo
 */

#include "VesselStructHandler.h"

VesselStructHandler::VesselStructHandler(){

}

double VesselStructHandler::getVesselAttribute(vessel *v, ATTRIBUTE att){
	switch (att) {
	case DIAMETER:
		return v->radius*2;
	case RADIUS:
		return v->radius;
	case FLOW:
		return v->flux;
	case PRESSURE:
		return v->pressure;
	case RESISTANCE:
		return v->resistance;
	case LENGTH:
		return v->length;
	case LEVEL:
		return v->nLevel;
	case BETA:
		return v->beta;
	case VOLUME:
		return v->treeVolume;
	default:
		cout << "Unrecognize type of attribute." << endl;
		return NAN;
	}
}
