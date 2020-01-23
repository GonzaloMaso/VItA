/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselHandler.h
 *
 *  Created on: Feb 14, 2018
 *      Author: gonzalo
 */

#ifndef STATS_VESSELSTRUCTHANDLER_H_
#define STATS_VESSELSTRUCTHANDLER_H_

#include "../structures/CCOCommonStructures.h"
/**
 * Adapter class that allows data extraction from vessel structure by attribute in a dynamic
 * manner.
 */
class VesselStructHandler{
public:

	/** Fields of vessel structure. */
	enum ATTRIBUTE {DIAMETER, RADIUS, FLOW, PRESSURE, RESISTANCE, LENGTH, LEVEL, BETA, VOLUME};

	/**
	 * Dummy constructor.
	 */
	VesselStructHandler();

	/**
	 * Returns the specific @p attribute of the @p v vessel.
	 * @param v	Vessel of interest.
	 * @param attribute	Field of interest.
	 * @return Value of the field @p attribute in vessel @p v.
	 */
	double getVesselAttribute(vessel *v, ATTRIBUTE attribute);
};

#endif /* STATS_VESSELSTRUCTHANDLER_H_ */
