/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * bloodViscosity.h
 *
 *  Created on: Jun 1, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef BLOODVISCOSITY_H_
#define BLOODVISCOSITY_H_

/**
 * Returns the same viscosity for all levels.
 * @param nLevel Bifurcation level.
 * @return Viscosity value.
 */
double constantViscosity(int nLevel);

inline double constantViscosity(int nLevel){
	return 0.036;
}

#endif /* BLOODVISCOSITY_H_ */
