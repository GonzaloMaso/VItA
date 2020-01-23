/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MemoryMonitor.h
 *
 *  Created on: Mar 9, 2018
 *      Author: Gonzalo Maso Talou
 */

#ifndef MEMORYMONITOR_H_
#define MEMORYMONITOR_H_

/**
 * Class that report system memory usage.
 */
class MemoryMonitor {
	/** Unit conversion factor. */
	double unitFactor;
public:

	/**	Unit used as output. */
	enum UNITS {BYTE, KILOBYTE, MEGABYTE, GIGABYTE};
	/**
	 * Constructor.
	 * @param unit Unit used as output.
	 */
	MemoryMonitor(UNITS unit);
	/**
	 * Returns the RAM memory consumed by the current process.
	 */
	long long getProcessMemoryConsumption();
	/**
	 * Returns the RAM memory consumed by the whole system.
	 */
	long long getSystemMemoryConsumption();

private:
	int parseLine(char* line);
	long long getValue();

};

#endif /* MEMORYMONITOR_H_ */
