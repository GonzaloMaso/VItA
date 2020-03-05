/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MemoryMonitor.cpp
 *
 *  Created on: Mar 9, 2018
 *      Author: Gonzalo D. Maso Talou
 *
 *  Code from parseLine and getValue developed by Calvin1602 ( https://stackoverflow.com/users/124038/calvin1602 )
 */

#include "MemoryMonitor.h"

#include <sys/types.h>
#include <sys/sysinfo.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct sysinfo memInfo;

MemoryMonitor::MemoryMonitor(UNITS unit){
	switch (unit) {
		case KILOBYTE:
			unitFactor = 1./1024.;
			break;
		case MEGABYTE:
			unitFactor = 1./(1024.*1024.);
			break;
		case GIGABYTE:
			unitFactor = 1./(1024.*1024.*1024.);
			break;
		default:
			unitFactor = 1.;
			break;
	}
}

long long MemoryMonitor::getSystemMemoryConsumption(){
	sysinfo (&memInfo);
	long long physMemUsed = memInfo.totalram - memInfo.freeram;
	physMemUsed *= memInfo.mem_unit * unitFactor;
	return physMemUsed;
}

long long MemoryMonitor::getProcessMemoryConsumption(){
	sysinfo (&memInfo);
	return getValue()*unitFactor;
}

int MemoryMonitor::parseLine(char* line){
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

long long MemoryMonitor::getValue(){
    FILE* file = fopen("/proc/self/status", "r");
    long long result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result * 1024;
}
