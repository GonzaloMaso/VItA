/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractSavingTask.h
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef ABSTRACTSAVINGTASK_H_
#define ABSTRACTSAVINGTASK_H_
#include "../../structures/tree/AbstractObjectCCOTree.h"

/**
 * Abstract class part of an executer design pattern. Each time a tree is saved a set of executers will write the tree content.
 * Different executers provide specific persistance of the tree, such as the generation of CCO files or VTK, depending if is a checkpoint
 * for further tree generation or just a visualisation temporal file.
 */
class AbstractSavingTask {
public:
	AbstractSavingTask();
	virtual ~AbstractSavingTask();

	virtual void execute(long long int terminals, AbstractObjectCCOTree *tree) = 0;
};

#endif /* ABSTRACTSAVINGTASK_H_ */
