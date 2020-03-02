/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractConstraintFunction.h
 *
 *  Created on: 18 de mar de 2018
 *      Author: gonzalo
 */

#ifndef ABSTRACTCONSTRAINTFUNCTION_H_
#define ABSTRACTCONSTRAINTFUNCTION_H_

/**
 * Abstract class to model tree constraints with different data type
 */
template <class T, class S> class AbstractConstraintFunction {
public:
	AbstractConstraintFunction();
	virtual ~AbstractConstraintFunction();
	virtual T getValue(S treeCondition) = 0;
};

template <class T, class S>
AbstractConstraintFunction<T,S>::AbstractConstraintFunction() {
}

template <class T, class S>
AbstractConstraintFunction<T,S>::~AbstractConstraintFunction() {
}

#endif /* ABSTRACTCONSTRAINTFUNCTION_H_ */
