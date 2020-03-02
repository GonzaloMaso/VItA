/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * ConstantConstraintFunction.h
 *
 *  Created on: 18 de mar de 2018
 *      Author: gonzalo
 */

#ifndef CONSTANTCONSTRAINTFUNCTION_H_
#define CONSTANTCONSTRAINTFUNCTION_H_

#include "AbstractConstraintFunction.h"

template <class T, class S> class ConstantConstraintFunction: public AbstractConstraintFunction<T,S> {
	/**	Constant value. */
	T value;
public:
	ConstantConstraintFunction(T value);
	virtual ~ConstantConstraintFunction();
	T getValue(S treeCondition);
};

template <class T, class S>
ConstantConstraintFunction<T,S>::ConstantConstraintFunction(T value) : AbstractConstraintFunction<T,S>() {
	this->value = value;

}

template <class T, class S>
ConstantConstraintFunction<T,S>::~ConstantConstraintFunction() {
}

template<class T, class S>
T ConstantConstraintFunction<T,S>::getValue(S treeCondition) {
	return value;
}

#endif /* CONSTANTCONSTRAINTFUNCTION_H_ */
