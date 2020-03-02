/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * PiecewiseBifurcationExponentEstimator.h
 *
 *  Created on: 18 de mar de 2018
 *      Author: gonzalo
 */

#ifndef CONSTANTPIECEWISECONSTRAINTFUNCTION_H_
#define CONSTANTPIECEWISECONSTRAINTFUNCTION_H_

#include "AbstractConstraintFunction.h"

#include <vector>

using namespace std;

/**
 * Constraint piecewise function with different constant values for each condition.
 */
template <class T, class S> class ConstantPiecewiseConstraintFunction: public AbstractConstraintFunction<T,S> {
	/**	Piecewise constant values for each function piece. */
	vector<T> values;
	/**	Lower bound condition for each function piece. The upper bounds are the next lower bound. S must define the operator >=. */
	vector<S> conditions;
public:
	/**
	 * Constructor.
	 * @param values Piecewise constant values for each function piece.
	 * @param conditions Lower bound condition for each function piece.
	 */
	ConstantPiecewiseConstraintFunction(vector<T> values, vector<S> conditions);
	/**
	 * Destructor.
	 */
	virtual ~ConstantPiecewiseConstraintFunction();
	/**
	 * Estimator based on the condition @p treeCondition.
	 * @param treeCondition Condition value.
	 * @return Constraint function value.
	 */
	T getValue(S treeCondition);
};

template <class T, class S>
ConstantPiecewiseConstraintFunction<T,S>::ConstantPiecewiseConstraintFunction(vector<T> values, vector<S> conditions) : AbstractConstraintFunction<T,S>(){
	this->values = values;
	this->conditions = conditions;
}

template <class T, class S>
ConstantPiecewiseConstraintFunction<T,S>::~ConstantPiecewiseConstraintFunction() {
}

template <class T, class S>
T ConstantPiecewiseConstraintFunction<T,S>::getValue(S treeCondition){
	int i = 0;
	for(vector<int>::iterator it = conditions.begin(); it != conditions.end(); ++it ){
		if(treeCondition >= *it)
			++i;
		else
			break;
	}
	--i;
	return values[i];
}

#endif /* CONSTANTPIECEWISECONSTRAINTFUNCTION_H_ */
