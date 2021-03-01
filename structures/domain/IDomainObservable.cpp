/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * IDomainObservable.cpp
 *
 *  Created on: Mar 15, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "IDomainObservable.h"

IDomainObservable::IDomainObservable() {
}

IDomainObservable::~IDomainObservable() {
	(this->observers).clear();
}

void IDomainObservable::registerObserver(IDomainObserver* observer) {
	observers.push_back(observer);
}

void IDomainObservable::notifyObservers() {
	for(vector<IDomainObserver *>::iterator it = observers.begin(); it != observers.end(); ++it) {
	    (*it)->observableModified(this);
	}
}