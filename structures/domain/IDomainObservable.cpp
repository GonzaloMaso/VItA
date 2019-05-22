/*
 * IDomainObservable.cpp
 *
 *  Created on: Mar 15, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "IDomainObservable.h"

IDomainObservable::IDomainObservable() {
}

void IDomainObservable::registerObserver(IDomainObserver* observer) {
	observers.push_back(observer);
}

void IDomainObservable::notifyObservers() {
	for(vector<IDomainObserver *>::iterator it = observers.begin(); it != observers.end(); ++it) {
	    (*it)->observableModified(this);
	}
}
