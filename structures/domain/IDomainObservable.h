/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * IDomainObservable.h
 *
 *  Created on: Mar 15, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_IDOMAINOBSERVABLE_H_
#define DOMAIN_IDOMAINOBSERVABLE_H_

#include "IDomainObserver.h"
#include <vector>

using namespace std;

/**
 * Observable class coplaint to observer/observable pattern.
 */
class IDomainObservable {
	/**	List of observer objects. */
	vector<IDomainObserver *> observers;
public:
	/**
	 * Constructor.
	 */
	IDomainObservable();
	/**
	 * Destructor
	 */
	~IDomainObservable();
	/**
	 * Registrate a new observer in order to notify each time that the observable class is updated.
	 * @param observer Observer class type IDomainObserver.
	 */
	void registerObserver(IDomainObserver *observer);
	/**
	 * Notify to all registrated observer executing the notified.
	 */
	void notifyObservers();
};

#endif /* DOMAIN_IDOMAINOBSERVABLE_H_ */
