/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * IObserver.h
 *
 *  Created on: Mar 15, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_IDOMAINOBSERVER_H_
#define DOMAIN_IDOMAINOBSERVER_H_

class IDomainObservable;

/**
 * Class that implements the observer end of an observer/observable pattern.
 * Observable object will call observableModified when its state changes. An observer
 * object must register in the observable object in order to be called.
 */
class IDomainObserver {
public:
	/**
	 * Constructor
	 */
	IDomainObserver();
	/**
	 * Destructor
	 */
	virtual ~IDomainObserver() = 0;
	/**
	 * Method executed when @p modifiedDomain has changed.
	 * @param observableInstance Observable domain.
	 */
	virtual void observableModified(IDomainObservable * observableInstance) = 0;
};

#endif /* DOMAIN_IDOMAINOBSERVER_H_ */
