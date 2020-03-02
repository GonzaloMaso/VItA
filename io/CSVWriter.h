/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CSVWriter.h
 *
 *  Created on: 4/10/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef CSVWRITER_H_
#define CSVWRITER_H_
#include <string>
#include <vector>

using namespace std;

/**
 * Writes a set of vectors of the same size as columns in a CVS file.
 */
class CSVWriter {
	/* Output CSV filename */
	string filename;
	/* Column headers */
	vector<string> headers;
	/* Data organised as column first */
	vector<vector<string>> data;
public:
	/**
	 * Constructor.
	 * @param filename Output CSV filename.
	 */
	CSVWriter(string filename);
	/**
	 * Destructor.
	 */
	virtual ~CSVWriter();

	/**
	 * Appends a new column to the current structure.
	 * @param header Header of the new column.
	 * @param values Data for the new column.
	 */
	void addColumn(string header, vector<double> values);
	/**
	 * Appends a new column to the current structure.
	 * @param header Header of the new column.
	 * @param values Data for the new column.
	 */
	void addColumn(string header, vector<int> values);
	/**
	 * Appends a new column to the current structure.
	 * @param header Header of the new column.
	 * @param values Data for the new column.
	 */
	void addColumn(string header, vector<long> values);
	/**
	 * Appends a new column to the current structure.
	 * @param header Header of the new column.
	 * @param values Data for the new column.
	 */
	void addColumn(string header, vector<string> values);

	/**
	 * Writes the CSV file. Execute this only after finalising of adding all columns.
	 */
	void write();
};

#endif /* CSVWRITER_H_ */
