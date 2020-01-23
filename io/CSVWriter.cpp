/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CSVWriter.cpp
 *
 *  Created on: 4/10/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "CSVWriter.h"
#include <iostream>
#include <fstream>

CSVWriter::CSVWriter(string filename) {
	this->filename = filename;
}

CSVWriter::~CSVWriter(){
}

void CSVWriter::addColumn(string header, vector<double> values){
	headers.push_back(header);
	vector<string> converted;
	for (std::vector<double>::iterator it = values.begin(); it != values.end(); ++it) {
		converted.push_back(to_string(*it));
	}
	data.push_back(converted);
}

void CSVWriter::addColumn(string header, vector<int> values){
	headers.push_back(header);
	vector<string> converted;
	for (std::vector<int>::iterator it = values.begin(); it != values.end(); ++it) {
		converted.push_back(to_string(*it));
	}
	data.push_back(converted);
}

void CSVWriter::addColumn(string header, vector<long> values){
	headers.push_back(header);
	vector<string> converted;
	for (std::vector<long>::iterator it = values.begin(); it != values.end(); ++it) {
		converted.push_back(to_string(*it));
	}
	data.push_back(converted);
}

void CSVWriter::addColumn(string header, vector<string> values){
	headers.push_back(header);
	data.push_back(values);
}

void CSVWriter::write(){
	ofstream csvFile;
	csvFile.open(filename);

	//	Writing headers of CSV file
	for (std::vector<string>::iterator it = headers.begin(); it != headers.end(); ++it) {
		if(it == headers.begin()){
			csvFile << *it;
		}
		else{
			csvFile << "," << *it;
		}
	}
	csvFile << "\n";

	//	Writing rows of CSV file
	long long int numRows = data[0].size();
	long long int numCols = headers.size();
	for (int idxRow = 0; idxRow < numRows; ++idxRow) {
		for (int idxCol = 0; idxCol < numCols-1; ++idxCol) {
			csvFile << data[idxCol][idxRow] << ",";
		}
		csvFile << data[numCols-1][idxRow] << "\n";
	}

	csvFile.flush();
	csvFile.close();
}
