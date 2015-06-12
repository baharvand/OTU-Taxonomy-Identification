/*  
 *    Copyright (C) 2015 Zohreh Baharvand Irannia, Rishvanth Prabakar
 *
 *    Authors: Rishvanth Prabakar
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//The OTU matrix class
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <mlpack/core.hpp>
#include "parseArgs.h"
#include "dictionary.h"
#include "correlation.h"
#include "mrf.h"
#include "debug.h"

using namespace std;
using namespace arma;

class OTUmatrix{

//public members
public:
	OTUmatrix();
	~OTUmatrix();
	void createMatrix(parseArgs*);
	//create the otu graph with the high correlated otus
	void createGraph(parseArgs*, mrf*);

//private members
private:
	//data structure to store the OTU matrix
	//columns are OTUs, rows are samples
	vector<vector<float> > matrix;
	
	//the correlation amd pval matrix
	Mat<float> corr;
	Mat<float> pvalue;
	int numPosCorr;
	int numNegCorr;

	int numOrigOTU;
	int numAddedOTU;

	//log files
	ofstream OTUlog;	

	//store the name of each OTU
	vector<string> otuName; 

	//store the taxonomic name for each OTU
	vector<string> taxon;
	//store the pobability of the classification
	vector<int> taxonPercent;

	//dictionary to map the sample numbers to a continuous int value
	dictionary<int>* dict;

	//process the rdp file
	void processRdpFile(parseArgs*);

	//add an element to the OTU matrix
	void addOTU(string, vector<int>);	

	//resize the otu matrix so that every otu has the same length
	void resizeMatrix();
	
	//print the otu matrix
	void printMatrix();

	//normalize the otu matrix
	void normalizeMatrix();	

	//create the correlation matrix beteen every otu pair
	void createCorrMatrix();

};
