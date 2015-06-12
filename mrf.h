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

//the markov randon field calss
//exdends the graph class
#pragma once
#include <iostream>
//#include <armadillo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "graph.h"
#include <math.h>
#include <mlpack/core.hpp>
#include <libxml/parser.h>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>

using namespace std;

class mrf: public graph{

//public members
public:
	mrf();
	~mrf();	

	//the method where all the classification happens
	void classify();

//private members
private:
	int unknownTaxonNum;
	int numKnowns;
	int numUnknowns; 

	//keep track of the taxon names
	//moving this to graph.h to make log files easier
	//vector<string> taxName;

	//check for an unknown node
	int checkUnknown();

	//preprocess the nodes
	void preprocessNodes();

	//classify the unknown nodes with taxon id
	void classifyUnknownTaxon(int);

	//calculate the pi value
	float calculatePi(int);
	void traversePi(int, int);
	int piCount;

	//set unknown post prob with probability bernoulli(pi)
    void setUnknownPostProb(int, float);		
	void traverseSetUnknownPostProb(int, int, float);
 	gsl_rng * r;

	//parameter estimation from the known part of the network using logistic regression
	void estimateParameter(int);
	void traverseEstimateParam(int, int);
	arma::mat data; // The dataset itself.
	arma::vec responses; // The responses, one row for each row in data.
	int rowCount;
	double alp, bet, gam;

	//gibbs sampler
	void gibbsSampler(int);
	void traverseGibbsSampler(int, int);
	int lagPeriod;

	//set the final predicted prob for the unknown nodes
	void setFinalUnknownProb(int);

	//classify the known nodes for verification
   	void classifyKnownTaxon(int, int);
	int classifyKnown;

	//log files
	ofstream paramLog;
};
