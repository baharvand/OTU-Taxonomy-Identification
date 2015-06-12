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

//the graph class
#pragma once
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include "dictionary.h"
#include "node.h"

using namespace std;

class graph{

//public members
public: 
	graph();
	~graph();
	//add an edge between 2 nodes
	void addNode(string, string, string, string, int, int);
	//print the graph
	void printGraph();
	//walk through the graph
	void dfs();

	//find the cluster coefficients
	void clusterCoeff();
	
	int numNodes;
	int numEdges;

//protected members
protected:
	//adjacency list for the grpah
	vector<list<int> > adjList;
	//store the nodes in this vectoe
	vector<node*> nodeList;
	//doctionary to map the input otu names
	dictionary<string>* otuDict;
	//dictionary to map the taxon names
	dictionary<string>*	taxonDict;
	//traverse the graph. the actual dfs procedure
	void traverse(int);
	//traverse provedure for cluster coeff
	void clusterCoeffTraverse(int); 
	//reset the visited values
	void resetVisited();

	//used in mrf. declaring this here to log easier
	vector<string> taxName;

	//to keep track of the distribution of nodes in the subgraph
	vector<int> taxonDistribution;
	int subGraphNumNodes;	

	//to calculate the cluster coefficients
	vector<float> taxonCoeff;
	vector<int> taxonCount;

	//log files
	ofstream nodeLog;
	ofstream edgeLogTaxon;
	ofstream edgeLogOTU;
	ofstream graphLog;
	ofstream unknownPostLog;
	ofstream knownPostLog; 
	ofstream ukMostProb;
	ofstream kMostProb;
	ofstream clustCoeff;
};
