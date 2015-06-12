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

#include "graph.h"

//constructor
graph::graph(){

	//create and initialize the needed structures
	numNodes = 0;
	numEdges = 0;
	otuDict = new dictionary<string>();
	taxonDict = new dictionary<string>();

	//create log files
	nodeLog.open("nodes.log");
	nodeLog << "**********************************************************" << endl;
	nodeLog << "**                 OTU Network Nodes                    **" << endl;
	nodeLog << "**    Node Number <tab> Node OTU <tab> Node Taxonomy    **" << endl;
	nodeLog << "**********************************************************" << endl;
	nodeLog << endl;

	edgeLogTaxon.open("edgesTaxon.log");
	edgeLogTaxon << "**********************************************************" << endl;
	edgeLogTaxon << "**      			OTU Network Edges				     **" << endl;
	edgeLogTaxon << "**		  Node Taxonomy <-------> Node Taxonomy  		 **" << endl;
	edgeLogTaxon << "**********************************************************" << endl;
	edgeLogTaxon << endl;

	edgeLogOTU.open("edgesOTU.log");
	edgeLogOTU << "**********************************************************" << endl;
	edgeLogOTU << "**                 OTU Network Edges                    **" << endl;
	edgeLogOTU << "**       	 Node Num <-------> Node Num  	           **" << endl;
	edgeLogOTU << "**********************************************************" << endl;
	edgeLogOTU << endl;

	graphLog.open("graph.log");
	graphLog << "**********************************************************" << endl;
	graphLog << "**                 OTU Network Details					 **" << endl;
	graphLog << "**********************************************************" << endl;
	graphLog << endl;

	unknownPostLog.open("unknownPosterior.log");
	unknownPostLog << "**************************************************************" << endl;
	unknownPostLog << "**	  	 Posterior Probabilities for unknown nodes 		   **" << endl;
	unknownPostLog << "** Node Nunber <tab> Node OTU <tab> [taxonomy: probability] **" << endl;
	unknownPostLog << "**************************************************************" << endl;
	unknownPostLog << endl;

	


}

//destructor
graph::~graph(){
	
	//walk through the node list and destroy each node
	for(unsigned int i=0; i<nodeList.size(); i++){
		delete nodeList[i];
	}	

	//close the log files
	nodeLog.close();
	edgeLogTaxon.close();
	edgeLogOTU.close();
	graphLog.close();
	unknownPostLog.close();

	//destroy ant created structures
	delete otuDict;
	delete taxonDict;
}

//add an edge between two nodes
void graph::addNode(string a, string b, string ta, string tb, int tap, int tbp){

	//get the number of the node from the dictionary
	int aVal = otuDict->getVal(a);
	int bVal = otuDict->getVal(b);
	aVal--;
	bVal--;

	//add the taxon names to the dict. 
	int taVal = taxonDict->getVal(ta);
	int tbVal = taxonDict->getVal(tb);
    taVal--;
    tbVal--;

	//check if the node already exists, it not create it and make the adjList bigger 
	if(numNodes <= aVal){
		node* newNode = new node();
		newNode->name = a;
		newNode->taxon = ta;
		newNode->taxonNum = taVal;
		newNode->taxonPercent = tap;	
		newNode->nodeNum = numNodes;
		nodeList.push_back(newNode);
		adjList.resize(aVal+1);
	
		//log data
		nodeLog << numNodes << "\t" << a << "\t" << ta << endl;
 
		numNodes++;
	}
	if(numNodes <= bVal){
        node* newNode = new node();
		newNode->name = b;
		newNode->taxon = tb;
		newNode->taxonNum = tbVal;
		newNode->nodeNum = numNodes;
        newNode->taxonPercent = tbp;
		nodeList.push_back(newNode);
		adjList.resize(bVal+1);

		//log data
		nodeLog << numNodes << "\t" << b << "\t" << tb << endl;
        numNodes++;
    }

	//increment the number of edges
	numEdges++;

	//add the nodes to the adjacency list
	adjList.at(aVal).push_back(bVal);
	adjList.at(bVal).push_back(aVal);

	//log data
	edgeLogTaxon << ta << "\t" << tb << endl;
	edgeLogOTU << aVal << "\t" << bVal << endl;

	//cout << "A value: " << aVal << endl;
	//cout << "B value: " << bVal << endl;
}

//print the graph
void graph::printGraph(){

	cout << "Graph: " << endl;
	for(unsigned int i=0; i<adjList.size(); i++){
		cout << i << ": ";
		for(list<int>::iterator j=adjList[i].begin(); j!=adjList[i].end(); j++){
			cout << *j << " ";
		}
		cout << endl;
	}
}

//depth first search
//just sets up the dfs and calls traverse on the first node
void graph::dfs(){
	
	//reset the visited variable in every node
	for(unsigned int i=0; i<nodeList.size(); i++){
		nodeList[i]->visited = 0;
	}

	//nothing to traverse
	if(nodeList.size() == 0){
		return;
	}

int numSubgraph = 0;
	while(1){
		int done = 1;
		for(unsigned int i=0; i<nodeList.size(); i++){
			if(nodeList[i]->visited == 0){
				numSubgraph++;
				done = 0;
				subGraphNumNodes = 0;
				taxonDistribution.clear();
				taxonDistribution.resize(taxonDict->size());
				traverse(i);
				//cout << "Subgraph at node: " << i << endl; 	
				cout << endl;
				//print to log
				graphLog << endl << "Number of nodes in subgraph " << numSubgraph << " :" << subGraphNumNodes << endl;
				for(unsigned int k=0; k<taxonDistribution.size(); k++){
					graphLog << taxName[k] << " :" << taxonDistribution[k] << endl;
				}
			}	
		}
		if(done == 1) break;
	}

	cout << "Numver of nodes in the graph: " << numNodes << endl;
	cout << "Number of subgraphs: " << numSubgraph << endl;

	//print to log file
	//graphLog << "Numver of nodes in the graph: " << numNodes << endl;
	//graphLog << "Number of edges in the graph: " << numEdges << endl;  
    graphLog << endl << "Number of subgraphs: " << numSubgraph << endl;
}


//actual dfs happens here
void graph::traverse(int node){

	nodeList[node]->visited = 1;
	//cout << "Visited: " << nodeList[node]->name << endl;

	//increment the number of nodes in the subgraph
	subGraphNumNodes++;
	//keep track of the taxon dist
	taxonDistribution[nodeList[node]->taxonNum]++;	

	//print to posterior probability log files
	if(nodeList[node]->known==0){
		unknownPostLog << nodeList[node]->nodeNum << "\t" << nodeList[node]->name << "\t";
		for(unsigned int i=0; i<nodeList[node]->finalProb.size(); i++){
        	unknownPostLog << taxName[i] << ": " << nodeList[node]->finalProb[i] << "\t";
    	}	
		unknownPostLog << endl;
	}

    if(nodeList[node]->known==1){
        knownPostLog << nodeList[node]->nodeNum << "\t" << nodeList[node]->name << "\t" << nodeList[node]->taxon << "\t";
        for(unsigned int i=0; i<nodeList[node]->finalProb.size(); i++){
            knownPostLog << taxName[i] << ": " << nodeList[node]->finalProb[i] << "\t";
        }
        knownPostLog << endl;
    }	


	//print to most probable log file
	if(nodeList[node]->known==0){
		ukMostProb << nodeList[node]->nodeNum << "\t" << nodeList[node]->name << "\t";
		float max=0.0;
		int maxi=0;
		for(unsigned int i=0; i<nodeList[node]->finalProb.size(); i++){
			if(nodeList[node]->finalProb[i] > max){
				max = nodeList[node]->finalProb[i];
				maxi = i;
			}
		}
		ukMostProb << taxName[maxi] << "\t" << max;
		ukMostProb << endl;
	}

	if(nodeList[node]->known==1){
		kMostProb << nodeList[node]->nodeNum << "\t" << nodeList[node]->name << "\t";
        float max=0.0;
        int maxi=0;
        for(unsigned int i=0; i<nodeList[node]->finalProb.size(); i++){
            if(nodeList[node]->finalProb[i] > max){
                max = nodeList[node]->finalProb[i];
                maxi = i;
            }
        }
		kMostProb << nodeList[node]->taxon << "\t" << nodeList[node]->taxonPercent << "\t";
        kMostProb << taxName[maxi] << "\t" << max;
        kMostProb << endl;
	}


/*	if (nodeList[node]->known==0){
		nodeList[node]->printNode();
	}
*/
	for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
		if(nodeList[*i]->visited == 0){
			traverse(*i);
		}
	}	
}


//find the cluster coefficients for each node
void graph::clusterCoeff(){

	//reset the visited
	resetVisited();

	//initialize the data structures
	taxonCoeff.resize(taxonDict->size());	
	taxonCount.resize(taxonDict->size());


	while(1){
        int done = 1;
        for(unsigned int i=0; i<nodeList.size(); i++){
            if(nodeList[i]->visited == 0){
                done = 0;
                clusterCoeffTraverse(i);
            }
        }
        if(done == 1) break;
    }

	//calculate the avg coeff for every taxon and print to file
	for(unsigned int i=0; i<taxonCoeff.size(); i++){
		float clusterCoeff = taxonCoeff[i] / taxonCount[i];
		clustCoeff << taxName[i] << "\t" << clusterCoeff << endl;	
	}	

}


//traverse procedure for finding cluster coefficients
void graph::clusterCoeffTraverse(int node){

	
	int numEdges = 0;
	int numNeigh = 0;
	vector<int> neigh;
	//get all the neigh nodes
	for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
		neigh.push_back(nodeList[*i]->nodeNum);
	}	
	numNeigh = neigh.size();

	//get the number of neigh nodes for the neigh nodes
	for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
		for(list<int>::iterator j=adjList[*i].begin(); j!=adjList[*i].end(); j++){
			int nodeNumber = nodeList[*j]->nodeNum;
			for(unsigned int k=0; k<neigh.size(); k++){
				if(nodeNumber == neigh[k]){
					numEdges++;
				}
			}
		}
	}

	//divide by 2 because its a undirected graph
	numEdges /= 2;

	//calculate the clustering coeff	
	float tmp;
	tmp = (numNeigh*(numNeigh-1))/2;
	float coeff = numEdges/tmp;

	taxonCoeff[nodeList[node]->taxonNum] += coeff;
	taxonCount[nodeList[node]->taxonNum]++;
	

	//cout << nodeList[node]->nodeNum << "\t" << numNeigh << "\t" << numEdges << "\t" << coeff << endl;
	//print to the log file
	//clustCoeff << nodeList[node]->nodeNum << "\t" << coeff << endl;

	nodeList[node]->visited = 1;
    for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
        if(nodeList[*i]->visited == 0){
            clusterCoeffTraverse(*i);
        }
    }

}

//reset the visited values
void graph::resetVisited(){
	
	//reset the visited variable in every node
    for(unsigned int i=0; i<nodeList.size(); i++){
        nodeList[i]->visited = 0;
    }
}

