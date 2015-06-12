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

#include "mrf.h"

//constructor
mrf::mrf(){

	//initialize the variables
	numUnknowns = 0;
	numKnowns = 0;
	classifyKnown = 0;

	//crate log files
	paramLog.open("parameters.log");
	paramLog << "*************************************************************************" << endl;  
	paramLog << "**                       Classification Parameters                     **" << endl;
	paramLog << "**  Taxonomy <tab> Number of nodes <tab> Alpha <tab> Beta <tab> Gamma  **" << endl;
	paramLog << "*************************************************************************" << endl;
	paramLog << endl;
}

//destructor
mrf::~mrf(){

	//close all log files
	paramLog.close();

}

void mrf::classify(){

	cout << "Starting classifier." << endl;

	//check if there is an unknown node in the taxonon names
	if(checkUnknown()){
		//nothing to do. return
		cout << "Nothing to classify. All nodes are known." << endl;
		return;
	}

	//preprocess the nodes
	//set the prior prob of known nodes, etc.
	preprocessNodes();
	clusterCoeff();
	//print the graph;
	//dfs();

	cout << "Number of unknown nodes: " << numUnknowns << endl; 
	//print to log
	graphLog << "Numver of nodes in the graph: " << numNodes << endl;
    graphLog << "Number of edges in the graph: " << numEdges << endl;
	graphLog << "Number of unknown nodes: " << numUnknowns << endl;
	cout << "Number of taxon labels: " << taxonDict->size() << endl;

	//set the size of the data and response mat
	data.set_size(2,numKnowns);
	responses.set_size(numKnowns);

	//classify every taxonomic label
	for(int i=0; i<taxonDict->size(); i++){
		//skip over the "unknown" 
		if(i != unknownTaxonNum){
			cout << "Classifying taxon: " << taxName[i] << " " << i << endl;
			paramLog << taxName[i] << "\t";
			//first classify the unknown nodes
			classifyUnknownTaxon(i);
			setFinalUnknownProb(i);
		
			//classify the known nodes for verification
			classifyKnown = 1;
			numKnowns--;
			numUnknowns++;
			for(unsigned int j=0; j<nodeList.size(); j++){
        		if(nodeList[j]->known == 1){
					classifyKnownTaxon(i, j);
				}
			}
			numUnknowns--;
			numKnowns++;
			classifyKnown = 0;
	
		}
	}
	
	//print the final graph
	dfs();
	//calculate the clust coeff
	//clusterCoeff();
}


//classify the known nodes for verification
void mrf::classifyKnownTaxon(int taxId, int nodeId){

	//clear the posterior prob
	//set the prior to post	
	for(unsigned int i=0; i<nodeList.size(); i++){
		nodeList[i]->post[taxId] = nodeList[i]->prior[taxId];	
	}

	//mark the node as unknown
	nodeList[nodeId]->known = 0;

	//classify everything again. but no need for parameter estimation
	classifyUnknownTaxon(taxId);

	//set the final prob
	nodeList[nodeId]->finalProb[taxId] = nodeList[nodeId]->finalTmp/(float)10.0;
	//cout << "Final Knowm: " << nodeList[nodeId]->finalProb[taxId] << endl;
	//clear the finalTmp
    for(unsigned int i=0; i<nodeList.size(); i++){
   		nodeList[i]->finalTmp = (float)0.0;
		//cout << "FINALTMP VAL: " << nodeList[i]->finalTmp << endl;
	}
	//make it back as known
	 nodeList[nodeId]->known = 1;
}



//set the final predicted prob for the unknown nodes
//take the avg of the lag-in period predictions and set is as the final prediction
void mrf::setFinalUnknownProb(int taxId){
	for(unsigned int i=0; i<nodeList.size(); i++){
		if(nodeList[i]->known == 0){
			//float tmp = nodeList[i]->finalTmp/(float)10.0;
			//cout << "FINAL PROB: " << taxId << " " << nodeList[i]->finalTmp << " " << tmp << endl;
			nodeList[i]->finalProb[taxId] = nodeList[i]->finalTmp/(float)10.0;
		}
		//clear the finalTmp
		nodeList[i]->finalTmp = (float)0.0;
	}
}


//classify the unknown nodes for taxon id
void mrf::classifyUnknownTaxon(int taxId){

    //clear the posterior prob
    //set the prior to post 
    for(unsigned int i=0; i<nodeList.size(); i++){
        nodeList[i]->post[taxId] = nodeList[i]->prior[taxId];
    }


	//calculate pi
    float pi = calculatePi(taxId);
    //cout << "Pi: " << pi << endl;
	alp = log((pi)/(1-pi));
	if(classifyKnown == 0){
		paramLog << alp << "\t";
	}
	//cout << "Number of nodes: " << pi << endl;


   	//set unknown post prob with probability bernoulli(pi)
	setUnknownPostProb(taxId, pi);		
	
	//needed only when processing unknown nodes
	if(classifyKnown == 0){
  		//estimate parameters
		estimateParameter(taxId);
	}

  	//gibbs sample - train
	//training period - 100
	for(int i=0; i<100; i++){
		//cout << "Gibbs: " << i << endl;
		lagPeriod = 0;
		gibbsSampler(taxId);
	}

  	//gibbs sample - calculate post prob
	//lag-in period - 10
    for(int i=0; i<10; i++){
        //cout << "Gibbs: " << i << endl;
        lagPeriod = 1;
        gibbsSampler(taxId);
    }

}


//gibbs sampler
void mrf::gibbsSampler(int taxId){

	//reset the visited nodes
    resetVisited();

	//traverse the unknown edges and find the posterior probability
	while(1){
        int done = 1;
        for(unsigned int i=0; i<nodeList.size(); i++){
            if(nodeList[i]->visited == 0){
                done = 0;
                traverseGibbsSampler(i, taxId);
            }
        }
        if(done == 1) break;
    }		
}

//traverse procedure for gibbs sampler
void mrf::traverseGibbsSampler(int node, int taxId){
	
	//mark the node as visited
	nodeList[node]->visited = 1;
	
	int m0=0;
	int m1=0;
	//estimate the posterior probability for the unknown nodes
	if(nodeList[node]->known == 0){
		//get m0 and m1 from the adjacent nodes
		for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
			if(nodeList[*i]->post[taxId] >= 0.5){
				m1++;
			}
			else{
				m0++;
			}
		}
	
		//calculate the prob and store it
		float tmp;
		tmp = exp(alp + (bet*m0) + (gam*m1));
		float postProb = tmp/(1 + tmp);		
	
		nodeList[node]->post[taxId] = postProb;

		//keep track of the prob in the lag-in period
		if(lagPeriod == 1){
			nodeList[node]->finalTmp += postProb;		
			//cout << "IN SAMPLER: " << taxId << " " << postProb << " " << nodeList[node]->finalTmp << endl;
		}
	
	}

	
	//recursive call the procedure on all its adjacent nodes
    for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
        if(nodeList[*i]->visited == 0){
            traverseGibbsSampler(*i, taxId);
        }
    }

}



//parameter estimation
void mrf::estimateParameter(int taxId){

	//zero out the mat
	data.zeros();
	responses.zeros();
	rowCount = 0; 

	//reset the visited nodes
	resetVisited();

	//traverse the known nodes and fill the data and responses
	while(1){
        int done = 1;
        for(unsigned int i=0; i<nodeList.size(); i++){
            if(nodeList[i]->visited == 0){
                done = 0;
                traverseEstimateParam(i, taxId);
            }
        }
        if(done == 1) break;
    }

	//data.print();
	//responses.print();
	//logistic regression to get the parameters
	mlpack::regression::LogisticRegression<> lr(data, responses);
	arma::vec parameters = lr.Parameters();
	//set the regression parameters
	//parameters.print();
	//alp = parameters(0);
	bet = parameters(1);
	gam = parameters(2);

	paramLog << bet << "\t" << gam << endl;		
}

//traverse procedure for estimate parameter
void mrf::traverseEstimateParam(int node, int taxId){
	
	nodeList[node]->visited = 1;

	//if it is a known node, then check its surrounding known nodes and get M0 and M1
	int m0=0;
	int m1=0;
	if(nodeList[node]->known == 1){
		for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
		//get the neighbours data	
		if(nodeList[*i]->known == 1){
				if(nodeList[*i]->prior[taxId] == 1){
					m1++;
				}		
				else{
					m0++;
				}	
			}
		}
		data(0, rowCount) = m0;
		data(1, rowCount) = m1;

		//my response
		if(nodeList[node]->prior[taxId] == 1){
        	responses(rowCount) = 1;
    	}
    	else{
        	responses(rowCount) = 0;
    	}

		//increment the count
		rowCount++;
	}

	
	
    //recursive call the procedure on all its adjacent nodes
    for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
        if(nodeList[*i]->visited == 0){
            traverseEstimateParam(*i, taxId);
        }
    }
}


//set unknown post prob with probability bernoulli(pi)
void mrf::setUnknownPostProb(int taxId, float pi){
   
	//reset the counter
    resetVisited();	

	//initialize the PRNG
	const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

	//traverse the unknown nodes and set their post prob with probability bernoulli(pi)
	while(1){
        int done = 1;
        for(unsigned int i=0; i<nodeList.size(); i++){
            if(nodeList[i]->visited == 0){
                done = 0;
                traverseSetUnknownPostProb(i, taxId, pi);
            }
        }
        if(done == 1) break;
    }

	//free the PRNG
	gsl_rng_free (r);

}

//traverse procedure for setUnknownPostProb
void mrf::traverseSetUnknownPostProb(int node, int taxId, float pi){

	nodeList[node]->visited = 1;

	//set the post prob if it is a unknown node
	if(nodeList[node]->known == 0){
		nodeList[node]->post[taxId] = gsl_ran_bernoulli(r, (double)pi);
	}	


    //recursive call the procedure on all its adjacent nodes
    for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
        if(nodeList[*i]->visited == 0){
            traverseSetUnknownPostProb(*i, taxId, pi);
        }
    }

}


//calculate the pi value
float mrf::calculatePi(int taxId){

	//traverse all the nodes and count the number of known nodes with taxId
	//reset the counter
	resetVisited();
	piCount = 0;

	while(1){
        int done = 1;
        for(unsigned int i=0; i<nodeList.size(); i++){
            if(nodeList[i]->visited == 0){
                done = 0;
                traversePi(i, taxId);
            }
        }
        if(done == 1) break;
    }

	if (classifyKnown == 0){
		cout << "Number of nodes: " << piCount << endl;
		paramLog << piCount << "\t";
	}	

	return ((float)piCount/numKnowns);
}

//traverse procedure to calculate pi
void mrf::traversePi(int node, int taxId){

    nodeList[node]->visited = 1;

	//inc count if the node is known and the id match
	if(nodeList[node]->known==1 && nodeList[node]->taxonNum==taxId){
		piCount++;
	}

	//recursive call the procedure on all its adjacent nodes
    for(list<int>::iterator i=adjList[node].begin(); i!=adjList[node].end(); i++){
        if(nodeList[*i]->visited == 0){
            traversePi(*i, taxId);
        }
    }
}



//preprocess the nodes
void mrf::preprocessNodes(){
	cout << "Preprocessing nodes." << endl;
	int taxonSize = taxonDict->size();
	taxName.resize(taxonSize);
	for(unsigned int i=0; i<nodeList.size(); i++){
		

		taxName[nodeList[i]->taxonNum] = nodeList[i]->taxon;

		//initialize the prior and posterior prob
		nodeList[i]->prior.resize(taxonSize, 0);	
		nodeList[i]->post.resize(taxonSize, 0);	
		nodeList[i]->finalProb.resize(taxonSize, 0);

		//set the known variable
		if(nodeList[i]->taxon.compare("unknown") == 0){
			nodeList[i]->known = 0;
			unknownTaxonNum = nodeList[i]->taxonNum;
			numUnknowns++;
		}
		else{
			nodeList[i]->known = 1;
			//set the prior prob in the known nodes
			nodeList[i]->prior.at(nodeList[i]->taxonNum) = 1;
			nodeList[i]->post.at(nodeList[i]->taxonNum) = 1;
			numKnowns++;
		}


	}
}



//ckeck for an unknown node
//return 1 if there is no unknown. else return 0
int mrf::checkUnknown(){
	int taxonSize = taxonDict->size();
	int val = taxonDict->getVal("unknown");
	//cout << "Taxon size: " << taxonSize << " unknown val: " << val << endl;
	if(val > taxonSize){
		return 1;
	} 
	else
	{
		return 0;
	}
}
