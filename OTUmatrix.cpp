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

#include "OTUmatrix.h"

//constructor
OTUmatrix::OTUmatrix(){

	cout << "Creating OTUmatrix class" << endl;
	//initialize variables and data structures
	dict = new dictionary<int>();

	numOrigOTU = 0;
	numAddedOTU = 0;
}

//destructor
OTUmatrix::~OTUmatrix(){
	
	cout << "Destroying OTUmatrix class" << endl;
	//destroy the data structures
	delete dict;
}

//create the OTU matrix
//gets the list file as input
void OTUmatrix::createMatrix(parseArgs* parse){

	cout << "Creating OTU matrix from: " << parse->listFile << endl;

	//open and read the input file
	ifstream f1;
	f1.open(parse->listFile, ifstream::in);
	if(f1.fail()){
		cerr << "Cannot open: " << parse->listFile << endl;
	}

	//read every line of file 
	while(1){
	
		//Read the name of the OTu
		string name;
		getline(f1, name, '\t');
		//cout << "OTU name: " << name << endl;	

		//read the rest of the line, the samples
		string line;
		getline(f1, line);
		vector<int> sampleTmp;
		//cout << line << endl;
		//read each sample one by one
		char* tok;
		tok = strtok((char*)line.c_str(), " ,");
		while (tok != NULL){
			//cout << tok << endl;
			//extract the sample number
			int a=0, sampleNum=0;
			char num[15];
			for (unsigned int i=0; i<strlen(tok); i++){
				if (tok[i]>=48 && tok[i]<=57){
					for (unsigned int j=i; j<strlen(tok); j++){
						if (tok[j] == '.') break;
						num[a] = tok[j];
						a++;
					}
					num[a] = '\0';
					break;
				}
			}
			sampleNum = atoi(num);
			//cout << sampleNum << endl;	
			sampleTmp.push_back(sampleNum);		

			tok = strtok(NULL, " ,");
		}

		//cout << "END OTU" << endl;

		numOrigOTU++;
		//if there are more than 5 samples for this OTU, then add the otu to the matrix
		//increment the sample counts. 
		//Discard the OTU if there are less than 5 samples
		if (sampleTmp.size() >= 5){
			addOTU(name, sampleTmp);	
			numAddedOTU++;
		}
	
		//break if eof is reached
		if(f1.eof()){
			break;
		}
	}

	//close the list file
	f1.close();
	
	//resize all the OTUs to a constant size
	resizeMatrix();

	//printMatrix();

	//add the taxonomic nane to each otu from the rdp file
	processRdpFile(parse);

/*	for(unsigned int i=0; i<otuName.size(); i++){
		cout << otuName[i] << " " << taxon[i] << endl;
	}
*/
	//print the otu matrix
	//printMatrix();
	cout << endl;

	//Normalize the otu matrix
	normalizeMatrix();
	//printMatrix();

	//Log file for OTU details
	OTUlog.open("OTU.log");
	OTUlog << "**********************************************************" << endl;
	OTUlog << "**                      OTU Details                     **" << endl;
	OTUlog << "**********************************************************" << endl;
	OTUlog << endl;
	OTUlog << "Number of original OTUs: " << numOrigOTU << endl;
	OTUlog << "Number of OTUs with more than 5 samples: " << numAddedOTU << endl;

	createCorrMatrix();

	OTUlog.close();
}


//process the rdp file
//reads the rdp file and assigns the taxonomic name to the otu 
//if it is above the specified threshold
void OTUmatrix::processRdpFile(parseArgs* parse){
	
	cout << "Reading rdp results from: " << parse->rdpFile << endl;
	//open and read the rdp file
    ifstream f2;
    f2.open(parse->rdpFile, ifstream::in);
    if(f2.fail()){
        cerr << "Cannot open: " << parse->rdpFile << endl;
    }

	string line;
	int count = 0;

	//skip the first 7 lines
	//change this to a generic method later
	for (int i=0; i<7; i++){
		getline(f2, line);
		//cout << line << endl;
	}

	int phy=0;
	if (parse->rdpPhylum == 1){
		cout << "Processing RDP phylum." << endl;
		phy = 1;
	} 
	else {
		cout << "Processing RDP class." << endl;
	}
	cout << "Using RDP threshold as: " << parse->rdpThresh << endl;

	//actual rdp results starts from here
	while(1){

		string rdpName, rdpTaxon, rdpPercent, garble;
		getline(f2, rdpName, ';');
		//cout << "RDP name: " << rdpName << endl;

		//break if eof is reached
        if(f2.eof()){
            break;
        }

		//the next 5 cols in the file are useless, so eat them
		for(int i=0; i<5; i++){
			getline(f2, garble, ';');
		}
		//if we want the phylum, then eat the next two columns
		if(phy == 0){
			getline(f2, garble, ';');
			getline(f2, garble, ';');
		}

		//get the rdp taxon name and its percent
		getline(f2, rdpTaxon, ';');
		getline(f2, rdpPercent, ';');
	
		//eat up the rest of the line
		getline(f2, line);
	
		//set the taxon name to unknown if it is below the cut-off
		int rdpTh = atoi(rdpPercent.c_str());
		if (rdpTh <= parse->rdpThresh){
			rdpTaxon = "unknown";
		}

		//cout << "RDP taxon: " << rdpTaxon << endl;
		//cout << "RDP percent: " << rdpTh << endl;

		//check if the otu and the rdp names match
		//if they match we have a correct entry
		//if not set the taxon as unknown
		int match = 0;
		int numUnknown = 0;
		int tmpCount = count;
		while(1){
			if(rdpName.compare(otuName[tmpCount]) == 0){
				match = 1;
			}
			else{
				numUnknown++;
				//taxon.push_back("unknown");
			}
			tmpCount++;
			if(match == 1){
				for(int i=0; i<numUnknown; i++){
					taxon.push_back("unknown");
					taxonPercent.push_back(-1);
				}
				taxon.push_back(rdpTaxon);
				taxonPercent.push_back(rdpTh);
				count = tmpCount;
				break;
			}
			else if((unsigned)tmpCount == otuName.size()){
				//count++;
				break;
			}
		}

		//make sure the count does not exceed the number of otus
		if((unsigned)count == otuName.size()){
			break;
		}		

		//break if eof is reached
        if(f2.eof()){
            break;
        }
	}

	//push in unknowns in the end to fill the gap
	for(unsigned int i=taxon.size(); i<otuName.size(); i++){
		taxon.push_back("unknown");
		taxonPercent.push_back(-1);
	}
	
	//check if the their sizes match
	if(otuName.size() != taxon.size()){
		cout << "ERROR: OTU name and OTU taxon sizes dont match." << endl;
	}
	
	//close the file
	f2.close();
	
}


//add an element to the OTU matrix
void OTUmatrix::addOTU(string name, vector<int> samples){

	//cout << "OTU name: " << name << endl;
	//store the otu name
	otuName.push_back(name);

	//store the samples
	vector<float> otutmp;
	for(vector<int>::iterator i = samples.begin(); i!=samples.end(); i++){
		int sampleNum;
		sampleNum = dict->getVal(*i);
		//cout << *i << " " << sampleNum <<endl;
		if (otutmp.size() < (unsigned)sampleNum){
			//resize the otutmp to have atleast sampleNum elements
			//intermediate elements are filled with 0
			otutmp.resize(sampleNum, 0);
			//increment the sample count
			otutmp[sampleNum-1]++;
		}	
		else{
			otutmp[sampleNum-1]++;
		}
	}
	
	//add it to the otu matrix
	matrix.push_back(otutmp);	
}


//resize the otu matrix so that every otu has the same length
void OTUmatrix::resizeMatrix(){

	int dictSize = dict->size();
	cout << "Number of keys in dict: " << dict->size() << endl;
	
	for(vector<vector<float> >::iterator i=matrix.begin(); i!=matrix.end(); i++){
		if(i->size() < (unsigned)dictSize){
			i->resize(dictSize, 0);
		}	
	}	
}


//print the otu matrix
void OTUmatrix::printMatrix(){

	cout << "OTU matrix: " << endl;
	int k=0;
	for(vector<vector<float> >::iterator i = matrix.begin(); i!=matrix.end(); i++){
		cout << otuName[k] << ":\t";
		k++;
		for(vector<float>:: iterator j = i->begin(); j!=i->end(); j++){
			cout << *j << " ";
		}
		cout << endl;
	}	
}


//normalize the otu matrix
//find the number of times each sample is present and divide every element
//of that column with that number
void OTUmatrix::normalizeMatrix(){
	
	int numSamples = dict->size();
	for(int i=0; i<numSamples; i++){
		int tot=0;
		//get the total for each column
		for(vector<vector<float> >::iterator j=matrix.begin(); j!=matrix.end(); j++){
			tot += j->at(i); 
		}
		//divide every element in the column with the total
		for(vector<vector<float> >::iterator j=matrix.begin(); j!=matrix.end(); j++){
			j->at(i) /= tot;
		}
	}
	
}


//form the correlation matrix
void OTUmatrix::createCorrMatrix(){

	int size = otuName.size();
	corr.set_size(size, size);
	pvalue.set_size(size, size);
	corr.zeros();
	pvalue.zeros();

	//cout << "Correlation: " << endl;
	int numSamples = dict->size();
	float cor, pval;
	for(int i=0; i<(signed)matrix.size(); i++){
		for(int j=0; j<=i; j++){
			correlation::corr(matrix[i], matrix[j], numSamples, cor, pval);
			corr(i,j) = cor;
			pvalue(i,j) = pval;
			//cout << "(" << cor << ", " << pval << ") ";
		}
		//cout << endl;
	}

	//corr.print();
	//pvalue.print();

}

//create the otu graph
void OTUmatrix::createGraph(parseArgs* parse, mrf* g){
	
	//check which corr to use to create the graph
	int bothCorr=0;
	if(parse->posCorr==1 && parse->negCorr==1){
		bothCorr=1;
	}	


	//get the pvalue threshold
	float pTh;
	pTh = (float)parse->pThresh/100;
	cout << "P-value threshold: " << pTh << endl;

	numPosCorr = 0;
	numNegCorr = 0;
	

	//we dont care about the correlation if it is above the threshold. so set them to 0
	int size = otuName.size();
	for(int i=0; i<size; i++){
		for(int j=0; j<i; j++){
			float absCorr = 0;
			if (corr(i,j) < 0){
				absCorr = (float)(-1.0) * corr(i,j);
			}
			else{
				absCorr = corr(i,j);
			}
			//Conside the correlation only if it is above 0.5 and less than a p-val threshold
			if ((pvalue(i,j) > pTh) || (absCorr < 0.5)){
				corr(i,j) = (float)0.0;
			}
			//count the number od +ve and -ve correlated values
			if(corr(i,j) > 0.0){
				numPosCorr++;
			}
			else if(corr(i,j) < 0.0){
				numNegCorr++;
			}
		}
		corr(i,i) = (float)0.0;
	}

	cout << "Number of OTU's: " << size << endl;
	cout << "Number of +ve correlation less than Pval cutoff: " << numPosCorr << endl;
	cout << "Number of -ve correlation less than Pval cutoff: " << numNegCorr << endl; 

	OTUlog << "Number of +ve correlation(>0.5) less than Pval cutoff: " << numPosCorr << endl;
	OTUlog << "Number of -ve correlation(>0.5) less than Pval cutoff: " << numNegCorr << endl;

	//cout << "After p cutoff: " << endl;
	//corr.print();

	//we we want pos and neg corr, take the absolute value of the corr matrix then take the max 
	if(bothCorr == 1){
		corr = abs(corr);
	}
	if(parse->posCorr == 1){ //take only the max values
		int numAdded = 0;	
		while(1){
			//get the max corr and add an edge between those otus
			uword row, col;
			double maxVal = corr.max(row, col);
			//break if max nodes is reached or if there are no more nodes to add
			if(numAdded == parse->maxNodes || maxVal == 0){
				cout << "Number of edges in the graph: " << numAdded << endl;
				break;
			}		
			g->addNode(otuName[row], otuName[col], taxon[row], taxon[col], taxonPercent[row], taxonPercent[col]);
			//cout << "Adding edge between: " << otuName[row] << " " << row << " and " << otuName[col] << " " << col << " " << numAdded << " "<< g->numNodes << " Corr: " << maxVal << endl;
			
/*			for(int i=0; i<matrix[row].size(); i++){
				cout << matrix[row][i] <<  " ";
			}
			cout << endl;
			for(int i=0; i<matrix[col].size(); i++){
                cout << matrix[col][i] << " ";
            }
			cout << endl;
*/
			corr(row, col) = (float)0.0;
			numAdded++;
		}		
	}
	else if(parse->negCorr == 1){ //take only the min values
        int numAdded = 0;
        while(1){
            //get the max corr and add an edge between those otus
            uword row, col;
            double minVal = corr.min(row, col);
            //break if max nodes is reached or if there are no more nodes to add
            if(numAdded == parse->maxNodes || minVal == 0){
                break;
            }
            g->addNode(otuName[row], otuName[col], taxon[row], taxon[col], taxonPercent[row], taxonPercent[col]);
            //cout << "Adding edge between: " << otuName[row] << " " << row << " and " << otuName[col] << " " << col << endl;
            corr(row, col) = (float)0.0;
            numAdded++;
        } 

	}

	//g->printGraph();
	//g->dfs();
}
