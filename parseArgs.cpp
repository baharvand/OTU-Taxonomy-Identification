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

#include "parseArgs.h"

//constructor
parseArgs::parseArgs(int argc, char** argv){
#ifdef DEBUG
	cout << "Crating parseArgs class" << endl;
	cout << "Number of arguments: " << argc-1 << endl;
#endif

	//initialize the default values
	rdpPhylum = 0;
	rdpClass = 0;
	rdpThresh = 5;
	posCorr = 0;
	negCorr = 0;
	maxNodes = 10000000;
	pThresh = 5;

	//parse the args
	readArgs(argc, argv);
}

//destructor
parseArgs::~parseArgs(){
#ifdef DEBUG
	cout << "Destroying parseArgs class" << endl;
#endif
}

//parse the input args and set the variables
void parseArgs::readArgs(int argc, char** argv){

	//set the default values

	//read the input args one by one
	int i=1;
	while(i<argc){
		//cout << argv[i] << endl;
		if(strcmp(argv[i], "--help") == 0){
            printHelp();
            i++;
			exit(1);
        }
		else if(strcmp(argv[i], "--list") == 0){
			//check if the filename is given
			checkNextArg(i, argc, argv);	
			strcpy(listFile, argv[i+1]);
			cout << "Input file 1: " << listFile << endl;
			i=i+2;
		}
		else if(strcmp(argv[i], "--rdp") == 0){
			checkNextArg(i, argc, argv);
            strcpy(rdpFile, argv[i+1]);
            cout << "RDP file: " << rdpFile << endl;
			i=i+2;
        }
		else if(strcmp(argv[i], "--RDPthresh") == 0){
			checkNextArg(i, argc, argv);
			rdpThresh = atoi(argv[i+1]);
			i += 2;
		}
		else if(strcmp(argv[i], "--RDPclass") == 0){
			rdpClass = 1;
			i++;
		}
		else if(strcmp(argv[i], "--RDPphylum") == 0){
            rdpPhylum = 1;
            i++;
        }
		else if(strcmp(argv[i], "--PosCorr") == 0){
            posCorr = 1;
            i++;
        }
		else if(strcmp(argv[i], "--NegCorr") == 0){
            negCorr = 1;
            i++;
        }
		else if(strcmp(argv[i], "--MaxEdges") == 0){
            checkNextArg(i, argc, argv);
            maxNodes = atoi(argv[i+1]);
            i += 2;
        }
		else if(strcmp(argv[i], "--pThresh") == 0){
            checkNextArg(i, argc, argv);
            pThresh = atoi(argv[i+1]);
            i += 2;
        }

		else{ //unknown option
			cerr << "Unrecognized option." << endl;
			i++;
			//print help and exit
			printHelp();
			exit(1);
		}
		
	}

	//check the input args
	checkArgs();
}


//check the correctness of the input args
void parseArgs::checkArgs(){
	
	if(strlen(listFile) == 0){
        cerr << "The list file is not given as input." << endl;
        //print help and exit
        printHelp();
        exit(1);
    }

	if(strlen(rdpFile) == 0){
        cerr << "The RDP file is not given as input." << endl;
        //print help and exit
        printHelp();
        exit(1);
    }

	if(rdpPhylum == 0 && rdpClass == 0){
		cout << "RDP class/phylum not set. Using default as phylum" << endl;
		rdpPhylum = 1;
	}
	else if(rdpPhylum == 1 && rdpClass == 1) {
		cout << "Both RDP class and phylum are set. Using default as only phylum" << endl;
		rdpClass = 0;
		rdpPhylum = 1;
	}

	if(posCorr==0 && negCorr==0){
		cout << "+/- correlation is not set. Using default as + correlation" << endl;
		posCorr = 1;
	}

}


//check that current argument is not the last argument
// and if the next argument is valid
void parseArgs::checkNextArg(int i, int argc, char** argv){
	//there is no more parameters
	if (i+1 == argc){
		//print help and exit
		cerr << "Invalid argument" << endl;
        printHelp();
        exit(1);
	}

	//next arg starts with '-'
	if (argv[i+1][1] =='-'){
		//print help and exit
		cerr << "No input file." << endl;
        printHelp();
        exit(1);
	}
}


//print TACO help
void parseArgs::printHelp(){

	cout << "Usage: taco [arguments]" << endl;
	cout << endl;
	cout << "Arguments: " << endl;
	cout << "--list \t\t\t A file containing list of clusters" << endl;
	cout << "--rdp \t\t\t A file with RDP taxonomy for each cluster" << endl;
	cout << "--rdpThresh \t\t The threshold below which a cluster will be considered unknown"\
																<< endl; 
	cout << "--rdpClass \t\t Classify using RDP class" << endl;
	cout << "--rdpPhylum \t\t Classify using RDP phylum" << endl; 
	cout << "--posCorr \t\t Use positive correlated OTU's" << endl;
	cout << "--negCorr \t\t Use Negative correlated OTU's" << endl;
	cout << "--pThresh \t\t The maximum p-value for correlation between OTU's" << endl;
	cout << "--maxEdges \t\t The maximum number of edges in the OTU network" << endl;
}
