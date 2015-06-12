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

//The parse arguments class
#pragma once
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "debug.h"

using namespace std;

class parseArgs{

//public members
public:
	parseArgs(int, char**);
	~parseArgs();
	void readArgs(int, char**);
	void checkNextArg(int, int, char**);
	void printHelp();

	//input list file
	char listFile[100];
	//input rdp file
	char rdpFile[100];
	//use the class or phylum from the rdp file
	int rdpPhylum;
	int rdpClass;
	//threashold below which to classify as unknown
	int rdpThresh;
	//use positive or negative correlation
	int posCorr;
	int negCorr;
	//max number of nodes in the graph
	int maxNodes;
	//pval threshold
	int pThresh;

private:
	//check the correctness of the args
	void checkArgs();
};
