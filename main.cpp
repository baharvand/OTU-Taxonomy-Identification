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

#include <iostream>
#include "debug.h"
#include "parseArgs.h"
#include "OTUmatrix.h"
#include "mrf.h"

using namespace std;

int main(int argc, char* argv[]) {

#ifdef DEBUG
	cout << "Starting taco!" << endl;
#endif

	//parse the input arguments
	parseArgs* parse = new parseArgs(argc, argv);
	//construct the OTUmatrix
	OTUmatrix* otu = new OTUmatrix();
	//create the OTU matrix from the list file
	otu->createMatrix(parse);
	//create the high-correlated OTU graph
	mrf* otuGraph = new mrf();
	otu->createGraph(parse, otuGraph);
	//classify the unknown taxonomic classes
	otuGraph->classify();

	//cleanup and exit successfully
	delete otuGraph;
	delete otu;
	delete parse;
	return 0;
}
