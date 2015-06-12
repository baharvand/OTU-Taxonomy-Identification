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

//the node of a graph
#pragma once
#include <iostream>
#include <vector>

using namespace std;

class node{

//public members
public:
	node();
	~node();

	string name;
	string taxon;
	int taxonNum;
	int taxonPercent;
	int nodeNum;
	bool visited;
	int known;
	vector<float> prior;
	vector<float> post;
	vector<float> finalProb;
	float finalTmp;

	void printNode();


//private members
private:

};
