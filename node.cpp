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

#include "node.h"

//constructor
node::node(){
	//initialize the variables
	visited = 0;
	known = 0;
	finalTmp = 0.0;
}

//destructor
node::~node(){

}

void node::printNode(){
	cout << "Name: " << name << " Taxon: " << taxon << " Taxon num: " << taxonNum << " Known: " << known << endl;
	cout << "Prior: ";
	for(unsigned int i=0; i<prior.size(); i++){
		cout << prior[i] << " ";
	}
	cout << endl;
	cout << "Post: ";
    for(unsigned int i=0; i<post.size(); i++){
        cout << post[i] << " ";
    }
	cout << endl;
	cout << "Final: ";
    for(unsigned int i=0; i<finalProb.size(); i++){
        cout << finalProb[i] << " ";
    }
    cout << endl;
}

