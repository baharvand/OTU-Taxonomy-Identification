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

#include "dictionary.h"

//constructor
template <class Type>
dictionary<Type>::dictionary(){

	cout << "Creating dictionary class" << endl;
	//initialize the variables
	numKeys = 1;
}


//destructor
template <class Type>
dictionary<Type>::~dictionary(){
	
	cout << "Destroying dictionary class" << endl;
}

//return the value associated with a key
//if the key is already in the map, then just return the value
//if not, add the key to the map and return the value
template <class Type>
int dictionary<Type>::getVal(Type key){

	typename map<Type, int>::iterator it;
	
	//search the map for the key
	it = dict.find(key);
	
	//add the key if it is not found
	if (it == dict.end()){
		dict.insert(pair<Type, int>(key, numKeys));
		numKeys++;
		return (numKeys-1);
	}

	//if the key is found return its value
	return it->second;
}

//return the number of elements
template <class Type>
int dictionary<Type>::size(){
	return numKeys-1;
}

template class dictionary<int>;
template class dictionary<string>;
