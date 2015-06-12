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

#include "correlation.h"

//calculate the correlation and pvalue between 2 vectors
void correlation::corr(vector<float>& a, vector<float>& b, int n, float& correl, float& pval){

	//convert the input vector into a doubel array and store it here
	double* da = new double[n];
	double* db = new double[n];

	//copy the vector into the double arrays
	for(int i=0; i<n; i++){
		da[i] = a[i];
		db[i] = b[i];
	}

	//print the double array
/*	cout << "Double array: " << endl;
	for(int i=0; i<n; i++){
		cout << da[i] << " " << db[i] << endl;
	}
*/	

	
	correl = gsl_stats_correlation(da, 1, db, 1, n);	

	double t = (correl*(sqrt(n-2)) / sqrt(1-(correl*correl)));
	if(t<0){
		pval = 2*(1-gsl_cdf_tdist_P(-t, n-2));
	}
	else {
		pval = 2*(1-gsl_cdf_tdist_P(t, n-2));
	}
	//double pv0=t0<0?2*(1-gsl_cdf_tdist_P(-t0,n-2)):2*(1-gsl_cdf_tdist_P(t0,n-2));

	//cleanup
	delete da;
	delete db;
}
