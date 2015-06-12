Taxonomic identification of unknown OTUs through OTU co-abundance networks
========================================================


Library Dependencies
====================
  GSL - GNU Scientific Library  
  MLPACK - open-source scalable c++ machine learning library  

All of those should be available in your distribution's package  
manager.  If not, you will have to compile each of them by hand.   
See the documentation for each of those packages for more information.  


Installation Instructions
=========================
You can then build it by typing,  

  make  

If you encounter problems, remove any existing compiled  
files with,  

  make clean  

before running "make" again  


Running
============
In the MRF directory,  
  ./ [arguments]  

Arguments:  
--list      A file containing list of clusters, output of the CROP clustering program                                         
--rdp       A file with RDP taxonomy for each cluster  
--rdpThresh 	  The threshold below which a cluster will be considered unknown                                                 
--rdpClass 		  Classify using RDP class level result                                                                          
--rdpPhylum 	  Classify using RDP phylum level result                                                                         
--posCorr 		  Use positive correlated OTU's                                                                                 
--negCorr 		  Use Negative correlated OTU's                                                                                 --pThresh 		  The maximum p-value for correlation between OTU's                                                             
--maxEdges 		  The maximum number of edges in the OTU network                                                                 
